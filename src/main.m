clc
clear
close all

%% add EPANET matlab toolkit function
addpath(genpath( '../' ))
disp('All files and folder Paths Loaded.');

%% Read Premise Plumbing Features
% load premise plumbing .INP file and save it in a PP variable
PP = epanet( 'HouseholdNDandIRnode.inp');
% Read all pattern values and Indices
k = 0;
PP.setNodeDemandPatternIndex(1,0); % remove the first pattern

% get node properties
for i  = PP.getNodeDemandPatternIndex{:}
    k = k + 1;
    if i == 0
        patternID{ k } = [];
        patternValue { k } = [];
    else
        patternID{ k } = PP.PatternNameID{i};
        patternValue{ k } = PP.Pattern(i,:);
    end
end
% nodes information
pp_node = struct(...
    'ID',PP.NodeNameID',...
    'X', num2cell( PP.NodeCoordinates{1} + (0 - min(PP.NodeCoordinates{1}) ) ),... % shift it to origin
    'Y', num2cell( PP.NodeCoordinates{2} + (0 - min(PP.NodeCoordinates{2}) ) ),... % shift it to origin
    'BDemand', num2cell(PP.NodeBaseDemands{:}'),...
    'Pattern_ind', num2cell( PP.getNodeDemandPatternIndex{:}' ),...
    'PatternID', patternID',...
    'Pattern', patternValue' ,...
    'Elv', num2cell(PP.NodeElevations'),...
    'Type', num2cell(PP.NodeTypeIndex') );


% get pipes information
pp_pipe = struct(...
    'ID', PP.LinkNameID',...
    'Xver', PP.NodeCoordinates{3},...
    'Yver', PP.NodeCoordinates{4},...
    'i', PP.NodesConnectingLinksID(:,1),...
    'j', PP.NodesConnectingLinksID(:,2),...
    'L', num2cell(PP.LinkLength') ,...
    'D', num2cell(PP.LinkDiameter'),...
    'C', num2cell(PP.LinkRoughnessCoeff'),...
    'loss', num2cell(PP.LinkMinorLossCoeff'),...
    'state', num2cell(PP.LinkInitialStatus'),...
    'type', num2cell( PP.LinkTypeIndex' ),...
    'bulk', num2cell(PP.LinkBulkReactionCoeff'),...
    'wall', num2cell(PP.LinkWallReactionCoeff') );

% PP.saveInputFile(fullfile('../data/network/','temp_PP.inp')); % save the premis plumbing
PP.unload; %unload the premise plumbing file


%% Read Network file
%  Read network .INP file
WDN = epanet(fullfile('../data/network','WDN_2019_zeroDemandIntermediateNode.inp'));
% get ALL nodes base demands
nodes_basedemand = WDN.getNodeBaseDemands{1};
% Find all terminal Nodes = non-zero basedemand
target_node_ind = WDN.getNodeIndex(...
    WDN.getNodeNameID( find( nodes_basedemand > 0 ) ) );
% set node's demand pattern to zero (remove all nodes' pattern)
WDN.setNodeDemandPatternIndex(WDN.getNodeIndex,zeros(1,numel(nodes_basedemand) ));
% remove all nodes' base demand
WDN.setNodeBaseDemands(WDN.getNodeIndex,zeros(1,numel(nodes_basedemand) ));
% Change demand pattern of IR and ND to WDN pattern
sample_pattern = WDN.Pattern( WDN.getPatternIndex('MaxDayDemand') ,:);
pattern_temp = num2cell( repmat(sample_pattern(1:24),[2,1]),2);
[pp_node( find( [pp_node.Pattern_ind] ) ).Pattern] = pattern_temp{:};


%% create households porperties
total_houses = sum(nodes_basedemand)/.2;
num_hh = ceil( total_houses/numel(target_node_ind) )-1; % number of houses for at each TN
num_hh = 18 % manually modify the number of houses
angles = conv( linspace( 0, 360, num_hh+1 ), [ 0.5 0.5 ], 'valid' );


scale_factor = .02; % resize the PP to fit in the WDN
demand = {0.2;0};
Ndemand = {-0.2;0};
IRonly = true;
IRND = false;
for tn = target_node_ind
    
    for theta = angles
        if IRonly
            [pp_node( [2,3] ).BDemand] = demand{:};
        elseif IRND
            if rand(1) <= 0.5
                [pp_node( [2,3] ).BDemand] = demand{:};
                
            else
                [pp_node( [3,2] ).BDemand] = Ndemand{:};
            end
        end
        
        addHousehold( WDN, tn, pp_node, pp_pipe, scale_factor, theta )
    end
    
end
WDN.saveInputFile(fullfile('../data/network/',['sample_IRonly', num2str(num_hh),'houses.inp']))
WDN.unload



%%
wdn = epanet(fullfile('../data/network','sample_IRonly18houses_Modified Pump.inp'));
nodes_id = wdn.getNodeNameID;
nodes_index = wdn.getNodeIndex;
irs_id = nodes_id(startsWith(nodes_id, {'IR'}));
nds_id = nodes_id(startsWith(nodes_id, {'ND' }));
meters_id = nodes_id(startsWith(nodes_id, {'Meter' }));
writecell({irs_id{:}; nds_id{:}; meters_id{:}}','NodeIDs.txt','Delimiter','\t')
ir_nd_ids = {irs_id{:}, nds_id{:}};
ir_nd_indices = wdn.getNodeIndex(ir_nd_ids);
ir_nd_meter_indices = wdn.getNodeIndex({irs_id{:}, nds_id{:}, meters_id{:}});
% remove basedemand of all nodes
wdn.setNodeBaseDemands(wdn.getNodeIndex,zeros(1,numel(wdn.getNodeIndex) ));
% Assign 1 as a base demand of all IR and ND nodes
wdn.setNodeBaseDemands(ir_nd_indices, ones(1,numel(ir_nd_indices) ) );
patterns_id = wdn.getPatternNameID;
patterns_index = wdn.getPatternIndex;
[ ~, loc ]=ismember(ir_nd_ids, patterns_id);
ir_nd_pattern_indices = patterns_index(loc);

% change the pattern to the length of simulation
pattern_length = wdn.getTimeSimulationDuration/wdn.getTimePatternStep + 1;
for i = ir_nd_pattern_indices
    wdn.setPattern(i, ones(1,pattern_length));
end

% assign the negative node a negative pattern
for id = nds_id
    wdn.setNodeDemandPatternIndex( wdn.getNodeIndex(id{:}), wdn.getPatternIndex(id{:}) )
end

patterns_id = wdn.getPatternNameID;
patterns_index = wdn.getPatternIndex;
[ ~, loc ]=ismember(ir_nd_ids, patterns_id);
ir_nd_pattern_indices = patterns_index(loc);

%
% [ ~, loc ]=ismember(ir_ids, pattern_ids);
% ir_patterns_index = pattern_indices(loc);
%
intinal_demand = .2*1.8* [ones(1, numel(ir_nd_pattern_indices)/2) zeros(1, numel(ir_nd_pattern_indices)/2)];
for node = 1 : numel(ir_nd_pattern_indices)
    wdn.setPatternValue( ir_nd_pattern_indices(node) , 1 , intinal_demand(node));
end
a = wdn.getPattern;
wdn.saveInputFile('temp.inp')

wdn.openHydraulicAnalysis;
wdn.initializeHydraulicAnalysis;
tstep=1;P=[];T_H=[];D=[];H=[];F=[];
counter = 1;

t = 1


while (tstep>0)
    
    if mode(t, wdn.getTimePatternStep) == 0
        counter = counter + 1;
        abm_demand = readmatrix('ABM Nodal Demand.txt','Delimiter','\t');
        abm_demand(:,3) = []; %remove the meter demands
        abm_demand = abm_demand(:);
        
        for node = 1 : numel(ir_nd_ids)
            wdn.setPatternValue( ir_nd_pattern_indices(node) , counter , abm_demand(node));
        end
    end
    t=wdn.runHydraulicAnalysis;
    P=[P; wdn.getNodePressure];
    D=[D; wdn.getNodeActualDemand];
    H=[H; wdn.getNodeHydaulicHead];
    F=[F; wdn.getLinkFlows];
    T_H=[T_H; t];
    tstep=wdn.nextHydraulicAnalysisStep;
    writematrix(reshape( P(end, ir_nd_meter_indices), [], 3) , 'Nodal Pressures.txt','Delimiter','\t');
end
wdn.closeHydraulicAnalysis;
wdn.unload;

%%
rng('default')
rng(1)
list = {'TRUE', 'FALSE'};
ind = datasample(1:2016, 200,'Replace',false);
log_ind = 2*ones(2016, 1);
log_ind(ind) = 1;
true_list = list(log_ind)';
writecell(true_list, fullfile('../data/input','producer boolean.txt'))
true_list2 = true_list;



%% simulate negative and irrigation files
wdn = epanet(fullfile('../data/network','sample_IRonly18houses_Modified Pump.inp'));
nodes_id = wdn.getNodeNameID;
nodes_index = wdn.getNodeIndex;
irs_id = nodes_id(startsWith(nodes_id, {'IR'}));
irs_ind = wdn.getNodeIndex(irs_id);
nds_id = nodes_id(startsWith(nodes_id, {'ND' }));
nds_ind = wdn.getNodeIndex(nds_id);
meters_id = nodes_id(startsWith(nodes_id, {'Meter' }));
writecell({irs_id{:}; nds_id{:}; meters_id{:}}', ...
    fullfile('../data/output','NodeIDs.txt'),'Delimiter','\t')
ir_nd_ids = {irs_id{:}, nds_id{:}};
ir_nd_indices = wdn.getNodeIndex(ir_nd_ids);
ir_nd_meter_indices = wdn.getNodeIndex({irs_id{:}, nds_id{:}, meters_id{:}});
% remove basedemand of all nodes
wdn.setNodeBaseDemands(wdn.getNodeIndex,zeros(1,numel(wdn.getNodeIndex) ));
% Assign 1 and -1 as a base demand of all IR and ND nodes
wdn.setNodeBaseDemands(irs_ind, ones(1,numel(irs_ind) ) );
wdn.setNodeBaseDemands(nds_ind, -1*ones(1,numel(nds_ind) ) );

patterns_id = wdn.getPatternNameID;
patterns_index = wdn.getPatternIndex;
[ logical, loc ]=ismember(patterns_id,ir_nd_ids );
ir_nd_pattern_indices = patterns_index(loc(logical));

% change the pattern to the length of simulation
folder = 'Produce 60 GPH_75';
if ~exist(fullfile('../data/input',folder), 'dir')
    mkdir(fullfile('../data/input',folder))
end
abm_demand_ir = readmatrix( fullfile('../data/input',folder, 'Irrigation Demand Output.txt'), 'Delimiter','\t');
abm_demand_nd = abs(readmatrix( fullfile('../data/input',folder, 'Negative Demand Output.txt'), 'Delimiter','\t') );
abm_demand = [abm_demand_ir; abm_demand_nd];
%  row = 2016;
% abm_demand(row,  find(abm_demand(row+2016,:)) )
abm_demand( :,isnan(abm_demand(1,:))) = [];
wdn.setTimeSimulationDuration(size(abm_demand,2)*3600);
pattern_length = wdn.getTimeSimulationDuration/wdn.getTimePatternStep + 1 ;
wdn.setPattern(1, ones(1,pattern_length)); pattern_matrix = wdn.getPattern;
wdn.setPatternMatrix( [ pattern_matrix(1:2,:) ; abs( [abm_demand(ir_nd_pattern_indices,24), abm_demand( ir_nd_pattern_indices,: )] )] );

% assign the negative node a negative pattern
for id = nds_id
    wdn.setNodeDemandPatternIndex( wdn.getNodeIndex(id{:}), wdn.getPatternIndex(id{:}) )
end
if ~exist(fullfile('../data/output',folder), 'dir')
    mkdir(fullfile('../data/output',folder))
end
wdn.saveInputFile(fullfile('../data/output',folder,'IRandND_18_ModifiedPump.inp') )
% remove negative ND nodes basedemand
wdn.setNodeBaseDemands(nds_ind, zeros(1, numel(nds_ind)))

abm_demand( find(log_ind == 1 ) , : ) = repmat(   abm_demand(2,:),200, 1);
wdn.setPatternMatrix( [ pattern_matrix(1:2,:) ; abs( [abm_demand(ir_nd_pattern_indices,24), abm_demand( ir_nd_pattern_indices,: )] )] );
wdn.saveInputFile(fullfile('../data/output',folder,'IR_only_18_ModifiedPump.inp') )
% abm_demand(row, : )
wdn.unload;

%%


wdn.openHydraulicAnalysis;
wdn.initializeHydraulicAnalysis;
tstep=1;P=[];T_H=[];D=[];H=[];F=[];
counter = 1;

t = 1


while (tstep>0)
    
    if mode(t, wdn.getTimePatternStep) == 0
        counter = counter + 1;
        abm_demand = readmatrix('ABM Nodal Demand.txt','Delimiter','\t');
        abm_demand(:,3) = []; %remove the meter demands
        abm_demand = abm_demand(:);
        
        for node = 1 : numel(ir_nd_ids)
            wdn.setPatternValue( ir_nd_pattern_indices(node) , counter , abm_demand(node));
        end
    end
    t=wdn.runHydraulicAnalysis;
    P=[P; wdn.getNodePressure];
    D=[D; wdn.getNodeActualDemand];
    H=[H; wdn.getNodeHydaulicHead];
    F=[F; wdn.getLinkFlows];
    T_H=[T_H; t];
    tstep=wdn.nextHydraulicAnalysisStep;
    writematrix(reshape( P(end, ir_nd_meter_indices), [], 3) , 'Nodal Pressures.txt','Delimiter','\t');
end
wdn.closeHydraulicAnalysis;
wdn.unload;

