clc
clear
close all

%% add EPANET matlab toolkit function
addpath(genpath( '../' ))
disp('All files and folder Paths Loaded.');


%% Read the Network
folder = 'Produce 60 GPH_75';
wdn = epanet(fullfile('../data/output',folder,'IRandND_18_ModifiedPump.inp'));

% find the procumers' ID from the the network
patterns_value = wdn.getPattern();
node_pattern_index = wdn.getNodeDemandPatternIndex();
node_pattern_index = node_pattern_index{:};
nodes_basedemand = wdn.getNodeBaseDemands();
nodes_basedemand = nodes_basedemand{:};




procumers = readcell( 'producer boolean.txt' ); % 
nodes_id_read = readcell('NodeIDs.txt', 'Delimiter','\t' );
temp_ind = ismember( procumers,'TRUE' );
procumers_id = nodes_id_read(temp_ind,2);

pattern_ind = wdn.getPatternIndex();


for h = 1 : numel(procumers_id)
    h_id = procumers_id{ h };
    pattern_ind = wdn.getPatternIndex(h_id);
    wdn.getPatternValue(h_id);
end

nodes_base_demand = wdn.getNodeBaseDemands();
procumer = find(nodes_base_demand{1} == -1 );



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
