function [ nodes, pipes ] = readWDN( water_network )

% Read all pattern values and Indices
k = 0;
for i  = water_network.NodePatternIndex
    k = k + 1;
    if i == 0
        patternID{ k } = [];
        patternValue { k } = [];
    else
        patternID{ k } = water_network.PatternNameID{i};
        patternValue{ k } = water_network.Pattern(i,:);
    end
end
% nodes information
nodes = struct(...
    'ID',water_network.NodeNameID',...
    'X', num2cell( water_network.NodeCoordinates{1} ),...
    'Y', num2cell( water_network.NodeCoordinates{2} ),...
    'BDemand', num2cell(water_network.NodeBaseDemands{:}'),...
    'Pattern_ind', num2cell( water_network.NodePatternIndex' ),...
    'PatternID', patternID',...
    'Pattern', patternValue' ,...
    'Elv', num2cell(water_network.NodeElevations'),...
    'Type', num2cell(water_network.NodeTypeIndex') );

% pipes information
pipes = struct(...
    'ID', water_network.LinkNameID',...
    'Xver', water_network.NodeCoordinates{3},...
    'Yver', water_network.NodeCoordinates{4},...
    'i', water_network.NodesConnectingLinksID(:,1),...
    'j', water_network.NodesConnectingLinksID(:,2),...
    'L', num2cell(water_network.LinkLength') ,...
    'D', num2cell(water_network.LinkDiameter'),...
    'C', num2cell(water_network.LinkRoughnessCoeff'),...
    'loss', num2cell(water_network.LinkMinorLossCoeff'),...
    'state', num2cell(water_network.LinkInitialStatus'),...
    'type', num2cell( water_network.LinkTypeIndex' ),...
    'bulk', num2cell(water_network.LinkBulkReactionCoeff'),...
    'wall', num2cell(water_network.LinkWallReactionCoeff') );
end
