function addHousehold( WDN, tn_index, pp_node, pp_pipe, scale_factor, theta )

% get host node information tn
TN_ID = WDN.getNodeNameID{tn_index};
TN_ind = tn_index;
TN_coor = WDN.getNodeCoordinates( TN_ind );
TN_Elv = WDN.getNodeElevations( TN_ind );


% add ND and IR demand
for patt = find( [pp_node.Pattern_ind] )
    pattID = [pp_node(patt).ID '_' TN_ID '_' num2str(theta)];
    WDN.addPattern( pattID );
    pp_node(patt).WDN_PattIND = WDN.getPatternIndex( pattID );
    WDN.setPattern( pp_node(patt).WDN_PattIND, pp_node(patt).Pattern );
end


% add PP's junctions
for n =  find( [pp_node.Type] == 0 ) % run for nodes after the meter
    %     if ~strcmpi( pp_node(n).ID, 'RE1' )
    n_name = [ pp_node(n).ID '_' TN_ID  '_' num2str(theta)];
    WDN.addNodeJunction( n_name );
    n_ind = WDN.getNodeIndex( n_name );
    WDN.setNodeCoordinates( n_ind, TN_coor + scale_factor * rotation(pp_node(n).X, pp_node(n).Y, theta)' );
    WDN.setNodeElevations( n_ind, TN_Elv + pp_node(n).Elv );
    if pp_node(n).BDemand ~= 0
        WDN.setNodeBaseDemands( n_ind, pp_node(n).BDemand );
        WDN.setNodeDemandPatternIndex( n_ind, pp_node(n).WDN_PattIND );
    end
    %     end
end


% add pipes
source = pp_node( [pp_node.Type] == 1 ).ID; % reterive source ID
for p = 1 : length( pp_pipe )
    i = pp_pipe( p ).i;
    j = pp_pipe( p ).j;
    
    if strcmpi( source, i)
        i_ID = TN_ID;
        j_ID = [ j '_' TN_ID  '_' num2str(theta)];
    elseif strcmpi( source, j)
        i_ID = [ i '_' TN_ID  '_' num2str(theta)];
        j_ID = TN_ID;
    else
        i_ID = [ i '_' TN_ID  '_' num2str(theta)];
        j_ID = [ j '_' TN_ID  '_' num2str(theta)];
    end
    p_ID = [ pp_pipe( p ).ID '_' TN_ID '_' num2str(theta) ];
    
    if pp_pipe( p ).type <= 1 % pipe or CV pipes
        if pp_pipe( p ).type == 1
            WDN.addLinkPipe( p_ID, i_ID, j_ID );
        else
            WDN.addLinkPipeCV( p_ID, i_ID, j_ID );
        end
        p_ind = WDN.getLinkIndex( p_ID );
        WDN.setLinkLength( p_ind, pp_pipe( p ).L );
        WDN.setLinkDiameter( p_ind, pp_pipe( p ).D );
        WDN.setLinkRoughnessCoeff( p_ind, pp_pipe( p ).C )
        
    else % this is valve
        WDN.addLinkValvePRV( p_ID, i_ID, j_ID );
        p_ind = WDN.getLinkIndex( p_ID );
        PP.setLinkInitialSetting( p_ind, pp_valve( p ).Setting );
        WDN.setLinkDiameter( p_ind, pp_pipe( p ).D );
    end
end
end
