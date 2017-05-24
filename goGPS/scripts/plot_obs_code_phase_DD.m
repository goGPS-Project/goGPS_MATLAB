%----------------------------------------------------------------------------------------------
% REPRESENTATION OF CODE AND PHASE DOUBLE DIFFERENCES
%----------------------------------------------------------------------------------------------

if (goGNSS.isDD(mode))
    %code
    for i = 1 : nSatTot
        index = find(conf_sat_OUT(i,:) == 1)';
        index_pivot = intersect(index, find(pivot_OUT == i));
        if ~isempty(index) & (length(index) ~= length(index_pivot))
            figure
            title(['CODE DOUBLE DIFFERENCES between SATELLITE ',num2str(i),'and PIVOT']);
            hold on
            grid on
            for j = 1 : length(time_GPS)
                if (conf_sat_OUT(i,j) & i ~= pivot_OUT(j))
                    plot(j,(pr1_R(i,j)-pr1_M(i,j))-(pr1_R(pivot_OUT(j),j)-pr1_M(pivot_OUT(j),j)),'b.');
                end
            end
        end
    end

    %phase
    for i = 1 : nSatTot
        index = find(conf_sat_OUT(i,:) == 1)';
        index_pivot = intersect(index, find(pivot_OUT == i));
        if ~isempty(index) & (length(index) ~= length(index_pivot))
            figure
            title(['PHASE DOUBLE DIFFERENCES between SATELLITE ',num2str(i),'and PIVOT']);
            hold on
            grid on
            for j = 1 : length(time_GPS)
                if (conf_sat_OUT(i,j) & i ~= pivot_OUT(j))
                    plot(j,(ph1_R(i,j)-ph1_M(i,j))-(ph1_R(pivot_OUT(j),j)-ph1_M(pivot_OUT(j),j)),'b.');
                end
            end
        end
    end
end
