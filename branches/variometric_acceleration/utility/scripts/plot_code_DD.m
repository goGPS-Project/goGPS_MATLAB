%----------------------------------------------------------------------------------------------
% REPRESENTATION OF CODE DOUBLE DIFFERENCES
%----------------------------------------------------------------------------------------------

for i = 1 : nSatTot
    index = find(conf_sat(i,:) == 1)';
    index_pivot = intersect(index, find(pivot == i));
    if ~isempty(index) & (length(index) ~= length(index_pivot))
        figure
        title(['CODE DOUBLE DIFFERENCES between SATELLITE ',num2str(i),'and PIVOT']);
        hold on
        grid on
        for j = 1 : length(time_GPS)
            if (conf_sat(i,j) & i ~= pivot(j) & pivot(j) ~= 0)
                try
                plot(j,(pr1_R(i,j)-pr1_M(i,j))-(pr1_R(pivot(j),j)-pr1_M(pivot(j),j)),'b.');
                catch
                    keyboard
                end
            end
        end
    end
end