%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE COMBINATIONS OF ESTIMATED AMBIGUITIES
%----------------------------------------------------------------------------------------------

if (mode == goGNSS.MODE_PP_KF_CP_DD) | (mode == goGNSS.MODE_PP_KF_CP_SA)

    pivot_change = find(diff(pivot_OUT)~=0);

    for i = 1 : nSatTot
        index = find(conf_sat_OUT(i,:) == 1)';
        index_cs = find(conf_cs(i,:) == 1)';
        if ~isempty(index)
            index_pc = intersect(pivot_change(pivot_change>index(1)), pivot_change(pivot_change<index(end)));
            index_pc = intersect(index,index_pc);
            j = [1; find(index(2:end) - index(1:end-1) > 1)+1];
            figure
            %combination of estimated ambiguities
            plot(index,estim_amb(i,index),'b.-'); grid on;
            hold on
            %cycle-slip
            plot(index_cs,estim_amb(i,index_cs),'go');
            %combination of estimated ambiguities for new satellites
            plot(index(j),estim_amb(i,index(j)),'g.');
            %pivot change
            plot(index_pc,estim_amb(i,index_pc),'mo');
            %acceptance interval
            plot(index, estim_amb(i,index) + sigma_amb(i,index),'r:');
            plot(index, estim_amb(i,index) - sigma_amb(i,index),'r:');
            hold off
            title(['Combination of estimated ambiguities between PIVOT and SATELLITE ',num2str(i)]);
        end
    end

end
