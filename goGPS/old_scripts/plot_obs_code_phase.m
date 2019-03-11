%----------------------------------------------------------------------------------------------
% REPRESENTATION OF PSEUDORANGE AND PHASE MEASUREMENT
%----------------------------------------------------------------------------------------------

%rover pseudorange
for i = 1 : nSatTot
    index = find(pr1_R(i,:) ~= 0)';
    if ~isempty(index)
        figure
        plot(index,pr1_R(i,index),'b.-'); grid on;
        title(['ROVER: PSEUDORANGE for SATELLITE ',num2str(i)]);
    end
end
%rover phase measurement
for i = 1 : nSatTot
    index = find(ph1_R(i,:) ~= 0)';
    if ~isempty(index)
        figure
        plot(index,ph1_R(i,index),'b.-'); grid on;
        title(['ROVER: PHASE MEASUREMENT for SATELLITE ',num2str(i)]);
    end
end

if (goGNSS.isDD(mode))
    %master pseudorange
    for i = 1 : nSatTot
        index = find(pr1_M(i,:) ~= 0)';
        if ~isempty(index)
            figure
            plot(index,pr1_M(i,index),'b.-'); grid on;
            title(['MASTER: PSEUDORANGE for SATELLITE ',num2str(i)]);
        end
    end
    %master phase measurement
    for i = 1 : nSatTot
        index = find(ph1_M(i,:) ~= 0)';
        if ~isempty(index)
            figure
            plot(index,ph1_M(i,index),'b.-'); grid on;
            title(['MASTER: PHASE MEASUREMENT for SATELLITE ',num2str(i)]);
        end
    end
end
