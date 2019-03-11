%----------------------------------------------------------------------------------------------
% REPRESENTATION OF LAMBDA1*L1-P1
%----------------------------------------------------------------------------------------------

%rover
for i = 1 : nSatTot
    index = find(conf_sat_OUT(i,:) == 1)';
    if ~isempty(index)
        figure
        plot(index,goGNSS.LAMBDA1*ph1_R(i,index)-pr1_R(i,index),'b.-'); grid on;
        title(['ROVER: lambda1*L1-P1 for SATELLITE ',num2str(i)]);
    end
end

%master
if (goGNSS.isDD(mode))
    for i = 1 : nSatTot
        index = find(conf_sat_OUT(i,:) == 1)';
        if ~isempty(index)
            figure
            plot(index,goGNSS.LAMBDA1*ph1_M(i,index)-pr1_M(i,index),'b.-'); grid on;
            title(['MASTER: lambda1*L1-P1 for SATELLITE ',num2str(i)]);
        end
    end
end
