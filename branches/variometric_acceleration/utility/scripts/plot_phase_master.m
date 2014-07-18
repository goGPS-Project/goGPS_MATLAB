%----------------------------------------------------------------------------------------------
% REPRESENTATION OF MASTER PHASE OBSERVATIONS
%----------------------------------------------------------------------------------------------

%master phase measurement
for i = 1 : nSatTot
    index = find(ph1_M(i,:) ~= 0)';
    if ~isempty(index)
        figure
        plot(index,ph1_M(i,index),'b.-'); grid on;
        title(['MASTER: PHASE MEASUREMENT for SATELLITE ',num2str(i)]);
    end
end