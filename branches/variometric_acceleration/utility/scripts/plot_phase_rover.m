%----------------------------------------------------------------------------------------------
% REPRESENTATION OF ROVER PHASE OBSERVATIONS
%----------------------------------------------------------------------------------------------

%rover phase measurement
for i = 1 : nSatTot
    index = find(ph1_R(i,:) ~= 0)';
    if ~isempty(index)
        figure
        plot(index,ph1_R(i,index),'b.-'); grid on;
        title(['ROVER: PHASE MEASUREMENT for SATELLITE ',num2str(i)]);
    end
end
