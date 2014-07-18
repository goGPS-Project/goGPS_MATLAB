%----------------------------------------------------------------------------------------------
% REPRESENTATION OF MASTER PSEUDORANGE OBSERVATIONS
%----------------------------------------------------------------------------------------------

%master pseudorange
for i = 1 : nSatTot
    index = find(pr1_M(i,:) ~= 0)';
    if ~isempty(index)
        figure
        plot(index,pr1_M(i,index),'b.-'); grid on;
        title(['MASTER: PSEUDORANGE for SATELLITE ',num2str(i)]);
    end
end
