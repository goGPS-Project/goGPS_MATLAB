%----------------------------------------------------------------------------------------------
% REPRESENTATION OF ROVER PSEUDORANGE OBSERVATIONS
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
