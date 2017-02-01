function [flagIntervals] = getOutliers(flags)
    % convert flags into array
    if isstruct(flags)
        flagArray = int8(struct2arrayFlags(flags, max(flags.pos+1)));
        % add padding to avoid problem with flags on the borders
    else
        flagArray = flags;
    end
    flagArray = flagArray(:);
    flagArray = [0; flagArray; 0];
    flagArray(flagArray ~= 0) = 1;
    diff = flagArray(1:end-1) - flagArray(2:end);
    clear flagArray;
    flagIntervals = [find(diff<0), find(diff>0)-1];
end

function [flagArray] = struct2arrayFlags(flags, maxSize)
flagArray = int8(zeros(maxSize,1));
flagArray(flags.pos) = flags.val;
end