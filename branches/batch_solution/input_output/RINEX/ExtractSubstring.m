function [ Substring ] = ExtractSubstring( String, BeginningIndex, EndingIndex )
%EXTRACTSUBSTRING Summary of this function goes here
%   Detailed explanation goes here
    
    BeginningIndex = max(BeginningIndex, 1);
    EndingIndex = min(EndingIndex, length(String));
    
    if (EndingIndex-BeginningIndex>0)
        Substring = String(BeginningIndex:EndingIndex);
    else
        Substring = [];
    end %if

end

