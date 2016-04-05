function [ String ] = RemoveUnwantedTrailingSpaces( String )
%REMOVEUNWANTEDTRAILINGSPACES Summary of this function goes here
%   Detailed explanation goes here
    
    LengthOfString = length(String);
    IndexOfLatestCharInString = LengthOfString;
            
    while (IndexOfLatestCharInString >= 1) && isspace(String(IndexOfLatestCharInString))
        IndexOfLatestCharInString = IndexOfLatestCharInString-1;
    end %while
    
    if IndexOfLatestCharInString == 0
        String = [];
    else
        String = String(1:IndexOfLatestCharInString);
    end %if

end

