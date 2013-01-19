function [ String ] = PadStringWithSpaces( String, RequestedLengthOfPaddedString )
    NumberOfSpacesToAdd = RequestedLengthOfPaddedString - length(String);
    if NumberOfSpacesToAdd > 0
        String = [String, repmat(' ',1,NumberOfSpacesToAdd)];
    end %if
end
