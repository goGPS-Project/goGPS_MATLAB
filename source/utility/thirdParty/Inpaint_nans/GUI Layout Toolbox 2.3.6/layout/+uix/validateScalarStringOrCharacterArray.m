function value = validateScalarStringOrCharacterArray( value, propertyName )
%VALIDATESCALARSTRINGORCHARACTERARRAY Verify that the given value is a 
%scalar string or a character array.

% Check if we have a string scalar.
if isa( value, 'string' ) && isscalar( value )
    % Convert the string to a character array, handling the missing case
    % separately.
    if ismissing( value )
        value = '';
    else
        value = char( value );
    end % if
end % if

% Check that we have a character array.
try
    assert( ischar( value ), ...
        'uix:InvalidPropertyValue', ...
        ['Property ''', propertyName, ''' must be a scalar ', ...
        'string or a character array.'] )
catch e
    e.throwAsCaller()
end % try/catch

end % validateScalarStringOrCharacterArray