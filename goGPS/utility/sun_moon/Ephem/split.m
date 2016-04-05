function [ fr1, fr2 ] = split( tt )
% breaks up a double number into a double integer part and a fractional part
% 
% void split (double tt,
% 
%             double *fr)
% 
% ------------------------------------------------------------------------
% 
%    PURPOSE:
%       This function breaks up a double number into a double integer
%       part and a fractional part.
% 
%    REFERENCES:
%       Standish, E.M. and Newhall, X X (1988). "The JPL Export
%          Planetary Ephemeris"; JPL document dated 17 June 1988.
% 
%    INPUT
%    ARGUMENTS:
%       tt (double)
%          Input number.
% 
%    OUTPUT
%    ARGUMENTS:
%       *fr (double)
%          2-element output array;
%             fr[0] contains integer part,
%             fr[1] contains fractional part.
%          For negative input numbers,
%             fr[0] contains the next more negative integer;
%             fr[1] contains a positive fraction.
% 
%    RETURNED
%    VALUE:
%       None.
% 
%    GLOBALS
%    USED:
%       None.
% 
%    FUNCTIONS
%    CALLED:
%       None.
% 
%    VER./DATE/
%    PROGRAMMER:
%       V1.0/06-90/JAB (USNO/NA): CA coding standards
%       V1.1/03-93/WTH (USNO/AA): Convert to C.
%       V1.2/07-93/WTH (USNO/AA): Update to C standards.
%       V1.3/10-10/WKP (USNO/AA): Renamed function to lowercase to
%                                 comply with coding standards.
% 
%    NOTES:
%       None.
% 
% ------------------------------------------------------------------------

%   Get integer and fractional parts.

   ir = floor(tt);
   fr2 = tt - ir;
   fr1 = ir;

end

