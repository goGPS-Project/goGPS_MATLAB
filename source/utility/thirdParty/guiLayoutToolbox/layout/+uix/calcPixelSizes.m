function pSizes = calcPixelSizes( pTotal, mSizes, pMinima, pPadding, pSpacing )
%calcPixelSizes  Calculate child sizes in pixels
%
%  pSizes = uix.calcPixelSizes(total,mSizes,minSizes,padding,spacing)
%  computes child sizes (in pixels) given total available size (in pixels),
%  child sizes (in pixels and/or relative), minimum child sizes (in
%  pixels), padding (in pixels) and spacing (in pixels).
%
%  Notes:
%  * All children are at least as large as the minimum specified size
%  * Relative sizes are respected for children larger than then minimum
%  specified size
%  * Children may extend beyond the total available size if the minimum
%  sizes, padding and spacing are too large

%  Copyright 2009-2015 The MathWorks, Inc.
%  $Revision: 1182 $ $Date: 2015-12-07 14:27:30 -0500 (Mon, 07 Dec 2015) $

% Initialize
pSizes = NaN( size( mSizes ) ); % output
n = numel( mSizes ); % need this later

% Apply absolute sizes
a = mSizes >= 0; % absolute
pSizes(a) = max( mSizes(a), pMinima(a) );

while true
    
    u = isnan( pSizes ); % unsolved
    pUnsolvedTotal = pTotal - max( (n-1), 0 ) * pSpacing ...
        - 2 * sign( n ) * pPadding - sum( pSizes(~u) );
    pUnsolvedSizes = mSizes(u) / sum( mSizes(u) ) * pUnsolvedTotal;
    pUnsolvedMinima = pMinima(u);
    s = pUnsolvedSizes < pUnsolvedMinima; % small
    if any( s )
        pUnsolvedSizes(s) = pUnsolvedMinima(s);
        pUnsolvedSizes(~s) = NaN;
        pSizes(u) = pUnsolvedSizes;
        % repeat
    else
        pSizes(u) = pUnsolvedSizes;
        break % done
    end
    
end

end % calcPixelSizes