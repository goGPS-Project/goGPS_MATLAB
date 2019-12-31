function cdata = loadIcon( filename, bgcol )
%loadIcon  Load an icon and set the transparent color
%
%   cdata = uix.loadIcon(filename) loads the icon from the specified
%   filename.  For PNG files with transparency, the transparent pixels are
%   set to NaN.  For other files, pixels that are pure green are set to
%   transparent (i.e., "green screen").  The resulting cdata is an RGB
%   double array.
%
%   cdata = uix.loadIcon(filename,bgcol) tries to merge the color data with
%   the specified background colour bgcol.  Fully transparent pixels are
%   still set to NaN, but partially transparent pixels are merged with the
%   background.
%
%   See also: imread

%  Copyright 2009-2016 The MathWorks, Inc.
%  $Revision: 1436 $ $Date: 2016-11-17 17:53:29 +0000 (Thu, 17 Nov 2016) $

% Check inputs
narginchk( 1, 2 )
if nargin < 2
    bgcol = get( 0, 'DefaultUIControlBackgroundColor' );
end

% First try normally
thisDir = fileparts( mfilename( 'fullpath' ) );
iconDir = fullfile( thisDir, 'Resources' );
if exist( filename, 'file' )
    [cdata, map, alpha] = imread( filename );
elseif exist( fullfile( iconDir, filename ), 'file' )
    [cdata, map, alpha] = imread( fullfile( iconDir, filename ) );
else
    error( 'uix:FileNotFound', 'Cannot open file ''%s''.', filename )
end

% Convert indexed images to RGB
if ~isempty( map )
    cdata = ind2rgb( cdata, map );
end

% Convert to double before applying transparency
cdata = convertToDouble( cdata );

% Handle transparency
[rows, cols, ~] = size( cdata );
if ~isempty( alpha )
    
    % Transparency specified
    alpha = convertToDouble( alpha );
    f = find( alpha==0 );
    if ~isempty( f )
        cdata(f) = NaN;
        cdata(f + rows*cols) = NaN;
        cdata(f + 2*rows*cols) = NaN;
    end
    % Now blend partial alphas
    f = find( alpha(:)>0 & alpha(:)<1 );
    if ~isempty(f)
        cdata(f) = cdata(f).*alpha(f) + bgcol(1)*(1-alpha(f));
        cdata(f + rows*cols) = cdata(f + rows*cols).*alpha(f) + bgcol(2)*(1-alpha(f));
        cdata(f + 2*rows*cols) = cdata(f + 2*rows*cols).*alpha(f) + bgcol(3)*(1-alpha(f));
    end
    
else
    
    % Do a "green screen", treating anything pure-green as transparent
    f = find( cdata(:,:,1)==0 & cdata(:,:,2)==1 & cdata(:,:,3)==0 );
    cdata(f) = NaN;
    cdata(f + rows*cols) = NaN;
    cdata(f + 2*rows*cols) = NaN;
    
end

end % uix.loadIcon

% -------------------------------------------------------------------------

function cdata = convertToDouble( cdata )
%convertToDouble  Convert image data to double in the range [0,1]
%
%  cdata = convertToDouble(cData)

switch lower( class( cdata ) )
    case 'double'
        % do nothing
    case 'single'
        cdata = double( cdata );
    case 'uint8'
        cdata = double( cdata ) / 255;
    case 'uint16'
        cdata = double( cdata ) / 65535;
    case 'int8'
        cdata = ( double( cdata ) + 128 ) / 255;
    case 'int16'
        cdata = ( double( cdata ) + 32768 ) / 65535;
    otherwise
        error( 'uix:InvalidArgument', ...
            'Image data of type ''%s'' is not supported.', class( cdata ) )
end

end % convertToDouble