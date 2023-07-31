function varargout=dem(x,y,z,varargin)
%DEM Shaded relief image plot
%
%	DEM(X,Y,Z) plots the Digital Elevation Model defined by X and Y 
%	coordinate vectors and elevation matrix Z, as a lighted image using
%	specific "landcolor" and "seacolor" colormaps. DEM uses IMAGESC 
%	function which is much faster than SURFL when dealing with large 
%	high-resolution DEM. It produces also high-quality and moderate-size 
%	Postscript image adapted for publication.
%
%	DEM(X,Y,Z,'Param1',Value1,'Param2',Value2,...) specifies options or
%	parameter/value couple (case insensitive):
%
%	[H,I] = DEM(...); returns graphic handle H and optional illuminated 
%	image as I, a MxNx3 matrix (if Z is MxN and DECIM is 1).
%
%	I = DEM(...,'noplot') returns a structure I containing fields x, y, z,
%	and illuminated image .rgb without producing a graph on current figure.
%
%
%	--- Lighting options ---
%
%	'Azimuth',A
%		Light azimuth in degrees clockwise relative to North. Default is
%		A = -45 for	a natural northwestern illumination.
%
%	'Contrast',C
%		Light contrast, as the exponent of the gradient value:
%			C = 1 for linear contrast (default),
%			C = 0 to remove lighting,
%			C = 0.5 for moderate lighting,
%			C = 2 or more for strong contrast.
%
%	'LCut',LC
%		Lighting scale saturation cut with a median-style filter in % of 
%	    elements, such as LC% of maximum gradient values are ignored:
%			LC = 0.2 is default, 
%			LC = 0 for full scale gradient.
%
%	'km'
%		Stands that X and Y coordinates are in km instead of m (default).
%		This allows correct lighting. Ignored if LATLON option is used.
%
%
%	--- Elevation colorscale options ---
%
%	'ZLim',[ZMIN,ZMAX]
%		Fixes min and max elevation values for colormap. Use NaN to keep 
%		real min and/or max data values.
%
%	'ZCut',ZC
%		Median-style filter to cut extremes values of Z (in % of elements),
%		such that ZC% of most min/max elevation values are ignored in the
%		colormap application:
%			ZC = 0.5 is default, 
%			ZC = 0 for full scale.
%
%
%	--- "No Value" elevation options ---
%
%	'NoValue',NOVALUE
%		Defines the values that will be replaced by NaN. Note that values 
%		equal to minimum of Z class are automatically detected as NaN 
%		(e.g., -32768 for int16 class).
%
%	'NaNColor',[R,G,B]
%		Sets the RGB color for NaN/NoValue pixels (default is a dark gray).
%		Note that your must specify a valid 3-scalar vector (between 0 and
%		1);	color characters like 'w' or 'k' are not allowed, use [1,1,1]
%		or [0,0,0] instead.
%
%	'Interp'
%		Interpolates linearly all NaN values (fills the gaps using linear 
%		triangulation), using an optimized algorithm.
%
%
%	--- Colormap options ---
%
%	'LandColor',LMAP
%		Uses LMAP colormap instead of default (landcolor, if exists or 
%		jet) for Z > 0 elevations.
%
%	'SeaColor',SMAP
%		Sets the colormap used for Z <= 0 elevations. Default is seacolor 
%		(if exists) or single color [0.7,0.9,1] (a light cyan) to simulate
%		sea color.
%
%	'ColorMap',CMAP
%		Uses CMAP colormap for full range of elevations, instead of default 
%		land/sea. This option overwrites LANDCOLOR/SEACOLOR options.
%
%	'Lake'
%		Detects automatically flat areas different from sea level (non-zero 
%		elevations) and colors them as lake surfaces.
%
%	'LakeZmin',ZMIN
%		Activates the 'lake' option only above ZMIN elevations. For 
%		example, use 'lakezmin',0 to limit lake detection on land.
%
%	'Saturation',N
%		Changes the whole image color saturation by a factor of N.
%
%	'GrayScale'
%		Converts the used colormap(s) to grayscale (colour-blind).
%
%	'Watermark',N
%		Makes the whole image lighter by a factor of N.
%
%
%	--- Basemap and scale options ---
%
%	'Legend'
%		Adds legends to the right of graph: elevation scale (colorbar)
%		and a distance scale (in km).
%
%	'Cartesian'
%		Plots classic basemap-style axis, considering coordinates X and Y 
%		as cartesian in meters. Use parameter "km' for X/Y in km.
%
%	'LatLon'
%		Plots geographic basemap-style axis in deg/min/sec, considering 
%		coordinates X as longitude and Y as latitude. Axis aspect ratio 
%		will be adjusted to approximatively preserve distances (this is  
%		not a real projection!). This overwrites ZRatio option.
%
%	'AxisEqual', 'auto' (default) | 'manual' | 'off'
%		When 'Cartesian' or 'LatLon' option is used, automatic axes scaling
%		is applied to respect data aspect ratio. Default mode is 'auto' and
%		uses AXIS EQUAL and DASPECT functions. The 'manual' mode modifies
%		axes width or height with respect to the paper size in order to
%		produce correct data scaling at print (but not necessarily at 
%		screen). The 'off' mode disables any scaling.
%
%	Additionnal options for basemap CARTESIAN, LATLON, and LEGEND:
%
%	'BorderWidth',BW
%		Border width of the basemap axis, in % of axis height. Default is
%		BW = 1%.
%
%	'XTick',DX
%	'YTick',DY
%		X and Y Tick length (same unit as X and Y). Default is automatic.
%		Tick labels are every 2 ticks.
%
%	'FontSize',FS
%		Font size for X and Y tick labels. Default is FS = 10.
%
%	'FontBold'
%		Font weight bold for tick labels.
%
%	'Position',P
%		Position of the tick labels: 'southwest' (default), 'southeast',
%		'northwest','northeast'
%
%	'ZUnit',ZU
%		Sets the elevation unit as string ZU, used when min/max values are
%		indicated on colorbar legend (default is 'm').
%
%	--- Decimation options ---
%
%	For optimization purpose, DEM will automatically decimate data to limit
%	to a total of 1500x1500 pixels images. To avoid this, use following
%	options, but be aware that large grids may require huge computer 
%	ressources or induce disk swap or memory errors.
%
%	'Decim',N
%		Decimates matrix Z at 1/N times of the original sampling.
%       If N < 0, oversamples at -N rate.
%
%	'NoDecim'
%		Forces full resolution of Z, no decimation (N =1).
%
%
%
%	--- Informations ---
%
%	Colormaps are Mx3 RGB matrix so it is easy to modify contrast 
%	(CMAP.^N), set darker (CMAP/N), lighter (1 - 1/N + CMAP/N), inverse
%	it (flipud(CMAP)), etc...
%
%	To get free worldwide topographic data (SRTM), see READHGT function.
%
%	For backward compatibility, the former syntax is still accepted:
%	DEM(X,Y,Z,OPT,CMAP,NOVALUE,SEACOLOR) where OPT = [A,C,LC,ZMIN,ZMAX,ZC],
%	also option aliases DEC, DMS and SCALE, but there is no argument 
%	checking. Please prefer the param/value syntax.
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	Created: 2007-05-17 in Guadeloupe, French West Indies
%	Updated: 2021-01-08

%	History:
%	[2021-01-08] v2.11
%		- minor optimization for 'interp' option
%	[2020-11-29] v2.10
%		- new option 'grayscale'
%	[2020-11-25] v2.9
%		- fix a possible issue with saturation under Matlab < 2014
%	[2020-06-01] v2.8
%	    - new option 'saturation' to change colormaps
%	    - tick label decimals with smaller font size
%	[2019-06-17] v2.7
%		- fix an issue for single RGB color in sea/land color option
%	[2017-03-29] v2.6
%		- fix in 'lakezmin' option (thanks to Mustafa Çomo?lu)
%	[2017-01-09] v2.5
%		- new option 'lakezmin' to limit lake detection 
%	[2016-12-21] v2.4
%		- improves the colormap splitting between land and sea 
%	[2016-04-19] v2.3
%		- major update (thanks to mas Wiwit)
%	[2016-01-31] v2.2
%		- adds option 'Position' for tick labels
%	[2015-08-22] v2.1
%		- minor fix (former versions of Matlab compatibility)
%	[2015-08-19] v2.0
%		- image is now 100% true color (including the legend colorbar), 
%	      thus completely independent from the figure colormap
%	[2014-10-14]
%		- 'decim' option allows oversampling (negative value)
%	[2014-06-06]
%		- improves backward compatibility (adds strjoin subfunction)
%	[2014-03-18]
%		- adds new axisequal option
%	[2013-03-11]
%		- new options: 'km', 'watermark', 'fontsize', 'bordersize'
%		- improve legend colorbar
%		- all options now passed as param/value
%	[2013-01-14]
%		- improved light rendering (using surface normals instead of gradient)
%		- improved 'lake' detection algorithm
%		- new 'nancolor' option to set NaN color
%		- adds a length scale with 'dec' option
%		- minor code improvements
%	[2013-01-07]
%		- adds 'interp' option (fill the gaps)
%		- adds 'seacolor' colormap for negative elevations (bathymetry)
%	[2013-01-02]
%		- adds a 'lake' option
%		- minor bug correction
%	[2012-09-26]
%		- now accepts row/column vectors for X and/or Y.
%	[2012-05-29]
%		- adds basemap-style axis in decimal or lat/lon modes
%		- adds elevation and distance scales
%	[2012-05-18]
%		- new landcolor.m colormap function
%		- new arguments to control colormap scaling
%		- median-style filters for light and colormap
%	[2012-04-26]
%		- Optimizations: adds a decimation for large DEM grids.

%
%	Copyright (c) 2016, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

if nargin < 3
	error('Not enough input arguments.');
end

holdon = ishold;

degkm = 6378*pi/180; % one latitude degree in km
sea_color = [.7,.9,1]; % default sea color (light cyan)
grey = 0.2*[1,1,1]; % a dark gray


% -------------------------------------------------------------------------
% --- Manage input arguments

% number of arguments param/value
nargs = 0;

if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z)
	error('X,Y and Z must be numeric.')
end

if all(size(x) ~= 1) || all(size(y) ~= 1)
	error('X and Y must be vectors, not matrix.')
end

if length(x) ~= size(z,2) || length(y) ~= size(z,1)
	error('If Z has a size of [M,N], X must have a length of N, and Y a length of M.')
end

if size(z,3) == 3
	rgb = true;
else
	rgb = false;
end

% OPTIONS and PARAM/VALUE arguments
			
% AZIMUTH param/value
[s,az] = checkparam(varargin,'azimuth',@isscalar);
nargs = nargs + 2;
if s==0
	az = -45; % default
end

% ELEVATION param/value
[s,el] = checkparam(varargin,'elevation',@isscalar);
nargs = nargs + 2;
if s==0
	el = 0; % default
end

% CONTRAST param/value
[s,ct] = checkparam(varargin,'contrast',@isscalar);
nargs = nargs + 2;
if s
	ct = abs(ct);
else
	ct = 1; % default
end

% LCUT param/value
[s,lcut] = checkparam(varargin,'lcut',@isperc);
nargs = nargs + 2;
if s==0
	lcut = .2; % default
end

% NOVALUE param/value
[s,novalue] = checkparam(varargin,'novalue',@isscalar);
nargs = nargs + 2;
if s==0
	% default: min value for integer class / NaN for float
	S = whos('z');
	if strfind(S.class,'int')
		novalue = intmin(S.class);
	else
		novalue = NaN;
	end
end

% NANCOLOR param/value
[s,novalue_color] = checkparam(varargin,'nancolor',@isrgb);
nargs = nargs + 2;
if s==0
	novalue_color = grey; % default
end

% LANDCOLOR param/value
[s,cland] = checkparam(varargin,'landcolor',@isrgb);
nargs = nargs + 2;
if s==0
	% default: landcolor or jet
	if exist('landcolor','file')
		cland = landcolor.^1.3;
	else
		cland = jet(256);
	end
end

% SEACOLOR param/value
[s,csea] = checkparam(varargin,'seacolor',@isrgb);
nargs = nargs + 2;
if s==0
	% default: seacolor or single color
	if exist('seacolor','file')
		csea = seacolor;
	else
		csea = sea_color;
	end
end

% COLORMAP param/value
[s,cmap] = checkparam(varargin,'colormap',@isrgb);
nargs = nargs + 2;
if s
	cland = [];
	csea = [];
else
	% default
	cmap = cland;
end

% ZLIM param/value
[s,zmm] = checkparam(varargin,'zlim',@isvec);
nargs = nargs + 2;
if s
	zmin = min(zmm);
	zmax = max(zmm);
else
	zmin = NaN; % default
	zmax = NaN; % default
end

% ZCUT param/value
[s,zcut] = checkparam(varargin,'zcut',@isperc);
nargs = nargs + 2;
if s==0
	zcut = .5; % default
end

% ZRATIO param/value
[s,zratio] = checkparam(varargin,'zratio',@isscalar);
nargs = nargs + 2;
if s==0
	zratio = 1; % default
end

% SATURATION param/value
[s,csat] = checkparam(varargin,'saturation',@isscalar);
nargs = nargs + 2;
if s
	csat = abs(csat);
else
	csat = 1; % default
end

% WATERMARK param/value
[s,wmark] = checkparam(varargin,'watermark',@isscalar);
nargs = nargs + 2;
if s
	wmark = abs(wmark);
else
	wmark = 0; % default
end

% DECIM param/value and NODECIM option
[s,decim] = checkparam(varargin,'decim',@isscalar);
if s
	decim = round(decim);
	nargs = nargs + 2;
else
	decim = any(strcmpi(varargin,'nodecim')); % default
	nargs = nargs + 1;
end

% LAKEZMIN param/value option
[s,lakezmin] = checkparam(varargin,'lakezmin',@isscalar);
if s
    lake = 1;
	nargs = nargs + 2;
else
    lake = 0;
    lakezmin = NaN;
end

% FONTSIZE param/value
[s,fs] = checkparam(varargin,'fontsize',@isscalar);
nargs = nargs + 2;
if s==0
	fs = 10; % default
end

% BORDERWIDTH param/value
[s,bw] = checkparam(varargin,'borderwidth',@isperc);
nargs = nargs + 2;
if s==0
	bw = 1; % default
end

% XTICK param/value
[s,ddx] = checkparam(varargin,'xtick',@isscalar);
nargs = nargs + 2;
if s==0
	ddx = 0; % default (automatic)
end

% YTICK param/value
[s,ddy] = checkparam(varargin,'ytick',@isscalar);
nargs = nargs + 2;
if s==0
	ddy = 0; % default (automatic)
end

% POSITION param/value
[s,tpos] = checkparam(varargin,'position',@ischar,{'southwest','southeast','northwest','northeast'});
nargs = nargs + 2;
if s==0
	tpos = 'southwest'; % default
end

% ZUNIT param/value
[s,zunit] = checkparam(varargin,'zunit',@ischar);
nargs = nargs + 2;
if s==0
	zunit = 'm'; % default
end

% AXISEQUAL param/value
[s,axeq] = checkparam(varargin,'axisequal',@ischar,{'auto','manual','off'});
nargs = nargs + 2;
if s==0 || ~any(strcmpi(axeq,{'manual','off'}))
	axeq = 'auto'; % default (automatic)
end

% CROP param/value
[s,crop] = checkparam(varargin,'crop',@isvec,4);
nargs = nargs + 2;

% CLRGB param/value
[s,clrgb] = checkparam(varargin,'CLRGB',@isrgb);
nargs = nargs + 2;
if s==0
	clrgb = .5*ones(1,3); % default (mid-grey)
end

% CLLEVEL param/value
[s,cllevel] = checkparam(varargin,'CLLevel',@isvec,1:2);
nargs = nargs + 2;
if s==0
	cllevel = [0 0]; % default (auto)
end

% options without argument value
km = any(strcmpi(varargin,'km'));
dec = any(strcmpi(varargin,'cartesian') | strcmpi(varargin,'dec'));
dms = any(strcmpi(varargin,'latlon') | strcmpi(varargin,'dms'));
kmscale = any(strcmpi(varargin,'kmscale'));
scale = any(strcmpi(varargin,'legend') | strcmpi(varargin,'scale'));
inter = any(strcmpi(varargin,'interp'));
lake = any(strcmpi(varargin,'lake')) || lake;
fbold = any(strcmpi(varargin,'fontbold'));
noplot = any(strcmpi(varargin,'noplot'));
clines = any(strcmpi(varargin,'contourlines'));
gscale = any(strcmpi(varargin,'grayscale'));


% for backward compatibility (former syntax)...
nargs = nargs + dec + dms + scale + kmscale + inter + lake + km + fbold + noplot + clines + gscale;

if (nargin - nargs) > 3 && ~isempty(varargin{1})
	opt = varargin{1};
	if ~isnumeric(opt)
		error('OPT = [A,C,S,ZMIN,ZMAX,ZCUT] argument must be numeric.');
	end
	if ~isempty(opt)
		az = opt(1);
	end
	if length(opt) > 1
		ct = opt(2);
	end
	if length(opt) > 2
		lcut = opt(3);
	end
	if length(opt) > 4
		zmin = opt(4);
		zmax = opt(5);
	end
	if length(opt) > 5
		zcut = opt(6);
	end
end

if (nargin - nargs) > 4 && ~isempty(varargin{2})
	cmap = varargin{2};
	csea = [];
end

if (nargin - nargs) > 5 && ~isempty(varargin{3})
	novalue = varargin{3};
end

if (nargin - nargs) > 6 && ~isempty(varargin{4})
	csea = varargin{4};
end


% further test of input arguments
if dms && any(abs(y) > 91)
	error('With LATLON option Y must be in valid latitudes interval (decimal degrees).')
end

if km
	zratio = 1000;
end


% -------------------------------------------------------------------------
% --- Pre-process DEM data

% crops data if needed
if numel(crop)==4
	fprintf('DEM: crops original data from [%g,%g,%g,%g] to [%g,%g,%g,%g]...\n', ...
		min(x(:)),max(x(:)),min(y(:)),max(y(:)),crop);
	kx = find(x >= crop(1) & x <= crop(2));
	ky = find(y >= crop(3) & y <= crop(4));
	x = x(kx);
	y = y(ky);
	z = z(ky,kx,:);
end

% decimates data to avoid disk swap/out of memory...
nmax = 1500;
if decim
	n = decim;
else
	n = ceil(sqrt(numel(z))/nmax);
end
if n > 1
	x = x(1:n:end);
	y = y(1:n:end);
	z = z(1:n:end,1:n:end,:);
	fprintf('DEM: data has been decimated by a factor of %d...\n',n);
end

z = double(z); % necessary for most of the following calculations...
z(z==novalue) = NaN;

if isempty(csea)
	k = (z~=0 & ~isnan(z));
else
	k = ~isnan(z);
end

if isnan(zmin)
	zmin = nmedian(z(k),zcut/100);
end
if isnan(zmax)
	zmax = nmedian(z(k),1 - zcut/100);
end
dz = zmax - zmin;

if decim && n < 0
	xi = linspace(x(1),x(end),-n*length(x));
	yi = linspace(y(1),y(end),-n*length(y))';
	[xx,yy] = meshgrid(xi,yi);
	z = interp2(x,y,z,xx,yy,'*cubic');
	x = xi;
	y = yi;
	fprintf('DEM: data has been oversampled by a factor of %d...\n',-n);
end

if ~rgb && inter
	z = fillgap(x,y,z);
end

% -------------------------------------------------------------------------
% --- Process lighting

if ~rgb && dz > 0
	% first check if colormaps have the minimum required size
	if size(csea,1) < 2
		csea = repmat(csea,256,1);
	end
	if size(cland,1) < 2
		cland = repmat(cland,256,1);
	end
	% builds the colormap: concatenates seacolor and landcolor around 0
	% after interpolation to have exactly one color level per meter.
	if ~isempty(csea)
% 		l = size(csea,1);
% 		if zmin < 0 && zmax > 0
% 			r = size(cland,1)*abs(zmin)/zmax/l;
% 			cmap = cat(1,interp1(1:l,csea,linspace(1,l,ceil(l*r)),'*linear'),cland);
		if zmin < 0 && zmax > 0
			lcs = size(csea,1);
			lcl = size(cland,1);
			cmap = cat(1,interp1(1:lcs,csea,linspace(1,lcl,abs(zmin)+1),'*linear'), ...
				interp1(1:lcl,cland,linspace(1,lcl,abs(zmax)),'*linear'));
		elseif zmax <=0
			cmap = csea;
		end
	end
	
	% normalisation of Z using CMAP and convertion to RGB
	I = ind2rgb(uint16(round((z - zmin)*(size(cmap,1) - 1)/dz) + 1),cmap);
	
	if ct > 0
		% computes lighting from elevation gradient
		%[fx,fy] = gradient(z,x,y);
		if dms
			ryz = degkm*1000;
			rxz = degkm*1000*cosd(mean(y));
		else
			rxz = zratio;
			ryz = zratio;
		end
		[xx,yy] = meshgrid(x*rxz,y*ryz);
		[fx,fy,fz] = surfnorm(xx,yy,z);
		[ux,uy,uz] = sph2cart((90-az)*pi/180,el*pi/180,1);
		fxy = fx*ux + fy*uy + fz*uz;
		clear xx yy fx fy fz	% free some memory...
		
		fxy(isnan(fxy)) = 0;

		% computes maximum absolute gradient (median-style), normalizes,
		% saturates and duplicates in 3-D matrix
		li = 1 - abs(sind(el)); % light amplitude (experimental)
		r = repmat(max(min(li*fxy/nmedian(abs(fxy),1 - lcut/100),1),-1),[1,1,3]);
		rp = (1 - abs(r)).^ct;
	
		% applies contrast using exponent
		I = I.*rp;
	
		% lighter for positive gradient
		I(r>0) = I(r>0) + (1 - rp(r>0));
				
	end

	% set novalues / NaN to nancolor
	[i,j] = find(isnan(z));
	if ~isempty(i)
		I(sub2ind(size(I),repmat(i,1,3),repmat(j,1,3),repmat(1:3,size(i,1),1))) = repmat(novalue_color,size(i,1),1);
	end
	
	% lake option
	if lake
        klake = islake(z);
        if ~isnan(lakezmin)
            klake(z < lakezmin) = false; % removes indexes below ZMIN
        end
	else
		klake = 0;
	end
	
	% set the seacolor (upper color) for 0 values
	if ~isempty(csea)
		[i,j] = find(z==0 | klake);
		if ~isempty(i)
			I(sub2ind(size(I),repmat(i,1,3),repmat(j,1,3),repmat(1:3,size(i,1),1))) = repmat(csea(end,:),size(i,1),1);
		end
	end

	txt = '';
	
elseif rgb
	I = z;
	txt = '';

else
	I = repmat(shiftdim(sea_color,-1),size(z));
	cmap = repmat(sea_color,[256,1]);
	txt = 'Mak Byur!';	% Splash !
end

% -------------------------------------------------------------------------
% --- applies saturation, grayscale and watermark to image and cmap (for legend)

if gscale
	I = rgb2gray(I);
	cmap = rgb2gray(cmap);
end

if csat~=1
	I = saturation(I,csat);
	cmap = saturation(cmap,csat);
end

if wmark
	I = watermark(I,wmark);
	cmap = watermark(cmap,wmark);
end


% -------------------------------------------------------------------------
% --- ends the function when 'noplot' option is on
if noplot
	varargout{1} = struct('x',x,'y',y,'z',z,'rgb',I,'cmap',cmap);
	return
end

% -------------------------------------------------------------------------
% --- plots the RGB image
hh = imagesc(x,y,I);

if ~isempty(txt)
	text(mean(x),mean(y),txt,'Color',sea_color/4, ...
		'FontWeight','bold','HorizontalAlignment','center')
end

orient tall; axis xy
if strcmpi(axeq,'auto')
	axis equal
end
axis tight
xlim = [min(x),max(x)];
ylim = [min(y),max(y)];
zlim = [min([z(z(:) ~= novalue);zmin]),max([z(z(:) ~= novalue);zmax])];

if dms
	% approximates X-Y aspect ratio for this latitude (< 20-m precision for 1x1° grid)
	xyr = cos(mean(y)*pi/180);
else
	xyr = 1;
end

bw0 = max(diff(xlim)*xyr,diff(ylim))/100;
bwy = bw*bw0; % Y border width = 1%
bwx = bwy/xyr; % border width (in degree of longitude)


% -------------------------------------------------------------------------
% --- Axis basemap style
if dec || dms
	axis off

	if strcmpi(axeq,'manual')
		ppos = get(gcf,'PaperPosition');
		apos = get(gca,'Position');
		xyf = (xyr*diff(xlim)/apos(3)/ppos(3))/(diff(ylim)/apos(4)/ppos(4));
		if xyf >= 1
			set(gca,'Position',[apos(1),apos(2),apos(3),apos(4)/xyf]);
		else
			set(gca,'Position',[apos(1),apos(2),apos(3)*xyf,apos(4)]);
		end
	end
	if strcmpi(axeq,'auto')
		if diff(xlim)*xyr <= diff(ylim)
			set(gca,'DataAspectRatio',[1,xyr,1])
		else
			set(gca,'DataAspectRatio',[1/xyr,1,1])
		end
	end

	if bw > 0
		% transparent borders
		patch([xlim(1)-bwx,xlim(2)+bwx,xlim(2)+bwx,xlim(1)-bwx],ylim(1) - bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
		patch([xlim(1)-bwx,xlim(2)+bwx,xlim(2)+bwx,xlim(1)-bwx],ylim(2) + bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
		patch(xlim(1) - bwx*[0,0,1,1],[ylim(1)-bwy,ylim(2)+bwy,ylim(2)+bwy,ylim(1)-bwy],'k','FaceColor','none','clipping','off')
		patch(xlim(2) + bwx*[0,0,1,1],[ylim(1)-bwy,ylim(2)+bwy,ylim(2)+bwy,ylim(1)-bwy],'k','FaceColor','none','clipping','off')
	end
	dlon = {'E','W'};
	dlat = {'N','S'};
	if fbold
		fw = 'bold';
	else
		fw = 'normal';
	end
	
	if ddx == 0
		ddx = dtick(diff(xlim),dms);
		ddxn = 0;
	else
		ddxn = double(ddx<0);
		ddx = abs(ddx);
	end
	if ddy == 0
		ddy = dtick(diff(ylim),dms);
		ddyn = 0;
	else
		ddyn = double(ddy<0);
		ddy = abs(ddy);
	end
   	xtick = (ddx*ceil(xlim(1)/ddx)):ddx:xlim(2);
	for xt = xtick(1:2:end)
		dt = ddx - max(0,xt + ddx - xlim(2));
		patch(repmat(xt + dt*[0,1,1,0]',[1,2]),[ylim(1) - bwy*[0,0,1,1];ylim(2) + bwy*[0,0,1,1]]','k','clipping','off')
		if fs > 0
			if ~isempty(regexp(tpos,'north','once'))
				text(xt + dt*ddxn,ylim(2) + 1.2*bwy,deg2dms(xt + dt*ddxn,dlon,dec,fs),'FontSize',fs,'FontWeight',fw, ...
					'HorizontalAlignment','center','VerticalAlignment','bottom');
			else
				text(xt + dt*ddxn,ylim(1) - 1.2*bwy,deg2dms(xt + dt*ddxn,dlon,dec,fs),'FontSize',fs,'FontWeight',fw, ...
					'HorizontalAlignment','center','VerticalAlignment','top');
			end
		end
	end

	ytick = (ddy*ceil(ylim(1)/ddy)):ddy:ylim(2);
	for yt = ytick(1:2:end)
		dt = ddy - max(0,yt + ddy - ylim(2));
		patch([xlim(1) - bwx*[0,0,1,1];xlim(2) + bwx*[0,0,1,1]]',repmat(yt + dt*[0,1,1,0]',[1,2]),'k','clipping','off')
		if fs > 0
			if ~isempty(regexp(tpos,'east','once'))
				text(xlim(2) + 1.2*bwx,yt + dt*ddyn,deg2dms(yt + dt*ddyn,dlat,dec,fs),'FontSize',fs,'FontWeight',fw, ...
					'HorizontalAlignment','center','VerticalAlignment','top','rotation',90);
			else
				text(xlim(1) - 1.2*bwx,yt + dt*ddyn,deg2dms(yt + dt*ddyn,dlat,dec,fs),'FontSize',fs,'FontWeight',fw, ...
					'HorizontalAlignment','center','VerticalAlignment','bottom','rotation',90);
			end
		end
	end
end

% -------------------------------------------------------------------------
% --- Contour lines

% contour lines
if clines
	dz = diff(zlim);
	% empirical ratio between horizontal extent and elevation interval (dz)
	%rzh = dz/min(diff(x([1,end]))*cosd(mean(dlat)),diff(y([1,end])))/degkm/4e2;
	% major level lines
	if cllevel(1)==0
		dd = dtick(dz);
	else
		dd = abs(cllevel(1));
	end
	dz0 = ceil(zlim(1)/dd)*dd:dd:floor(zlim(2)/dd)*dd;
	dz0(ismember(0,dz0)) = [];	% eliminates 0 value
	% minor level lines
	if numel(cllevel)<2 || cllevel(2)==0
		dd = dtick(dz/5);
	else
		dd = abs(cllevel(2));
	end
	dz1 = ceil(zlim(1)/dd)*dd:dd:floor(zlim(2)/dd)*dd;
	dz1(ismember(dz1,dz0)) = [];	% eliminates minor ticks in major ticks
	hold on
	[~,h] = contour(x,y,z,[0,0,dz1],'-','Color',clrgb);
	set(h,'LineWidth',0.1);
	[cs,h] = contour(x,y,z,[0,0,dz0],'-','Color',clrgb);
	set(h,'LineWidth',1);
	if ~isempty(dz0)% && clineslabel
		clabel(cs,h,dz0,'Color',clrgb,'FontSize',fs/2,'FontWeight','bold', ...
			'LabelSpacing',288,'Margin',fs)
	end
	hold off
end

% -------------------------------------------------------------------------
% --- Scales legend

%wsc = diff(xlim)*0.01;
wsc = bw0;
xsc = xlim(2) + wsc*2 + bwx;

if scale

	% -- elevation scale (colorbar)
	zscale = linspace(zmin,zmax,length(cmap))';
	yscale = linspace(0,diff(ylim)/2,length(cmap));
	ddz = dtick(dz*max(0.5*xyr*diff(xlim)/yscale(end),1));
	ztick = (ddz*ceil(zscale(1)/ddz)):ddz:zscale(end);
 	rgbscale = ind2rgb(uint16(round((zscale - zmin)*(size(cmap,1) - 1)/dz) + 1),cmap);
	ysc = ylim(1);
	hold on
	imagesc(xsc + wsc*[-1,1]/2,ysc + yscale,repmat(rgbscale,1,2),'clipping','off');
	patch(xsc + wsc*[-1,1,1,-1],ysc + yscale(end)*[0,0,1,1],'k','FaceColor','none','Clipping','off')
	text(xsc + 2*wsc + zeros(size(ztick)),ysc + (ztick - zscale(1))*0.5*diff(ylim)/diff(zscale([1,end])),num2str(ztick'), ...
		'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',fs*.75)
	% indicates min and max Z values
	text(xsc,ysc - bwy/2,sprintf('%g %s',roundsd(zlim(1),3),zunit),'FontWeight','bold', ...
		'HorizontalAlignment','left','VerticalAlignment','top','FontSize',fs*.75)
	text(xsc,ysc + .5*diff(ylim) + bwy/2,sprintf('%g %s',roundsd(zlim(2),3),zunit),'FontWeight','bold', ...
		'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs*.75)
	
	% frees axes only if not hold on
	if ~holdon
		hold off
	end
	
end

if scale || kmscale
	% -- distance scale (in km)
	if dms
		fsc = degkm;
	else
		fsc = zratio/1e3;
	end
	dkm = dtick(diff(ylim)*fsc);
	ysc = ylim(2) - 0.5*dkm/fsc;
	if dkm > 1
		skm = sprintf('%g km',dkm);
	else
		skm = sprintf('%g m',dkm*1000);
	end
	if kmscale
		xsc = xlim(1) + wsc*2;
		ysc = ylim(1) + wsc*2;
		patch(xsc + dkm*[0,1,1,0]/fsc,ysc + wsc*[-1,-1,0,0],'k','FaceColor','w')
		for n = 0:2:(dkm-1)
			patch(xsc + (n + [0,1,1,0])/fsc,ysc + wsc*[-1,-1,0,0],'k','FaceColor','k')
		end
		text(xsc + .5*dkm/fsc,ysc,skm,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
			'Color','k','FontWeight','bold','FontSize',fs)
	else
		patch(xsc + wsc*[-1,-1,0,0],ysc + dkm*0.5*[-1,1,1,-1]/fsc,'k','FaceColor',grey,'clipping','off')
		text(xsc,ysc,skm,'rotation',-90,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
			'Color',grey,'FontWeight','bold','FontSize',fs*.75)
	end
end


if nargout > 0
	varargout{1} = hh;
end
if nargout > 1
	varargout{2} = I;
end
if nargout > 2
	varargout{3} = z;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nmedian(x,n)
%NMEDIAN Generalized median filter
%	NMEDIAN(X,N) sorts elemets of X and returns N-th value (N normalized).
%	So:
%	   N = 0 is minimum value
%	   N = 0.5 is median value
%	   N = 1 is maximum value

if nargin < 2
	n = 0.5;
end
y = sort(x(:));
y = interp1(sort(y),n*(length(y)-1) + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dd = dtick(dlim,deg)
%DTICK Tick intervals

if nargin < 2
	deg = 0;
end

if deg && dlim <= 2/60
	% less than 2 minutes: base 36
	m = 10^floor(log10(dlim*36))/36;
elseif deg && dlim <= 2
	% less than 2 degrees: base 6
	m = 10^floor(log10(dlim*6))/6;
else
	% more than few degrees or not degrees: decimal rules
	m = 10^floor(log10(dlim));
end
p = ceil(dlim/m);
if p <= 1
	dd = .1*m;
elseif p == 2
	dd = .2*m;
elseif p <= 5
	dd = .5*m;
else
	dd = m;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = deg2dms(x,ll,dec,fs)
%DEG2DMS Degree/minute/second display

% smaller font for decimals
fs2 = sprintf('\\fontsize{%d}',round(fs/1.25));

if dec
	s = regexprep(sprintf('%7.7g',x),'\.(.*)',sprintf('.{%s$1}',fs2));
else
	xa = abs(x) + 1/360000;
	%sd = sprintf('%d%c',floor(xa),176);	% ASCII char 176 is the degree sign
	sd = sprintf('%d°',floor(xa));
	sm = '';
	ss = '';
	if mod(x,1)
		sm = sprintf('%02d''',floor(mod(60*xa,60)));
		sa = floor(mod(3600*xa,60));
		if sa
			ss = sprintf('%02d"',sa);
		else
			if strcmp(sm,'00''')
				sm = '';
			end
		end
	end
	s = sd;
	if ~isempty(sm) || ~isempty(ss)
		s = [s,'{',fs2,sm,ss,'}'];
	end
	s = [s,ll{1+int8(x<0)}];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = watermark(x,n)
% makes colormap x watermark of ratio n
if nargin < 2
	n = 2;
end

if n == 0
    y = x;
else
    y = (x/n + 1 - 1/n);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=saturation(x,n)
% changes the color saturation by ratio n

% first needs to convert from RGB to HSV
if ndims(x) == 3
	r = x(:,:,1);
	g = x(:,:,2);
	b = x(:,:,3);
	mx = max(x,[],3); % max(r,g,b)
	mn = min(x,[],3); % min(r,g,b)
else
	r = x(:,1);
	g = x(:,2);
	b = x(:,3);
	mx = max(x,[],2); % max(r,g,b)
	mn = min(x,[],2); % min(r,g,b)
end

s = zeros(size(r));
h = zeros(size(r));
dt = mx - mn;
v = mx;

k = mx>0 & dt>0;
if any(k(:))
	s(k) = dt(k)./mx(k);
	kr = (k & mx==r);
	kg = (k & mx==g);
	kb = (k & mx==b);
	h(kr) = mod(60*(g(kr) - b(kr))./dt(kr),360);
	h(kg) = 60*(b(kg) - r(kg))./dt(kg) + 120;
	h(kb) = 60*(r(kb) - g(kb))./dt(kb) + 240;
end

% changes the saturation channel
s = s*n;
s(s>1) = 1;

% converting back from HSV to RGB
hi = mod(floor(h/60),6);
f = h/60 - hi;
l = v.*(1 - s);
m = v.*(1 - f.*s);
n = v.*(1 - (1 - f).*s);

% default for hi == 0
r = v; g = n; b = l;

k = hi==1;
r(k) = m(k); g(k) = v(k); b(k) = l(k);
k = hi==2;
r(k) = l(k); g(k) = v(k); b(k) =n(k);
k = hi==3;
r(k) = l(k); g(k) = m(k); b(k) = v(k);
k = hi==4;
r(k) = n(k); g(k) = l(k); b(k) = v(k);
k = hi==5;
r(k) = v(k); g(k) = l(k); b(k) = m(k);

y = cat(ndims(x),r,g,b);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=rgb2gray(x)
% removes color information

if ndims(x) == 3
	y = repmat(0.2989*x(:,:,1) + 0.5870*x(:,:,2) + 0.1140*x(:,:,3),1,1,3);
else
	y = repmat(0.2989*x(:,1) + 0.5870*x(:,2) + 0.1140*x(:,3),1,3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fillgap(x,y,z)
% This function reproduces the core of NANINTERP2

sz = size(z);
k = find(isnan(z));
k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
if ~isempty(k)
	[xx,yy] = meshgrid(x,y);
	mask = false(sz);
	k2 = ind90(sz,k); % k2 is linear index in the row order
	% sets to 1 every previous and next index, both in column and row order
	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = true; 
	mask(k) = false; % removes the novalue index
	z(k) = griddata(xx(mask),yy(mask),z(mask),xx(k),yy(k));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k2 = ind90(sz,k)

[i,j] = ind2sub(sz,k);
k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = islake(z)
% ISLAKE mask of zero gradient on 3x3 tiles
% We use diff matrix in row and column directions, and shift it to build
% a single vectorized test of surrounding pixels. To do this we must
% concatenate unit vectors in different combinations...

dx = diff(z,1,2);	% differences in X direction
dy = diff(z,1,1);	% differences in Y direction
u1 = ones(size(z,1),1);	% row unit vector 
u2 = ones(1,size(z,2));	% column unit vector
u2r = u2(2:end);

% index of the tiles center pixel
k = ( ...
	[u2;dy] == 0 & [dy;u2] == 0 & ...
	[u1,dx] == 0 & [dx,u1] == 0 & ...
	[u1,[dx(2:end,:);u2r]] == 0 & [[dx(2:end,:);u2r],u1] == 0 & ...
	[u1,[u2r;dx(1:end-1,:)]] == 0 & [[u2r;dx(1:end-1,:)],u1] == 0 ...
);

% now extends it to surrounding pixels
k(1:end-1,:) = (k(1:end-1,:) | k(2:end,:));
k(2:end,:) = (k(2:end,:) | k(1:end-1,:));
k(:,1:end-1) = (k(:,1:end-1) | k(:,2:end));
k(:,2:end) = (k(:,2:end) | k(:,1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isrgb(x,n)

if nargin < 2
	n = 0;
end
if isnumeric(x) && (n == 1 && all(size(x) == [1,3]) || n == 0 && size(x,2) == 3) ...
		&& all(x(:) >= 0 & x(:) <= 1)
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isperc(x)

if isnumeric(x) && isscalar(x) && x >= 0 && x <= 100
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isvec(x,n)

if nargin < 2
	n = 2;
end
if isnumeric(x) && (nargin > 1 && any(numel(x) == n) || numel(x) > 1)
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=roundsd(x,n)

og = 10.^(floor(log10(abs(x)) - n + 1));
y = round(x./og).*og;
y(x==0) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,v] = checkparam(arg,nam,func,val)

switch func2str(func)
	case 'isscalar'
		num = 1;
		mes = 'scalar value';
	case 'isperc'
		num = 1;
		mes = 'percentage scalar value';
	case 'isvec'
		num = 1;
		if nargin < 4
			val = 2;
		end
		mes = sprintf('%d-element vector',val);
	case 'isrgb'
		num = 1;
		mes = '[R,G,B] vector with 0.0 to 1.0 values';
	case 'ischar'
		num = 0;
		mes = 'string';
		if nargin > 3
			mes = sprintf('%s (%s)',mes,strjoin(val,' or '));
		end
	otherwise
		num = 1;
		mes = 'value';
end

s = 0;
v = [];
k = find(strcmpi(arg,nam));
if ~isempty(k)
	if (k + 1) <= length(arg) ...
			&& (~num || isnumeric(arg{k+1})) ...
			&& (nargin < 4 && func(arg{k+1}) ...
				|| (nargin > 3 && (strcmp(func2str(func),'ischar') && ismember(arg{k+1},val)) ...
					 || strcmp(func2str(func),'isvec') && func(arg{k+1},val)))
		v = arg{k+1};
		s = 1;
	else
		error('%s option must be followed by a valid %s.',upper(nam),mes)
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=strjoin(c,d)
%STRJOIN Join cell array of strings
%(this is for Matlab versions < 2013a backward compatibility)

if nargin < 2
	d = '';
end
n = numel(c);
ss = cell(2,n);
ss(1,:) = reshape(c,1,n);
ss(2,1:n-1) = {d};
s = [ss{:}];
