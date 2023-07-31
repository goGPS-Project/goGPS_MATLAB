function varargout = readhgt(varargin)
%READHGT Import/download NASA SRTM data files (.HGT).
%	READHGT(AREA) where AREA is a 4-element vector [LAT1,LAT2,LON1,LON2]
%	downloads the SRTM data and plots a map corresponding to the geographic
%	area defined by latitude and longitude limits (in decimal degrees).
%	Upper-right coordinates limits (LAT2,LON2) are excluded.
%	If the needed SRTM .hgt files are not found in the current directory (or 
%	in the path), they are downloaded from the USGS data server (needs an 
%	Internet connection and a companion file "readhgt_srtm_index.txt"). For
%	better plot results, it is recommended to install DEM personal function
%	available at author's Matlab page. 
%
%	READHGT(LAT,LON) reads or downloads the SRTM tiles corresponding to LAT
%	and LON (in decimal degrees) coordinates (lower-left corner). LAT and/or
%	LON can be vectors: in that case, tiles corresponding to all possible 
%	combinations of LAT and LON values will be downloaded, and optional 
%	output structure X will have as much elements as tiles.
%
%	READHGT(FILENAME) reads single .hgt data file FILENAME, must be in the
%	form "[N|S]yy[E|W]xxx.hgt[.zip]", as downloaded from SRTM data servers.
%
%	X=READHGT(...) returns a structure X containing: 
%		lat: coordinate vector of latitudes (decimal degree)
%		lon: coordinate vector of longitudes (decimal degree)
%		  z: matrix of elevations (meters, INT16 class)
%		hgt: downloaded filename(s)
%
%	X=READHGT(...,'plot') also plots the tile(s).
%
%
%	--- Additionnal options ---
%
%	'tiles'
%	   Imports and plots individual tiles instead of merging them (default 
%	   behavior if adjoining values of LAT and LON).
%
%	'interp'
%	   Linearly interpolates missing data.
%
%	'decim',N
%	   Decimates the tiles at 1/N times of the original sampling. Plot
%	   is automatically decimated if necessary. Use N=1 to force full
%	   resolution (might induces memory issue for large areas).
%
%	'crop'
%	   crops the resulting map around existing land (reduces any sea or 
%	   novalue areas at the borders).
%
%	'crop',[LAT1,lAT2,LON1,LON2]
%	   Former syntax that crops the map using latitude/longitude limits. 
%	   Prefer the new syntax READHGT(AREA).
%
%	'outdir',OUTDIR
%	   Specifies output directory OUTDIR to write downloaded files and/or 
%	   to search existing files. Former syntax READHGT(LAT,LON,OUTDIR) also
%	   accepted. Default is the current directory.
%
%	'url',URL
%	   Specifies the URL address to find HGT files (default is USGS). 
%	   Former syntax READHGT(LAT,LON,OUTDIR,URL) still accepted.
%
%	'srtm3'
%	   Forces SRTM3 download for all areas (by default, SRTM1 tiles are 
%	   downloaded only for USA territory, if exists).
%
%	'srtm1'
%	   Uses SRTM1 tiles which are 9 times bigger than default SRTM3. The
%	   tiles must be available in the OUTDIR directory or they will be
%	   downloaded from NASA/EarthDATA website. This action needs a valid
%	   user login (see LOGIN option below).
%	   ! Beware with large zones may lead to computer memory issues.
%	   ! SRTM1 and SRTM3 tiles hold the same filename, while they have 
%	   different size. Do not store them in the same directory to avoid
%	   errors when merging tiles.
%
%	'login',USER,PASSWORD
%	  Needed user authentification to download SRTM1 tiles from the
%	  NASA/EarthDATA center. See https://urs.earthdata.nasa.gov to register.
%
%	'wget'
%	   Will use external command wget to download the files (for Linux and
%	   MacOSX systems). Because access to SRTM data uses a secured https://
%	   URL, Matlab versions older than 2016a will fail to download the 
%	   tiles automatically due to unzip function and certificate problems. 
%	   This issue can be surrounded using this 'wget' option and installing
%	   wget command on your system. For MacOSX, wget must be installed 
%	   using homebrew, macports, fink or compiling from sources.
%	   Also with this option, the login USER/PASSWORD can be stored in your
%	   home .netrc instead of 'login' arguments, by adding these lines:
%	      machine urs.earthdata.nasa.gov
%	      user xxxxx
%	      password xxxxx
%
%
%	--- Examples ---
%
%	- to plot a map of the Paris region, France (single tile):
%		readhgt([48,49,2,3])
%
%	- to plot a map of Flores volcanic island, Indonesia (5 tiles):
%		readhgt([-9,-8,119,124])
%
%	- to plot a map of the Misti volcano, Peru (SRTM1 cropped tile, needs
%	   a valid login user and password at NASA/EarthDATA center):
%	   readhgt([-16.4,-16.2,-71.5,-71.3],'srtm1','login','xx','xx','interp')
%
%	- to download SRTM1 data of Cascade Range (27 individual tiles):
%		X=readhgt(40:48,-123:-121,'tiles');
%
%
%	--- Information ---
%
%	- each file corresponds to a tile of 1x1 degree of a square grid
%	  1201x1201 of elevation values (SRTM3 = 3 arc-seconds), and for USA  
%	  territory or when using the 'srtm1' option, at higher resolution 
%	  3601x3601 grid (SRTM1 = 1 arc-second). Note that SRTM1 and SRTM3 
%	  files have the same syntax names; only the size differs.
%
%	- elevations are of class INT16: sea level values are 0, unknown values
%	  equal -32768 (there is no NaN for INT class), use 'interp' option to
%	  fill the gaps.
%
%	- note that borders are included in each tile, so to concatenate tiles
%	  you must remove one row/column in the corresponding direction (this
%	  is made automatically by READHGT when merging tiles).
%
%	- downloaded file is written in the current directory or optional  
%	  OUTDIR directory, and it remains there. Take care that mixed SRTM1
%	  and SRTM3 files may lead to fail to merge. It is better to use
%	  different directories for SRTM1 and SRTM3 (see 'outdir' option).
%
%	- NASA Shuttle Radar Topography Mission [February 11 to 22, 2000] 
%	  produced a near-global covering on Earth land, but still limited to 
%	  latitudes from 60S to 60N. Offshore tiles will be output as flat 0
%	  value grid.
%
%	- if you look for other global topographic data, take a look to ASTER
%	  GDEM, worldwide 1 arc-second resolution (from 83S to 83N): 
%	  http://gdex.cr.usgs.gov/gdex/ (free registration required)
%
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%		Institut de Physique du Globe de Paris
%
%	References:
%		https://srtm.kurviger.de
%
%	Acknowledgments: Yves Gaudemer, Jinkui Zhu, Greg
%
%	Created: 2012-04-22 in Paris, France
%	Updated: 2021-08-11

%	Copyright (c) 2021, François Beauducel, covered by BSD License.
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

fidx = 'readhgt_srtm_index.txt';
% ATTENTION: this file must exist in the Matlab path to use default SRTM3 tiles
% since USGS delivers data continent-by-continent with nominative directories,
% this index file is needed to know the full path name of each tile.
sz1 = [3601,3601]; % SRTM1 tile size
sz3 = [1201,1201]; % SRTM3 tile size
novalue = intmin('int16'); % -32768
pixel = 1/sz1(1); % minimum pixel size
n = 1;

srtm1 = any(strcmpi(varargin,'srtm1'));
if srtm1
	% SRTM1 full resolution tiles available here with EarthData login (2020):
	url = 'https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11';
	%url = 'http://e4ftl01.cr.usgs.gov/MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11';
	%url = 'http://rmd.neoknet.com/srtm1'; % unavailable since May 2020...
else
	% SRTM3 tiles (and SRTM1 for USA):
	%url = 'https://dds.cr.usgs.gov/srtm/version2_1'; % unavailable since 2021...
	url = 'https://srtm.kurviger.de'; % USGS mirror site
end

srtm3 = any(strcmpi(varargin,'srtm3'));
makeplot = any(strcmpi(varargin,'plot'));
merge = any(strcmpi(varargin,'merge'));	% unused former option but needs to be considered as valid
tiles = any(strcmpi(varargin,'tiles'));
inter = any(strcmpi(varargin,'interp'));
wget = any(strcmpi(varargin,'wget'));

% --- option: 'crop' or 'crop',[LAT1,LAT2,LON1,LON2]
crop = [];
cropflag = 0;
kcrop = find(strcmpi(varargin,'crop'));
if ~isempty(kcrop)
	cropflag = 1;
	if (kcrop + 1) <= nargin && isnumeric(varargin{kcrop+1})
		crop = varargin{kcrop+1};
		if any(size(crop) ~= [1,4])
			error('CROP option arguments must be a 1x4 vector.')
		end
		cropflag = 2;
	end
end

% --- option: 'decim',N
decim = 0;
decimflag = 0;
kdecim = find(strcmpi(varargin,'decim'));
if ~isempty(kdecim)
	decimflag = 1;
	if (kdecim + 1) <= nargin && isnumeric(varargin{kdecim+1})
		decim = round(varargin{kdecim+1});
		if ~isscalar(decim) || decim < 1
			error('DECIM option argument must be a positive integer.')
		end
		decimflag = 2;
	end
end

% --- option: 'outdir',OUTDIR
out = '.';
outflag = 0;
koutdir = find(strcmpi(varargin,'outdir'));
if ~isempty(koutdir)
	if (koutdir + 1) <= nargin && ischar(varargin{koutdir+1})
		outflag = 2;
		out = varargin{koutdir+1};
		if ~exist(out,'dir')
			error('OUTDIR is not a valid directory.')
		end
	else
		error('''outdir'' option must be followed by OUTDIR string.')
	end
end

% --- option: 'url',URL
urlflag = 0;
kurl = find(strcmpi(varargin,'url'));
if ~isempty(kurl)
	if (kurl + 1) <= nargin && ischar(varargin{kurl+1})
		urlflag = 2;
		url = varargin{kurl+1};
	else
		error('''url'' option must be followed by URL string.')
	end
end

% --- option: 'login',USER,PASSWORD
loginflag = 0;
klogin = find(strcmpi(varargin,'login'));
if ~isempty(klogin)
	if (klogin + 2) <= nargin && ischar(varargin{klogin+1}) && ischar(varargin{klogin+2})
		loginflag = 3;
		usr = varargin{klogin+1};
		pwd = varargin{klogin+2};
	else
		error('''login'' option must be followed by USER and PASSWORD strings.')
	end
end

% needs to count the arguments to allow former syntaxes...
nargs = makeplot + merge + tiles + cropflag + srtm1 + srtm3 + inter ...
	+ decimflag + outflag + urlflag + loginflag + wget;

% syntax READHGT without argument: opens the GUI to select a file
if nargin == nargs
	[filename,pathname] = uigetfile('*.hgt;*.hgt.zip','Select a HGT file');
	f = {[pathname,filename]};
	if filename == 0
		error('Please select a HGT file or use function arguments.');
	end
end

% syntax READHGT(FILENAME, ...)
if nargin == (1 + nargs) && ischar(varargin{1})
	f = varargin{1};
	if ~exist(f,'file')
		error('FILENAME must be a valid file name')
	end
	[~,filename] = fileparts(f);
	f = {f};
end

if nargin < (2 + nargs) && exist('filename','var')
	lat = str2double(filename(2:3));
	if filename(1) == 'S'
		lat = -lat;
	end
	lon = str2double(filename(5:7));
	if filename(4) == 'W'
		lon = -lon;
	end
else
	if nargin < (2 + nargs)
		crop = varargin{1};
		if ~isnumeric(crop) || any(size(crop) ~= [1,4])
			error('Area must be a 4-element vector [LAT1,LAT2,LON1,LON2].')
		end
		lat = floor(min(crop(1:2))):floor(max(crop(1:2))-pixel);
		lon = floor(min(normlon(crop(3:4)))):floor(max(normlon(crop(3:4)))-pixel);
		cropflag = 2;
	else
		lat = floor(varargin{1}(:));
		lon = normlon(floor(varargin{2}(:)));	% longitudes are normilized to -180/+179 interval
	end
	if ~isnumeric(lon) || ~isnumeric(lat) || any(abs(lat) > 60) || any(lon < -180) || any(lon > 179) || isempty(lat) || isempty(lon)
		error('LAT and LON must be numeric and in valid SRTM interval (abs(LAT)<60).');
	end
	if ~tiles && (any(diff(lat) ~= 1) || any(diff(lon) ~= 1))
		fprintf('READHGT: Warning! LAT and LON vectors do not define adjoining tiles. Cannot merge and force TILES option.');
		tiles = 1;
	end

	% former syntax: readhgt(LAT,LON,OUTDIR)
	if nargin > (2 + nargs)
	if ~isempty(varargin{3})
			out = varargin{3};
			if ~exist(varargin{3},'dir')
				error('OUTDIR is not a valid directory.')
			end
		end
	end
	
	% wget option: automatic setting
	mrel = version('-release');
	if issorted({mrel,'2014b'}) && ~wget
		wget = 1;
		fprintf('** Warning ** Matlab release %s: ''wget'' option forced.\n',mrel);
	end
	
	% if LAT/LON are vectors, NDGRID makes a grid of corresponding tiles
	[lat,lon] = ndgrid(lat,lon);
	f = cell(size(lat));
	for n = 1:numel(f)
		if lat(n) < 0
			slat = sprintf('S%02d',-lat(n));
		else
			slat = sprintf('N%02d',lat(n));
		end
		if lon(n) < 0
			slon = sprintf('W%03d',-lon(n));
		else
			slon = sprintf('E%03d',lon(n));
		end
		f{n} = sprintf('%s/%s%s.hgt',out,slat,slon);
		
		if ~exist(f{n},'file')
			ff = '';
			% former syntax: readght(LAT,LON,OUTDIR,URL)
			if nargin > (3 + nargs)
				url = varargin{4};
				if ~ischar(url)
					error('URL must be a string.');
				end
			else
				if srtm1
					ff = sprintf('/%s%s.SRTMGL1.hgt.zip',slat,slon);
					%ff = sprintf('/%s%s.hgt.zip',slat,slon);
				else
					%fsrtm = sprintf('%s/%s',fileparts(mfilename('fullpath')),fidx);
					fsrtm = fidx;
					if exist(fsrtm,'file')
						fid = fopen(fsrtm,'rt');
						idx = textscan(fid,'%s');
						fclose(fid);
						k = find(~cellfun('isempty',strfind(idx{1},sprintf('%s%s',slat,slon))));
						if isempty(k)
							%fprintf('READHGT: Warning! Cannot find %s tile in SRTM database. Consider it offshore...\n',ff);
						else
							% forcing SRTM3 option: takes the first match in the list
							if srtm3
								ff = idx{1}{k(1)};
							else
								ff = idx{1}{k(end)};
							end
						end
					else
						error('Cannot find "%s" index file to parse SRTM database. Please download HGT file manually.',fsrtm);
					end
				end
			end
			if isempty(ff)
				f{n} = '';
			else
				fprintf('Download %s%s ... ',url,ff);
				if srtm1 && (~exist('usr','var') || ~exist('pwd','var'))
					error('SRTM1 new tiles download needs the LOGIN option.')
				end
				f{n} = '';
				try
					if wget
						if system('which wget')
							fprintf(' ** wget binary not found. Cannot download tiles.\n');
						else
							tmp = tempname;
							mkdir(tmp)
							if srtm1
								login = sprintf('--user %s --password %s',usr,pwd);
							else
								login = '';
							end
							ftmp = sprintf('%s/%s%s.hgt.zip',tmp,slat,slon);
							[s,w] = system(sprintf('wget %s -O %s %s%s',login,ftmp,url,ff));
							if s
								disp(w)
							end
							f(n) = unzip(ftmp,out);
							delete(ftmp)
						end
					else
						if exist('websave','file')
							ftmp = [tempname(out),'.zip'];
							if srtm1
								opt = weboptions('Username',usr,'Password',pwd);
							else
								opt = weboptions;
							end
							websave(ftmp,[url,ff],opt);
							f(n) = unzip(ftmp,out);
							delete(ftmp)
						else
							f(n) = unzip([url,ff],out);
						end
					end
					fprintf('done.\n');
				catch
					fprintf(' ** tile not found. Considering offshore.\n');
				end
			end
		end
	end
end

% pre-allocates X structure (for each file/tile)
X = repmat(struct('hgt',[],'lat',[],'lon',[]),[n,1]);

if n == 1
	tiles = 0;
end

for n = 1:numel(f)
	% unzips HGT file if needed
	if ~isempty(strfind(f{n},'.zip'))
		X(n).hgt = char(unzip(f{n}));
		funzip = 1;
	else
		X(n).hgt = f{n};
		funzip = 0;
	end

	if srtm1
		sz = sz1;
	else
		sz = sz3;
	end
	if isempty(f{n})
		% offshore: empty tile...
		X(n).z = [];
	else
		% loads data from HGT file
		fid = fopen(X(n).hgt,'rb','ieee-be');
			X(n).z = fread(fid,'*int16');
		fclose(fid);
		switch numel(X(n).z)
		case prod(sz1)
			% srtm3 option: decimates the tile...
			if srtm3
				z = reshape(X(n).z,sz1);
				X(n).z = z(1:3:end,1:3:end);
				sz = sz3;
			else
				sz = sz1;
			end
		case prod(sz3)
			sz = sz3;
		otherwise
			error('"%s" seems not a regular SRTM data file or is corrupted.',X(n).hgt);
		end
		X(n).z = rot90(reshape(X(n).z,sz));

		% erases unzipped file if necessary
		if (funzip)
			delete(f{n});
		end
	end

	% builds latitude and longitude coordinates
	X(n).lon = linspace(lon(n),lon(n)+1,sz(2));
	X(n).lat = linspace(lat(n),lat(n)+1,sz(1))';
	
	% interpolates NaN (if not merged)
	if inter && tiles
		X(n).z = fillgap(X(n).lon,X(n).lat,X(n).z,novalue);
	end
end

if ~tiles
	% NOTE: cannot merge mixted SRTM1 / SRTM3 or discontiguous tiles
	Y.lat = linspace(min(lat(:)),max(lat(:))+1,size(lat,1)*(sz(1)-1)+1)';
	Y.lon = linspace(min(lon(:)),max(lon(:))+1,size(lon,2)*(sz(2)-1)+1);
	Y.z = zeros(length(Y.lat),length(Y.lon),'int16');
	for n = 1:numel(X)
		if ~isempty(X(n).z)
			Y.z((sz(1)-1)*(X(n).lat(1)-Y.lat(1)) + (1:sz(1)),(sz(2)-1)*(X(n).lon(1)-Y.lon(1)) + (1:sz(2))) = X(n).z;
		end
	end

	if cropflag
		if cropflag == 1 || isempty(crop)
			klat = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,2));
			klon = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,1));
		else
			crop = [minmax(crop(1:2)),normlon(minmax(crop(3:4)))];
			klat = find(Y.lat >= crop(1) & Y.lat <= crop(2));
			klon = find(Y.lon >= crop(3) & Y.lon <= crop(4));
		end			
		Y.lat = Y.lat(klat);
		Y.lon = Y.lon(klon);
		Y.z = Y.z(klat,klon);
	end
	
	if inter
		Y.z = fillgap(Y.lon,Y.lat,Y.z,novalue);
	end
end

if nargout == 0 || makeplot
	if ~tiles
		fplot(Y.lon,Y.lat,Y.z,decim,url,novalue)
	else
		for n = 1:numel(X)
			fplot(X(n).lon,X(n).lat,X(n).z,decim,url,novalue)
		end
	end
end

if nargout == 3 % for backward compatibility...
	varargout{1} = X(1).lon;
	varargout{2} = X(1).lat;
	varargout{3} = X(1).z;
elseif nargout > 0
	if tiles
		if decim > 1
			for n = 1:numel(X)
				X(n).lon = X(n).lon(1:decim:end);
				X(n).lat = X(n).lat(1:decim:end);
				X(n).z = X(n).z(1:decim:end,1:decim:end);
			end
		end
		varargout{1} = X;
	else
		if decim > 1
			Y.lon = Y.lon(1:decim:end);
			Y.lat = Y.lat(1:decim:end);
			Y.z = Y.z(1:decim:end,1:decim:end);
		end
		varargout{1} = Y;
	end
	if nargout == 2
		varargout{2} = f{1}; % for backward compatibility...
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fplot(x,y,z,decim,url,novalue)
%FPLOT plot the data using DEM function if exists, or IMAGESC

demoptions = {'latlon','legend','lake','nodecim'};

figure
if decim
	n = decim;
else
	n = ceil(sqrt(numel(z))/1201);
end
if n > 1
	x = x(1:n:end);
	y = y(1:n:end);
	z = z(1:n:end,1:n:end);
	fprintf('READHGT: In the figure data has been decimated by a factor of %d...\n',n);
end

if exist('dem','file')
	dem(x,y,z,demoptions{:})
else
	warning('For better results you might install the function dem.m from http://www.ipgp.fr/~beaudu/matlab.html#DEM')
	z(z==novalue) = 0;
	imagesc(x,y,z);
	if exist('landcolor','file')
		colormap(landcolor(256).^1.3)
	else
		colormap(jet)
	end
	% aspect ratio (lat/lon) is adjusted with mean latitude
	xyr = cos(mean(y)*pi/180);
	set(gca,'DataAspectRatio',[1,xyr,1])

	orient tall
	axis xy, axis tight
end

title(sprintf('Data SRTM/NASA from %s',url),'FontSize',10,'Interpreter','none')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = firstlast(x)

k = find(x);
y = k(1):k(end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = minmax(x)

y = [min(x(:)),max(x(:))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fillgap(x,y,z,novalue)
% GRIDDATA is not efficient for large arrays, but has great advantage to be
% included in Matlab core functions! To optimize interpolation, we
% reduce the number of relevant data by building a mask of surrounding
% pixels of novalue areas... playing with linear index!

sz = size(z);
k = find(z == novalue);
k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
if ~isempty(k)
	[xx,yy] = meshgrid(x,y);
	mask = zeros(sz,'int8');
	k2 = ind90(sz,k); % k2 is linear index in the row order
	% sets to 1 every previous and next index, both in column and row order
	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = 1; 
	mask(k) = 0; % removes the novalue index
	kb = find(mask); % keeps only border values
	z(k) = int16(griddata(xx(kb),yy(kb),double(z(kb)),xx(k),yy(k)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k2 = ind90(sz,k)

[i,j] = ind2sub(sz,k);
k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = normlon(x)
% normalize longitude between -180 and 180
y = mod(x+180,360) - 180;
