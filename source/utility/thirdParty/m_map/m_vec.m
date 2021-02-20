function [hp, ht] = m_vec(s, x, y, varargin)  % and color, etc
% M_VEC Draws fancy arrows/quiverplots on a map.
%    [HP, HT]=M_VEC(S,LONG,LAT,VARARGIN) draws arrows as a patch 
%    object on a map created by the M_Map package.
%      HP: handle to the patch object
%      HT: handle to text object if this is a key; see below.
%      S:  scale factor, arrow data units per inch.
%          It is based on the figure PaperPosition, so be
%          sure to set this (e.g. via "landscape" or "portrait")
%          before calling m_vec.
%      LONG,LAT: longitude, latitude of the arrows
%          LONG and LAT must have the same dimensions, but can be
%          scalars or vectors; if scalars, multiple arrows can
%          be plotted at a single location.
%      VARARGIN can consist of any of the following:
%        U,V   U,V,C   Z,U,V,C
%      followed by optional arrow parameters,
%      followed by optional patch parameters.
%        U,V are vectors containing the east and north components
%                 of the arrows.
%        C is an optional colorspec for all arrows, or an
%                 array of CData, one value per arrow.
%                 Defaults to black.
%        Z is a height in axes data units: this is subject
%                 to future modification or omission, and
%                 is probably not useful now as-is.
%      optional arrow parameters: keyword-value pairs, shown
%                 here with default values:
%         'headangle',60     degrees: angle of arrow tip
%         'headwidth',NaN    points: direct specification of
%                                width, instead of headangle
%         'headlength',5     points: length of tip; set to 0
%                                to omit arrowhead entirely
%         'shaftwidth',1     points: width of arrow shaft
%         'centered', 'no'   'yes' to make x,y the arrow
%                                center instead of its tail
%         'key', ''          make a labelled horizontal arrow
%                                if the string is not empty;
%                                then the string labels the
%                                arrow, and the second argument
%                                returned, ht, is the handle of
%                                the string.
%          'edgeclip', 'off'  If 'on' then arrows IN the axes
%                             are clipped if their heads are
%                             OUT of the axes.
%
%      optional patch parameters: any valid patch properties
%         may be specified here; they are passed directly
%         to the patch function.
%
%    M_VEC called without any parameters generates a demonstration plot.

%
%  Mon  98/02/16 Eric Firing, efiring@soest.hawaii.edu
%
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 2/May/01 - small bug fix (thx to Pierre Jaccard)
% 7/jun/06 - arrows near boundaries were not done correctly - now fixed.
% 12/Dec/12 - added clipping for arrows at boundaries ('edgeclip')

global MAP_PROJECTION

% Have to have initialized a map first


if nargin==0,  % demo

    % demvec.m
    % This is a demonstration of m_vec.
    % Sun  98/02/22 Eric Firing

    % Set up the figure and axes at the start, so that the vector lengths
    % will come out right when the figure is printed:

  %  orient tall
    % Main axes, for the map.
    ha1 = axes;
    pos = get(ha1,'position');
    pos1 = pos;
    pos1(4) = pos1(4) - 0.1;
    pos1(2) = pos1(2) + 0.1;
    set(ha1,'position',pos1)

    % Colormap axes.  Make them smaller than for the default colormap.
    pos2 = [0.3  0  0.4  0];
    pos2(2) = pos(2);
    pos2(4) = 0.02;
    % This positioning of the colormap is actually not very good as
    % printed, but I am not going to fiddle with it any more now.
    % This is just a typical difficulty with normalized coordinates
    % combined with fixed DataAspectRatio; if
    % something looks reasonable on the screen, it is likely to look
    % wrong when printed.

    %ha2 = axes('position', pos2);

    axes(ha1); % Back to the main axes.

    m_proj('ortho','lat',48','long',-123', 'rad', 10, ...
        	'rec', 'off' );
    m_coast('patch',[0.9 0.95 0.9]);
    m_grid('linestyle','-','xtick',[-135:5:-110],'linewi',2);
    title('Demonstration of m\_vec')


    %% The following form, centered and without heads, is suitable for
    %% the square roots of the principle axes of variance ellipses,
    %% for example.  Color is given as a single RGB triplet; if it were
    %% given as an array of triplets, the two lines could be different
    %% colors; unfortunately, as of Matlab 5.1, if you do this the plot
    %% will be forced into Zbuffer mode with the warning:
    %%  "RGB color data not yet supported in Painter's mode"
    %% This example also illustrates specification of a black patch
    %% boundary, via direct patch property specifications following any
    %% vector parameter options.
    hpv1 = m_vec(100, [-133 -133], [49 49], [0 50], [100 0.0],...
	  [0.7 0.8 0.9],'centered','yes', ...
	  'shaftwidth', 5, 'headlength', 0,...
	  'EdgeColor','k');

    %% Here are three vectors, no color specified, all properties defaulted.
    hpv2 = m_vec(100, [-128 -128 -128], [46 46 46], ...
     [0 25*sqrt(2) 50], [50 25*sqrt(2) 0]);

    incs = (1:20)/20;
    vlat = 42 + incs*3;
    vlon = -127 - incs*2;
    uu = 50*sin(incs*2*pi);
    vv = 50*cos(incs*2*pi);

    %% This simulates a set of vectors such as currents measured
    %% with an ADCP along a cruise track.  Color is specified by letter,
    %% in this case. The arrows and heads are about as skinny as they
    %% can reasonably be.
    hpv3 = m_vec(100, vlon, vlat+2, uu, vv, 'm', ...
     'shaftwidth', 0.2, 'headlength', 2.5);

    %% Now let's make a similar vectors, but with colors based on the
    %% colormap, simulating SST, for example:
    vlon = vlon - 3;
    sst = 12 - incs * 4;
    hpv4 = m_vec(100, vlon-4, vlat, uu, vv, sst);
    ha2=colorbar('southoutside');
    set(get(ha2,'xlabel'),'string','SST');
    %axes(ha2); xlabel('SST'); axes(ha1);
    %% These vectors will not show up well on the screen because they
    %% have no edges, but they will print adequately.  There does not
    %% seem to be any easy way to get around this; I would have to
    %% completely change the way the patches are specified.  Ordinarily,
    %% one would probably use fatter than default vectors with this
    %% method anyway, so that the colors are clear when printed; and in
    %% this case the screen display will look OK also.

    % Here I show the edgeclip property
    hpv4 = m_vec(100, vlon+4, vlat-5, uu, vv, sst,'edgeclip','on');


    %% Key: Note that it can be on or off the map. The vector and
    %% text are the same color by default. Handles to both the vector
    %% (patch) and the label (text) are returned so that you can
    %% change the font, color, etc.
    [hpv5, htv5] = m_vec(100, -115, 38, 50, 0, 'b', 'key', '50 cm s^{-1}');
    set(htv5,'FontSize',8);

    return

end











if isempty(MAP_PROJECTION),
  error('No Map Projection initialized - call M_PROJ first!');
end;



% Default arrow parameters:
centered = 0;
headlength = 5/72;
headwidth  = NaN;
headangle = 40;
shaftwidth = 1/72;
c = 'k';
key = '';

clip = 'on'; % except for the key
edgeclip = 'off';  % for arrows at the edge

if nargin < 5,   % Minimum: s,x,y,u,v
   help('mvec');
   error('not enough input arguments');
end

UVperIn = s;
x = x(:);
y = y(:);

%% Begin the slightly complicated parsing of arguments.  This
%% could probably be done more efficiently and elegantly.
% One cause of complexity is that the optional argument c
% could be either character or numeric.
keyvars = {};  % This will hold keyword/value pairs.
nvarargin = length(varargin);
% Find the index of the first string argument in varargin.
istr0 = 0;
for ii = 1:nvarargin
   if ischar(varargin{ii}), istr0 = ii; break, end
end

if istr0 > 0,    % There are strings.
   if rem(nvarargin - istr0, 2) == 0,  % an odd number
      c = varargin{istr0}; % should be a colorspec string
      istr0 = istr0 + 1;
   end
   % istr0: start of keyword/value pairs.
   n_numeric = istr0 - 1;  %Actually, numeric or colorspec.
   keyvars = varargin(istr0:nvarargin);
   finished = 0;
   while ~isempty(keyvars) && ~finished,
      kv = lower(keyvars{1});
      value = keyvars{2};
      if     strcmp(kv, 'headlength')
         headlength = value/72;
      elseif strcmp(kv, 'headwidth')
         headwidth = value/72;
      elseif strcmp(kv, 'headangle')
         headangle = value;

      elseif strcmp(kv, 'shaftwidth')
         shaftwidth = value/72;
      elseif strcmp(kv, 'centered')  && lower(value(1)) == 'y',
         centered = 1;
      elseif strcmp(kv, 'key')
         key = value;
         clip = 'off'; % Can put key outside the map.
      elseif strcmp(kv, 'edgeclip')
         edgeclip = value;
      else
         finished = 1;  % no match; break out
      end
      if ~finished,
         keyvars(1:2) = [];
      end
   end
else
   n_numeric = nvarargin;
end

% Calculate the headwidth if it is not given explicitly:
if isnan(headwidth) && headangle < 170 && headangle > 0
   headwidth = headlength * tan(headangle*pi/180);
end
headwidth = max([headwidth; shaftwidth]);

if n_numeric == 2 || n_numeric == 3,
   u = varargin{1}(:);
   v = varargin{2}(:);
   if n_numeric == 3,
      c = varargin{3};
   end
   z = zeros(size(u));
elseif n_numeric == 4,
   z = varargin{1}(:);
   u = varargin{2}(:);
   v = varargin{3}(:);
   c = varargin{4};
else
   help('mvec')
   error('not enough numeric arguments')
end

[nr,nc] = size(c);
if nr == 1 && nc == length(u) && (nc ~= 3 || (any(c<=1) || any(c>=0))),
   c = c(:);
end
% c could be a 1x3 colorspec

if (length(x) == 1 && length(y) == 1 && length(u) > 1)
   x = x(ones(size(u)));
   y = y(ones(size(u)));
end

%% End of input argument parsing.

OrigAxUnits = get(gca,'Units');
if OrigAxUnits(1:3) == 'nor'
   OrigPaUnits = get(gcf, 'paperunits');
   set(gcf, 'paperunits', 'inches');
   figposInches = get(gcf, 'paperposition');
   set(gcf, 'paperunits', OrigPaUnits);
   axposNor = get(gca, 'position');
   axWidLenInches = axposNor(3:4) .* figposInches(3:4);
else
   set(gca, 'units', 'inches');
   axposInches = get(gca, 'position');
   set(gca, 'units', OrigAxUnits);
   axWidLenInches = axposInches(3:4);
end

% Multiply inches by the following to get data units:
scX = diff(get(gca, 'XLim'))/axWidLenInches(1);
scY = diff(get(gca, 'YLim'))/axWidLenInches(2);
sc = max([scX;scY]);  %max selects the dimension limited by
                      % the plot box.

Width = shaftwidth*sc;
HeadWidth = headwidth*sc;
HeadLength = headlength*sc;

uvmag = abs(u + i*v);
% Arrow lengths in plot data units:
L = uvmag*sc/UVperIn;

% base of arrows. Don't plot if outside boundaries (except for keys for which clip is off)
[xs, ys] = m_ll2xy(x,y, 'clip', clip);

%[xsp, ysp] = m_ll2xy(x+0.1*u./uvmag, y+0.1*v.*cos(y*pi/180)./uvmag, 'clip', clip);
% Vector angles in the Cartesian data-unit system of m_map:
%Ang = angle( (xsp-xs) + i*(ysp-ys) );
% ABove angle calc could fail when arrows were near boundaries. Replace with this.
% Now we calculate angles. Use a small offset, and keep clipping off to prevent odd things
% happening.
%    - RP 7/Jun/06
[xsp, ysp] = m_ll2xy([x x+0.00001*u./uvmag]',[y y+0.00001*v.*cos(y*pi/180)./uvmag]' , 'clip', 'off');
Ang = angle( diff(xsp)' + i*diff(ysp)' );

if ~isempty(key), Ang = 0; end


nvec = length(L);
Zero = zeros(nvec,1);
One  = ones(nvec,1);

% Normal arrow dimensions:
HL = Zero+HeadLength;
HW = Zero+HeadWidth;
W =  Zero+Width;

% Distinguish zero-length vectors from non-zero:
mm = (L < 100*eps);
i_zero = find(mm);
i_nonzero = find(~mm);
% Don't plot if length is zero.
if ~isempty(i_zero)
   HL(i_zero) = NaN; HW(i_zero) = NaN; W(i_zero) = NaN;
end

if ~isempty(i_nonzero)
   ii = i_nonzero;
   if HeadLength == 0,   %% square end; no arrowhead
      W(ii)  = Width;
      HW(ii) = Width;
      HL(ii) = Zero(ii);  % Thanks for Pierre Jaccard for this fix
   else
      % If the arrow length is less than the headlength,
      % omit the arrow shaft and just plot a head scaled
      % to the length.
      i_short = ii( find(L(ii) < HeadLength) );
      W(i_short)  = 0;
      HL(i_short) = L(i_short);
      HW(i_short) = HL(i_short) * (HeadWidth/HeadLength);
   end
end


% Just a change of variable names for historical reasons:
Y = ys;
X = xs;
% It is not clear whether non-zero elevations will be useful;
% that depends on whether the entire m_map system will work
% with 3-D views, as in making a surface plot of topography.
% Then one might want to have a current profile represented
% as a stack of vectors at appropriate heights above the topog.
Z = z(:);


Xzero = Zero;
if centered,
   Xzero = -L/2;
end


nV = 7*nvec;  % number of Vertices: 7 per arrow.
Vert = zeros(nvec,3);

Vert(1:7:nV,:) =  [Xzero,  W/2, Z ];       % back corner of shaft
Vert(2:7:nV,:) = Vert(1:7:nV,:) + ...
                 [(L-HL), Zero,  Zero];    % shaft-head junction
Vert(3:7:nV,:) = Vert(2:7:nV,:) + ...
                 [Zero, (HW-W)/2,  Zero];  % point of barb
Vert(4:7:nV,:) = [Xzero + L, Zero, Z];     % tip of arrow

%% This could be done more efficiently with fancier indexing, but
%% it is probably not a bottleneck, hence not worth the trouble.

% Reflect the top half to get the bottom half.
% First replicate:
Vert(5:7:nV,:) = Vert(3:7:nV,:);
Vert(6:7:nV,:) = Vert(2:7:nV,:);
Vert(7:7:nV,:) = Vert(1:7:nV,:);
% Then negate y to reflect:
Vert(5:7:nV,2) = -Vert(5:7:nV,2);
Vert(6:7:nV,2) = -Vert(6:7:nV,2);
Vert(7:7:nV,2) = -Vert(7:7:nV,2);

% Make an index array for operating on all vertices of each vector:
ii = (1:nvec);
ii = ii(ones(7,1),:);
ii = ii(:);

%% Rotate:
Vxy = exp(i*Ang(ii)).*(Vert(:,1) + i*Vert(:,2));

%% Translate:
Vxy = Vxy + X(ii) + i*Y(ii);

%%
Vert(:,1) = real(Vxy);
Vert(:,2) = imag(Vxy);


Faces = [1:nV].';            %Top
Faces = reshape(Faces,7,nvec).';

% Extremely narrow patches don't show up on the screen (although they seem
% to be printed OK) when EdgeColor is 'none', so when the arrows are all
% the same color, set the EdgeColor to be the same as FaceColor.
% Set clip off here so arrows are complete - RP 7/Jun/06

% Request to clip arrows at the plot edge.

if strcmp(edgeclip,'off'),
 hp = patch('Faces', Faces, 'Vertices', Vert, 'tag', 'm_vec','clipping','off');
else,
  [LG,LN]=m_xy2ll(reshape(Vert(:,1),7,nvec),reshape(Vert(:,2),7,nvec)); % Converts vertices in 7 point lines (i.e. columns) in lat/long
  [X,Y]=m_ll2xy(LG,LN  ,'clip','patch');                                % Converts back to x/y, but does clipping on for each column
  hp = patch('Faces', Faces, 'Vertices', [X(:) Y(:)], 'tag', 'm_vec','clipping','off');
end;

if ischar(c) || (size(c,1) == 1 && size(c,2) == 3),
   set(hp, 'EdgeColor', c, 'FaceColor', c, 'LineWidth', 0.1);
else
   set(hp, 'EdgeColor', 'none', 'FaceColor','Flat', ...
     'FaceVertexCdata', c);
end
if ~isempty(keyvars)
   set(hp, keyvars{:});
end

if ~isempty(key)
   ht = text(X(1), Y(1)-0.5*HW, Z(1), key, ...
      'color', c, 'horizontalalignment','left','verticalalignment','top', 'tag','m_vec', ...
      'clipping','off');
   set(hp,'clipping','off')
else
   ht = [];
end

if nargout==0
  clear hp ht
end;
