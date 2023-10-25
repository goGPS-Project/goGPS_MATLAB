function [GI,xlm,ylm]=m_image(lon,lat,C,varargin)
% M_IMAGE An image replacement for data regularly gridded in LON/LAT
%   M_IMAGE(LON,LAT,C) transforms a dataset regularly gridded in geographic
%   coordinates into one regularly gridded in map coordinates and displays
%   this. LON and LAT can be two-element vectors, in which case it is
%   assumed that the data is C is evenly spaced between these limits, or
%   they can be monotonic vectors (i.e. no jumps from 180 to -180 in
%   increasing longitudes).  C should be a true-colour uint8 image of
%   size NxMx3 (in which case nearest-neighbour interpolation is used to 
%   convert to map coordinates) or a double precision NxM image in which 
%   case linear interpolation is used.
%
%   M_IMAGE(...,'resolution',N) where is a scalar or N=[NX NY] interpolates 
%   the image to N or NX/NY pixels in the horizontal and vertical
%   directions (default is 1000)
%
%   IM=M_IMAGE(...) returns the image object created.
%
%   [IM,X,Y]=M_IMAGE(...) does not create an image object, instead it
%   returns the interpolated object with map coordinates X/Y for further 
%   processing (e.g., a call to M_SHADEDRELIEF(X,Y,IM,'coords','map') )
%
%   See alse m_shadedrelief, m_pcolor

% R. Pawlowicz (rich@eos.ubc.ca) Jan/2018
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% Feb/2019 - added handling for straight image call using real elevations
%            (with alphamapping for out-of-map areas)
 
global MAP_PROJECTION MAP_VAR_LIST
% Have to have initialized a map first
if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

N=[1000 999];  % Default size of grid


k=1;
while k<length(varargin)
    switch lower(varargin{k}(1:3))
        case 'res'
            N=varargin{k+1};
            if length(N)==1
                N=[N N];
            end
            varargin([k k+1])=[];
            k=k-1;
            
     end
     k=k+1;
end


if ~isvector(lon)       % Can't be a matrix
    error('Input LON must be a VECTOR');
end
if ~isvector(lat)
    error('Input LAT must be a VECTOR');
end

% Even spacing if 2-element vectors
if length(lon)==2
    lon=linspace(lon(1),lon(2),size(C,2));
end
if length(lat)==2
    lat=linspace(lat(1),lat(2),size(C,1));
end

% ndgrid needs monotonically increasing parameters.
if lat(end)<lat(1)
   lat=lat(end:-1:1);
   C=C(end:-1:1,:,:);
end
if lon(end)<lon(1)
   lon=lon(end:-1:1);
   C=C(:,end:-1:1,:);
end
[GC1,GC2]=ndgrid(lat,lon);

% Get a grid to evenly cover the map area 
xlm=linspace(MAP_VAR_LIST.xlims(1),MAP_VAR_LIST.xlims(2),N(1));
ylm=linspace(MAP_VAR_LIST.ylims(1),MAP_VAR_LIST.ylims(2),N(2));
[XX,YY]=meshgrid(xlm,ylm);
% ...and get the lat/long for each point.
[HLG,HLT]=m_xy2ll(XX,YY);



% Try to fix any +/-360 problems.
if strcmp(MAP_PROJECTION.routine,'mp_azim') | strcmp(MAP_PROJECTION.routine,'mp_omerc') 
    ii=HLG<MAP_VAR_LIST.longs(1) & HLG+360<MAP_VAR_LIST.longs(2);
    if any(ii(:)) HLG(ii)=HLG(ii)+360; end
    ii=HLG>MAP_VAR_LIST.longs(2) & HLG-360>MAP_VAR_LIST.longs(1);
    if any(ii(:)) HLG(ii)=HLG(ii)-360; end
end

% Find pixels outside the limits of the actual map (if the boundary
% isn't a rectangle)
if strcmp(MAP_VAR_LIST.rectbox,'off')
    [I,J]=find(HLT<MAP_VAR_LIST.lats(1) | HLT>MAP_VAR_LIST.lats(2) | HLG<MAP_VAR_LIST.longs(1) | HLG>MAP_VAR_LIST.longs(2));
elseif strcmp(MAP_VAR_LIST.rectbox,'circle')
    R=(XX.^2 +YY.^2);
    [I,J]=find(R>MAP_VAR_LIST.rhomax.^2);
else
    I=[];J=[];
end

% Find pixels outside the limits of the input data
 [I1,J1]=find(HLT < lat(1) | HLT>lat(end) | HLG<lon(1) | HLG>lon(end) );
 
 
% This is for true-colour images

if strcmp(class(C),'uint8') && length(size(C))==3
    
    backcolor=uint8(get(gcf,'color')*255);
    
    GI=zeros(length(ylm),length(xlm),3,'uint8');
    for k=1:3
       F=griddedInterpolant(GC1,GC2,single(C(:,:,k)),'nearest');  % Set up interpolant

       GI(:,:,k)=uint8(F(HLT,HLG));                               % Interpolate each color

       if any(I1)                       % if some pixels are gridded area
           IJ=sub2ind(size(GI),I1,J1,repmat(k,size(I1)));   
           GI(IJ)=uint8(255);           % Set them to white
       end
       if any(I)                          % if some pixels are outside the map area
           IJ=sub2ind(size(GI),I,J,repmat(k,size(I)));   
           GI(IJ)=backcolor(k);           % Set them to the background colour
       end

    end

elseif strcmp(class(C),'double') && length(size(C))==2

    GI=zeros(length(ylm),length(xlm));
   
    F=griddedInterpolant(GC1,GC2,C,'linear','none');  % Set up interpolant

    GI=F(HLT,HLG);                     % Interpolate each color

    if any(I)                          % if some pixels are outside the map area 
        IJ=sub2ind(size(GI),I,J);   
        GI(IJ)=NaN;           % Set them to NaN
    end
   
    
end

if nargout<=1
   if strcmp(class(C),'uint8') && length(size(C))==3
     if any(I)  % Points outside the map
        alphadata=  ones(size(GI,1),size(GI,2),'logical');
        IJ=sub2ind(size(alphadata),I,J);   
        alphadata(IJ)=0;
        GI=image('XData',xlm,'YData',ylm,'CData',GI,'alphadata',alphadata,varargin{:},'tag','m_image');         
     else
        GI=image('XData',xlm,'YData',ylm,'CData',GI,varargin{:},'tag','m_image');
     end
   elseif strcmp(class(C),'double') && length(size(C))==2
     GI=image('XData',xlm,'YData',ylm,'CData',GI,varargin{:},'tag','m_image',...
               'CDataMapping','scaled','alphadata',~isnan(GI));
   end
end


