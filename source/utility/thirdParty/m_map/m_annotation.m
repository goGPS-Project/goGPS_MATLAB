function h = m_annotation( varargin )
% M_ANNOTATION Creates an annotation object
%   M_ANNOTATION Generally passes its arguments through to ANNOTATION,
%   bur accepts lon/lat arguments.
%
%   ANNOTATION(ANNOTATIONTYPE) creates a default annotation of type
%   ANNOTATIONTYPE in the current figure.  ANNOTATIONTYPE may be one of the
%   following:
%       'rectangle'
%       'ellipse'
%       'textbox'
%       'line'
%       'arrow'
%       'doublearrow' = two headed arrow
%       'textarrow' = arrow with text at tail end
%
%   ANNOTATION('rectangle',POSITION) creates a rectangle annotation at the
%   position specified in normalized figure units by the vector POSITION
%   ANNOTATION('ellipse',POSITION) creates an ellise annotation at the
%   position specified in normalized figure units by the vector POSITION
%   ANNOTATION('textbox',POSITION) creates a textbox annotation at the
%   position specified in normalized figure units by the vector POSITION
%
%   POSTITION = [LON_LEFT LAT_BOTTOM WIDTH HEIGHT] where WIDTH and HEIGHT
%   are normalized to the size of the axis.
%
%   ANNOTATION('line',LON,LAT) creates a line annotation with endpoints
%   specified in normalized figure coordinates by the vectors LON and LAT
%   ANNOTATION('arrow',LON,LAT) creates an arrow annotation with endpoints
%   specified in normalized figure coordinates by the vectors LON and LAT. 
%   LON(1) and LAT(1) specify the position of the tail end of the arrow 
%   and LON(2) and LAT(2) specify the position at the tip of the arrow head.
%   ANNOTATION('doublearrow',LON,LAT) creates a doublearrow annotation with
%   endpoints specified in normalized figure coordinates by the vectors LON
%   and LAT.
%   ANNOTATION('textarrow',LON,LAT) creates a textarrow annotation with
%   endpoints specified in normalized figure coordinates by the vectors LON
%   and LAT. LON(1) and LAT(1) specify the position of the tail end of the 
%   arrow and LON(2) and LAT(2) specify the position at the tip of the arrow 
%   head.
%
%   You should (optionally) call ORIENT first, and then WYSIWYG to see how
%   the printed version will look, because otherwise the location of the
%   annotation may not look right.
%
%   M_ANNOTATION(HANDLE,...) creates the annotation in the  AXES specified 
%   by HANDLE (note - ANNOTATION requires a FIGURE handle)
%
%   H=M_ANNOTATION(...) returns a handle to the annotation object
%
%   The arguments to ANNOTATION can be followed by parameter/value pairs to
%   specify additional properties of the annotation object. The X and Y or
%   POSITION arguments to ANNOTATION can be omitted entirely, and all
%   properties specified using parameter/value pairs.
%
%   Examples: rh=annotation('rectangle',[-123 48 .3 .3]); 
%             ah=annotation('arrow',[10 30],[-5 10],'Color','r');
%             th=annotation('textarrow',[-50 -30],[-89 -70],'String','ABC');
% 

% R. Pawlowicz Nov/2017
%
%

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

% Is first argument an axis handle?
if nargin>0 && length(varargin{1})==1 && ishandle(varargin{1}) 
   if  strcmp(get(varargin{1},'type'),'axes')
      useax=varargin{1};
      varargin(1)=[];
   else
      error(['map: ' mfilename ':invalidAxesHandle'],...
            ' First argument must be an axes handle ');
   end   
else
    useax=gca;
end
figH=get(useax,'parent');
 


annotype=[varargin{1} '   '];

switch lower(annotype(1:5)) 
   case {'recta','ellip','textb'}
       POS=varargin{2};
       LON=POS(1);
       LAT=POS(2);
       [X,Y]=m_ll2xy(LON,LAT,'clip','off');
       [fX,fY]=axistofig(useax,X,Y,MAP_VAR_LIST.xlims,MAP_VAR_LIST.ylims);
       varargin{2}=[fX fY POS(3) POS(4)];
       
   case {'line ','arrow','doubl','texta'}
       LON=varargin{2};
       LAT=varargin{3};
       [X,Y]=m_ll2xy(LON,LAT,'clip','off');
       [fX,fY]=axistofig(useax,X,Y,MAP_VAR_LIST.xlims,MAP_VAR_LIST.ylims);
       varargin{2}=fX;
       varargin{3}=fY;
   otherwise  % Just a pass through!    
end

 
h=annotation(figH,varargin{:});

set(h,'tag','m_annotation');

if nargout==0
    clear h
end

end

% Transform into normalized FIGURE coordinates
function [fX,fY]=axistofig(useax,X,Y,xlims,ylims)

axun = get(useax,'Units');
set(useax,'Units','normalized');
%%axpos = get(useax,'Position');  % This doesn't work if dataspectratio set
axpos =plotboxpos(useax);         %...but this does the right thing.

fX= (X-xlims(1))/diff(xlims)*axpos(3) + axpos(1);
fY= (Y-ylims(1))/diff(ylims)*axpos(4) + axpos(2);
 
%% Restore axes units
set(useax,'Units',axun);

end

% The following function came from MATLAB File exchange
% https://www.mathworks.com/matlabcentral/fileexchange/9615-plotboxpos
%
function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis limits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

set(hfig, 'currentaxes', currax);

end








