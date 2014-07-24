function [h] = m_legend(varargin)
% M_LEGEND Add a legend to a map. This function does not have the
%         full functionality of LEGEND but avoids some of the problems
%         encountered using that function in M_MAP. Feel free to
%         add code!
%
%    [legend_handle] = M_LEGEND(HANDLES,String1,String2,String3,...)
%
%    puts a legend on the plot containing the handles in the vector HANDLES 
%    using the specified strings as labels for the corresponding handles.
%    By holding the mouse button down, the legend can be moved around the
%    plot.

% Original Author:: Deirdre Byrne, dbyrne@umeoce.maine.edu 00/06/23

% This software is provided "as is" without warranty of any kind.

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if nargin < 2;
  help m_legend
  return
end

handles = varargin{1};
handles = handles(:);
varargin = varargin(:);
s = size(varargin,1);
if s-1 ~= size(handles,1)
  disp('must have same number of legends as handles')
  return
end

varargin = varargin(2:s);

% get axis location
ax=gca;
units = get(ax,'units');
p = get(ax,'pos');
c = get(ax,'color');
if strcmp(c,'none')
  c = get(gcf,'color');
end
if strcmp(c,'none')
  c = [0.702 0.702 0.702];
end

% calculate initial size for legend
xlen = (p(3))/5;
ylen = (p(4))/7;

% calculate position for legend
% default is to put it in bottom right corner
axpos(1) = p(1)+p(3) - xlen;
axpos(2) = p(2);
axpos(3:4) = [xlen ylen];

ButtonDownFcn = 'moveaxis';
Tag = 'legend';

h=axes('units',units,'position',axpos, ...
    'visible','on', ...
    'color',c, ...
    'xcolor','k', ...
    'ycolor','k', ...
    'tag',Tag, ...
    'Buttondownfcn',ButtonDownFcn, ...
    'xlim',[0 1], ...
    'xtick',[], ...
    'ylim',[0 1], ...
    'ytick',[], ...
    'nextplot','add', ...
    'box','on');

lh = length(handles);

% set width based on legend strings
axw = 4;
for i = 1:lh
  l = length(varargin{i});
  axw = max(axw,l);
end
axw = max(2*axw,axw+10);

sf = 1;
set(h,'units','characters');
axpos = get(h,'pos');
axh = (sf*lh+1);
set(h,'position',[(axpos(1)+axpos(3) - axw) axpos(2) axw axh])

for i = 1:lh
  xp = [0.05 0.2 0.35];
  iy = (i*sf-1)*0.8/lh + 0.25;
  yp = [iy iy iy];
  if strcmp(get(handles(i),'type'),'line')
    mark = get(handles(i),'marker');
    msize = get(handles(i),'markersize');
    linest = get(handles(i),'linestyle');
    linew = get(handles(i),'linewidth');
    clr = get(handles(i),'color');
    if ~strcmp(linest,'none')
      l1(i) = plot(xp,yp,'marker','none', ...
	  'color',clr,'linewidth',linew);
    end
    l2(i) = plot(xp(2),yp(2),'marker',mark, ...
	'linewidth',linew, 'markersize', msize, ...
	'color',clr);
    t(i) = text(0.5, iy, varargin{i});
  end
end


if nargout == 0
  clear h
end

