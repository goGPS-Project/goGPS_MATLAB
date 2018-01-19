function m_track(lon,lat,varargin)
% M_TRACK Draw a trackline on a map
%
%	M_TRACK draws navigation tracklines on a map. The trackline
%	itself is optionally annontated with tick marks, time labels, and
%	date labels.
%
%	M_TRACK (lon,lat) draws just a plain line. lon&lat are in decimal
%	degrees, +N, -S, +E, -W.
%
%	M_TRACK (lon,lat,navtimes) draws a line, with tick marks every
%	hour, time labels every four hours, and date labels every twelve
%	hours. navtimes is in MatLab "serial date numbers," representing the
%	number of days since January 1, 0000. By convention, ticks and 
%	time labels are drawn on the starboard side, dates on the port.
%
%	M_TRACK (lon,lat,navtime, 'string', property/value pairs) can be
%	used to change tick marks, time and date label properties. Allowable
%	combinations are:
%
%	'ticks'		tickmark interval in minutes, default=60
%	'times'		time label interval in minutes, default=240 (4hrs)
%	'dates'		date label interval in minutes, default=720 (12hrs)
%	'timef[ormat]'	format of time string, see datestr(), default=15
%	'datef[ormat]'	format of date string, see datestr(), default=2
%	'color'		color of track/labels, default is black
%	'linew[idth]'	width of trackline, default is current
%	'lines[tyle]'	style of line, default is solid line
%	'fonts[ize]'	size of font, default is current
%	'fontn[ame]'	name of font, default is current
%       'clip'		'on' or 'off' (clipping to map)
%	'orien[t]'	label orientation, 'true' or 'upright', default='true'
%
%	time labels need to be whole multiples of tick intervals. date
%	labels need to be whole multiples of time labels. using '0' for
%	any value will tick/label all points. using a negative number
%	will suppress tick/label.

% m_track.m, Peter Lemmond, peter@whoi.edu 13/Nov/98
%
% RP - 14/Nov/98 - corrected angle for first label, added tags, added CLIP
%                  options, made labels always face upwards
% PL - 04/Dec/98 - added orientation option, to allow labels to either
% 		   always follow heading ('true') or to always face upwards
%		   ('upright'). changed quite a bit internally to make
%		   faster and fix some bugs
%    - 23/Aug/01 - if numinputs = 2 an extra initialization statement
%                  was needed (thanks to Dan Lowen for this fix).


global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty (MAP_PROJECTION)
   disp ('No Map Projection initialized - call M_PROJ first!');
   return;
end

numinputs = nargin;			% save this

TICKS = 60;				% default of 60 minute ticks
TIMES = 240;				% default of 4 hour times, 240 mins
DATES = 720;				% default of 12 hour dates, 720 minutes
TIMEF = 15;				% default of HH:MM
DATEF = 2;				% default of mm/dd/yy
COLOR = 'k';				% default is black
LINES = '-';				% default is solid line
LINEW = get(gca,'linewidth');		% default is current width
FONTS = get(gca,'fontsize');		% default is current fontsize
FONTN = get(gca,'fontname');		% default is current fontname
CLIP  = 'on';				% default is to clip
ORIEN = 'true';				% default is always follow heading

MINSPERDAY = 1440;

% need at least lon & lat

if numinputs < 2
   disp ('Need at least lon & lat vectors!');
   return;
else
   l=length(lat);
   m=length(lon);
   if (l ~= m)
      disp ('long and lat vectors must be the same length');
      return;
   end
end

% check for time input. has to be the first varargin

if numinputs > 2
   if ischar(varargin{1})
      navTimes = [];
      k = 1;
   else
      navTimes=varargin{1};
      k = 2;
   end
else % numinputs = 2  % Extra case courtesy Dan Lowen.
        navTimes = [];
        k = 1;
end

% look at any remaining options

while k<length(varargin)
   optn=[lower(varargin{k}) '   '];
   switch optn(1:5)
      case 'ticks'
	 TICKS=varargin{k+1};
      case 'times'
	 TIMES=varargin{k+1};
      case 'dates'
	 DATES=varargin{k+1};
      case 'timef'
	 TIMEF=varargin{k+1};
      case 'datef'
	 DATEF=varargin{k+1};
      case 'color'
	 COLOR=varargin{k+1};
      case 'linew'
	 LINEW=varargin{k+1};
      case 'lines'
	 LINES=varargin{k+1};
      case 'fontn'
	 FONTN=varargin{k+1};
      case 'fonts'
	 FONTS=varargin{k+1};
      case 'clip '
	 CLIP=varargin{k+1};
      case 'orien'
	 ORIEN=varargin{k+1};
   end
   k=k+2;
end

% always want the line drawn at full resolution.

[x,y] = m_ll2xy(lon,lat,'clip',CLIP);

line(x,y,'clipping',CLIP,'linestyle',LINES,'linewidth',LINEW,'color',COLOR, ...
    'tag', 'm_track_line');


n=length(navTimes);
if ~n
   return				% no navtimes, all done
elseif (l ~= n)
   disp('long, lat, and navtimes vectors must be same length');
   return;
end

if TICKS < 0				% w/o ticks, nothing more
   return;
end

% which points to tick. will interpolate to exact, whole minute values,
% unless all requested

if (TICKS == 0)
   ty = y;
   tx = x;
   ttim = navTimes;
else
   tmp = TICKS/MINSPERDAY;
   i = tmp*ceil(min(navTimes)/tmp);
   j = tmp*floor(max(navTimes)/tmp);
   ttim = i:tmp:j;
   ty = interp1(navTimes,y,ttim,'linear');
   tx = interp1(navTimes,x,ttim,'linear');
end
nt = length(ttim);

% where do the time labels go?

if TIMES < 0				
   ltim=zeros(1,nt);			% no time lables
elseif TIMES == 0
   ltim=ones(1,nt);			% time label every tick
else
   tmp=TIMES/MINSPERDAY;
   ltim=~(mod(ttim,tmp));
end

% and date labels?

if DATES < 0
   dtim=zeros(1,nt);			% no date labels
elseif DATES == 0
   dtim=ltim;				% date label every time label
else
   tmp=DATES/MINSPERDAY;
   dtim=~(mod(ttim,tmp));
end

% since the 'rotation' attribute must be a scalar, have to loop thru
% and plot everything one at a time

for i=1:nt

   % angle for plotting. ticks are always perpendicular. will use the
   % point before and after this one to compute the angle
      
   dy  = ty(min(nt,i+1)) - ty(max(1,i-1));
   dx  = tx(min(nt,i+1)) - tx(max(1,i-1));
   angle = atan2(dy,dx)*180/pi - 90;
   
   % go ahead and tick here
   
   text(tx(i),ty(i),'-', 'verticalalignment', 'middle', 'horizontalalignment', 'left', ...
       'color', COLOR, 'fontsize', FONTS, 'fontname', FONTN, ...
       'clipping',CLIP,'rotation', angle,'tag','m_track_tick');

   % maybe time label here?
   
   if ltim(i)
      
      % make label at desired orientation
      
      if ((abs(angle) < 90) || ~strcmp(ORIEN,'upright'))
	 leadstr = ' ';
	 tailstr = '';
	 ang_off = 0;
	 tmlab   = 'left';
	 datlab  = 'right';
      else
	 leadstr ='';
	 tailstr =' ';
	 ang_off =180;
	 tmlab   ='right';
	 datlab  ='left';
      end
      
      % go ahead and plot the time
      
      text(tx(i),ty(i), ...
	  [leadstr leadstr datestr(ttim(i),TIMEF) tailstr tailstr], ...
	  'color', COLOR, 'fontsize', FONTS, 'fontname', FONTN, ...
	  'verticalalignment', 'middle', 'horizontalalignment', tmlab, 'rotation', ...
	  angle+ang_off, 'clipping', CLIP, 'tag','m_track_time');


      % maybe date label here?
      
      if dtim(i)

	 text(tx(i), ty(i), ...
	     [tailstr datestr(ttim(i),DATEF) leadstr], ...
	     'color', COLOR, 'fontsize', FONTS, 'fontname', FONTN, ...
	     'verticalalignment', 'middle', 'horizontalalignment', datlab, 'rotation', ...
	     angle-ang_off, 'clipping', CLIP, 'tag','m_track_date');
      
      end				% end of date labelling
   end					% end of time labelling
end					% end of ticking

