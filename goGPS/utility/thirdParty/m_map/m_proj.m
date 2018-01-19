function m_proj(proj,varargin)
% M_PROJ  Initializes map projections info, putting the result into a structure
%
%         M_PROJ('get') tells you the current state
%         M_PROJ('set') gives you a list of all possibilities
%         M_PROJ('set','proj name') gives info about a projection in the 
%                                   'get' list.
%         M_PROJ('proj name','property',value,...) initializes a projection.
%
%
%         see also M_GRID, M_LL2XY, M_XY2LL.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 20/Sep/01 - Added support for other coordinate systems.
% 25/Feb/07 - Swapped "get" and "set" at lines 34 and 47 
%		to make it consistent with the help 
%		(and common Matlab style)
%	    - Added lines 62-70 & 74 
%		to harden against error when no proj is set
%             (fixes thanks to Lars Barring)

global MAP_PROJECTION MAP_VAR_LIST MAP_COORDS

% Get all the projections
projections=m_getproj;

if nargin==0, proj='usage'; end

proj=lower(proj);

switch proj

  case 'set'              % Print out their names
    if nargin==1
      disp(' ');
      disp('Available projections are:'); 
      for k=1:length(projections)
        disp(['     ' projections(k).name]);
      end
    else
      k=m_match(varargin{1},projections(:).name);
      eval(['X=' projections(k).routine '(''set'',projections(k).name);']);
      disp(X);
    end

  case 'get'              % Get the values of all set parameters
    if nargin==1
      if isempty(MAP_PROJECTION)
         disp('No map projection initialized');
         m_proj('usage');
      else
         k=m_match(MAP_PROJECTION.name,projections(:).name);
         eval(['X=' projections(k).routine '(''get'');']);
         disp('Current mapping parameters -');
         disp(X);
      end
    else
      if isempty(MAP_PROJECTION)          
        k=m_match(varargin{1},projections(:).name);
        eval(['X=' projections(k).routine '(''set'',projections(k).name);']);
        X=strvcat(X, ...
                  ' ', ...
                  '**** No projection is currently defined      ****', ...
                  ['**** USE "m_proj(''' varargin{1} ''',<see options above>)" ****']);
        disp(X);
      else
	k=m_match(varargin{1},projections(:).name);
	eval(['X=' projections(k).routine '(''get'');']);
	disp(X);
      end	
    end

  case 'usage'
    disp(' ');
    disp('Possible calling options are:');
    disp('  ''usage''                    - this list');
    disp('  ''set''                      - list of projections');
    disp('  ''set'',''projection''         - list of properties for projection');
    disp('  ''get''                      - get current mapping parameters (if defined)');
    disp('  ''projection'' <,properties> - initialize projection\n');

 otherwise                % If a valid name, give the usage.
    k=m_match(proj,projections(:).name);
    MAP_PROJECTION=projections(k);
        
    eval([ projections(k).routine '(''initialize'',projections(k).name,varargin{:});']);

    % With the projection store what coordinate system we are using to define it.
    if isempty(MAP_COORDS)
      m_coord('geographic');
    end  
    MAP_PROJECTION.coordsystem=MAP_COORDS;
    
    % Save some other stuff that otherwise seems to take a while to run
    MAP_PROJECTION.version=ver('matlab');
    if isempty(MAP_PROJECTION.version)
       MAP_PROJECTION.version=ver('octave');
    end
    if strcmp(MAP_PROJECTION.version.Name,'Octave')
       MAP_PROJECTION.IsOctave=true;
       MAP_PROJECTION.newgraphics=false;
       MAP_PROJECTION.LARGVAL=bitmax;
    else
       MAP_PROJECTION.IsOctave=false;
       if verLessThan('matlab','8.4')
           MAP_PROJECTION.newgraphics=false;
       else
           MAP_PROJECTION.newgraphics=true;
       end   
      % I use bitmax in various places as 'a large number', but
      % as of 2014b this has been renamed
      if verLessThan('matlab','8.3')
           MAP_PROJECTION.LARGVAL=bitmax;
       else
           MAP_PROJECTION.LARGVAL=flintmax;
       end
    end
    
 


end

%---------------------------------------------------------
function projections=m_getproj
% M_GETPROJ Gets a list of the different projection routines
%           and returns a structure containing both their
%           names and the formal name of the projection.
%           (used by M_PROJ).

% Rich Pawlowicz (rich@ocgy.ubc.ca) 9/May/1997
%
% 9/May/97 - fixed paths for Macs (thanks to Dave Enfield)
%
% 7/05/98 - VMS pathnames (thanks to Helmut Eggers)

% Get all the projections

lpath=which('m_proj');
fslashes=findstr(lpath,'/');
bslashes=findstr(lpath,'\');
colons=findstr(lpath,':');
closparantheses=findstr(lpath,']');
if ~isempty(fslashes)
  lpath=[ lpath(1:max(fslashes)) 'private/'];
elseif ~isempty(bslashes)
  lpath=[ lpath(1:max(bslashes)) 'private\'];
elseif ~isempty(closparantheses)       % for VMS computers only, others don't use ']' in filenames
  lpath=[ lpath(1:max(closparantheses)-1) '.private]'];
else
  lpath=[ lpath(1:max(colons)) 'private:'];
end

w=dir([lpath 'mp_*.m']);

if isempty(w) % Not installed correctly
  disp('**********************************************************');
  disp('* ERROR - Can''t find anything in a /private subdirectory *');
  disp('*         m_map probably unzipped incorrectly - please   *');
  disp('*         unpack again, preserving directory structure   *');
  disp('*                                                        *');
  disp('*         ...Abandoning m_proj now.                      *');
  error('**********************************************************');
end  
	
l=1;
projections=[];
for k=1:length(w)
 funname=w(k).name(1:(findstr(w(k).name,'.'))-1);
 projections(l).routine=funname;
 eval(['names= ' projections(l).routine '(''name'');']);
 for m=1:length(names)
   projections(l).routine=funname;
   projections(l).name=names{m};
   l=l+1;
 end
end


%----------------------------------------------------------
function match=m_match(arg,varargin)
% M_MATCH Tries to match input string with possible options

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

match=strmatch(lower(arg),cellstr(lower(char(varargin))));

if length(match)>1
  error(['Projection ''' arg ''' not a unique specification']);
elseif isempty(match)
  error(['Projection ''' arg ''' not recognized']);
end
