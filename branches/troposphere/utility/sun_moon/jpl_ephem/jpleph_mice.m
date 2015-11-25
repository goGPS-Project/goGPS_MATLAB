function rrd = jpleph_mice (et, ntarg, ncent)

% reads the jpl planetary ephemeris and gives the position and velocity
% of the point 'ntarg' with respect to point 'ncent' using MICE routines

% input

%   et    = TDB julian date at which interpolation is wanted

%   ntarg = integer number of 'target' point

%   ncent = integer number of center point

%   the numbering convention for 'ntarg' and 'ncent' is:

%        1 = mercury           8 = neptune
%        2 = venus             9 = pluto
%        3 = earth            10 = moon
%        4 = mars             11 = sun
%        5 = jupiter
%        6 = saturn
%        7 = uranus

% output

%   rrd = output 6-word array containing position and velocity
%         of point 'ntarg' relative to 'ncent'. the units are
%         determined by the value of km passed via global.

% global

%   iephem  = initialization flag (1 = initialize)
%   ephname = name of ephemeris binary data file (de421.bsp, etc.)
%   km      = state vector units flag (1 = km & km/sec, 0 = au & au/day)
%   au      = numerical value of astronomical unit (kilometers)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iephem ephname km au

if (iephem == 1)
    
    % load binary ephemeris data file
    
    cspice_furnsh(ephname);
    
    % reset initialization flag
    
    iephem = 0;
    
end

% set name of target body

switch ntarg
    
    case (1)
        
        targ = 'mercury';
        
    case (2)
        
        targ = 'venus';
        
    case (3)
        
        targ = 'earth';
        
    case (4)
        
        targ = 'mars';
        
    case(5)
        
        targ = 'jupiter';
        
    case (6)
        
        targ = 'saturn';
        
    case (7)
        
        targ = 'neptune';
        
    case (8)
        
        targ = 'uranus';
        
    case (9)
        
        targ = 'pluto';
        
    case (10)
        
        targ = 'moon';
        
    case (11)
        
        targ = 'sun';
        
end

% set name of central body

switch ncent
    
    case (1)
        
        obs = 'mercury';
        
    case (2)
        
        obs = 'venus';
        
    case (3)
        
        obs = 'earth';
        
    case (4)
        
        obs = 'mars';
        
    case (5)
        
        obs = 'jupiter';
        
    case (6)
        
        obs = 'saturn';
        
    case (7)
        
        obs = 'neptune';
        
    case (8)
        
        obs = 'uranus';
        
    case (9)
        
        obs = 'pluto';
        
    case (10)
        
        obs = 'moon';
        
    case (11)
        
        obs = 'sun';
        
end

% compute time, expressed as TDB seconds past J2000 TDB (2451545.0)

etime = 86400.0d0 * (et - 2451545.0d0);

% compute position and velocity vectors in eme2000 system (no corrections)

starg = mice_spkezr(targ, etime, 'J2000', 'NONE', obs);

% provide output in user-requested units

if (km == 1)
    
    % state is kilometers and kilometers/second
    
    rrd = starg.state;
    
else
    
    % state is au's and au's/day
    
    rrd(1:3) = starg.state(1:3) / au;
    
    rrd(4:6) = 86400.0 * starg.state(4:6) / au;
    
    rrd = rrd';
    
end





