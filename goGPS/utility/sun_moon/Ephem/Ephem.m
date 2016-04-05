classdef Ephem < handle
% Class for opening JPL ephemeris files and calculating raw positions
% Copyright(c) 2014 by Jonathan Kipling Knight (<a href="matlab:
% web('mailto:drkipknight@aol.com')">drkipknight@aol.com</a>)
% DEPENDENCIES: This class requires that ephemeris files are available online by
% ftp or have been previously downloaded and placed in the Matlab path.
  properties (Constant)
    % body numbers for use with planet_ephemeris
    
    Mercury               = 1; % Mercury body number
    Venus                 = 2; % Venus body number
    Earth                 = 3; % Earth body number
    Mars                  = 4; % Mars body number
    Jupiter               = 5; % Jupiter body number
    Saturn                = 6; % Saturn body number
    Uranus                = 7; % Uranus body number
    Neptune               = 8; % Neptune body number
    Pluto                 = 9; % Pluto body number
    Moon                  = 10; % Geocentric Moon body number
    Sun                   = 11; % Sun body number
    SolarSystemBarycenter = 12; % Solarsystem Barycenter body number
    EarthMoonBarycenter   = 13; % Earth-Moon Barycenter body number
    Nutations             = 14; % nutations body number
    Librations            = 15; % librations body number
    LunarEulerRates       = 16; % lunar euler angel rates body number
    TTmTDB                = 17; % TT-TDB body number
  end
  properties (Constant)
    % state numbers for use with state or interpolate
    
    MercuryState               = 1; % Mercury state number
    VenusState                 = 2; % Venus state number
    EarthMoonBarycenterState   = 3; % Earth-Moon Barycenter state number
    MarsState                  = 4; % Mars state number
    JupiterState               = 5; % Jupiter state number
    SaturnState                = 6; % Saturn state number
    UranusState                = 7; % Uranus state number
    NeptuneState               = 8; % Neptune state number
    PlutoState                 = 9; % Pluto state number
    MoonState                  = 10; % Geocentric Moon state number
    SunState                   = 11; % Sun state number
    NutationsState             = 12; % nutations state number
    LibrationsState            = 13; % librations state number
    LunarEulerState            = 14; % lunar euler angle rates state number
    TTState                    = 15; % TT-TDB state number
  end
  properties (Constant)
    % list of available ASCII Development Ephemerides
    de_numbers = {'102','200','202','403','405','406','410','413','414','418','421','422','423','424','430','430t','431'};
    % list of minimum years for ASCII Development Ephemerides
    minyears = [-1400, 1600, 1900, 1600, 1600,-3000, 1960, 1900, 1600, 1900, 1900,-3000, 1800,-3000, 1550, 1550,-13000];
    % list of increments between files for ASCII Development Ephemerides
    incyears = [  300,   20,   50,  100,   20,  100,   20,   25,  100,  150,  150,  100,   50,  100,  100,  100,  1000];
    % list of maximum years for ASCII Development Ephemerides
    maxyears = [ 2800, 2160, 2050, 2100, 2200, 2900, 2000, 2025, 2100, 1900, 2050, 2900, 2150, 2900, 2550, 2550, 16000];
  end
  properties (Constant)
    TTLsize = 84; % max character length of the 3 TTL lines in the binary ephemeris
    CNAMsize = 6; % max character length of the Constant Name in the binary ephemeris
    OLDCONMAX = 400; % old maximum for number of constants
    JPLftpSite = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii'; % site for ASCII files
  end
  properties
    de_number = 405; % DE number of data stored here
    % flag defining physical units of the output states.
    %          = 1, km and km/sec
    %          = 0, au and au/day
    % files are native km and km/sec
    KM = 1;

    SS = [0.0,0.0,0.0]; % JD of first record, JD of last record, days between records
    JPLAU = 0.0; % definition of the astronomical unit (km)
    EM_RATIO = 0.0; % ratio of Earth's mass to Moon's mass
    RECORD_LENGTH = 0; % number of bytes in record
    IPT = zeros(3,Ephem.TTState); % table of indices into record
    scan = []; % linear array of values in records
    indexInc = 0; % increment between records in scan
    indexOffset = 1; % offset for buffer
    numRecs = 0; % number of records
    jd_begin = 0.0; % Julian Date of first record
    jd_end = 0.0; % Julian Date of last record
    jd_inc = 0.0; % days between records
    masses = zeros(Ephem.Mercury,Ephem.EarthMoonBarycenter); % G*mass (AU^3/Day^2) of body number 1-13
    output = 1; % message output
  end
  methods
    function s = Ephem( denum )
      % constructor
      if nargin < 1
        denum = s.de_number;
      end
      if isnumeric(denum)
        denum = sprintf('%d',denum);
      end
      s = s.setDEConstants( denum );
    end
    [ position, velocity, acceleration, error ] = planet_ephemeris( s, tjd, target, center );
    [ target_pos, target_vel, target_acc, error ] = state( s, jed, target );
    [ position, velocity, acceleration ] = interpolate( s, nrl, target, t );
    [ s, error ] = openDEasc( s, denum, fn, dontSave );
    [ s, error ] = openDEeph( s, denum, fn, dontSave );
    [ s, error ] = setDEConstants( s, denum );
    function [ rtn ] = hasNutations( s )
      % check if ephemeris file has nutations
      rtn = (s.IPT(2,Ephem.NutationsState) > 0);
    end
    function [ rtn ] = hasLibrations( s )
      % check if ephemeris file has librations
      [~,nc] = size(s.IPT);
      if nc >= Ephem.LibrationsState && s.IPT(2,Ephem.LibrationsState) > 0
        rtn = true;
      else
        rtn = false;
      end
    end
    function [ rtn ] = hasTTmTDB( s )
      % check if ephemeris file has TT-TDB
      [~,nc] = size(s.IPT);
      if nc >= Ephem.TTState && s.IPT(2,Ephem.TTState) > 0
        rtn = true;
      else
        rtn = false;
      end
    end
    function [ rtn ] = hasLunarEulerAngleRates( s )
      % check if ephemeris file has lunar Euler angle rates
      [~,nc] = size(s.IPT);
      if nc >= Ephem.LunarEulerState && s.IPT(2,Ephem.LunarEulerState) > 0
        rtn = true;
      else
        rtn = false;
      end
    end
  end
  methods (Static)
    function vector = newVector()
      % create a new position (row vector)
      vector = zeros(1,3);
    end
    [ fr1, fr2 ] = split( tt );
    [ error ] = asc2eph( denum, fn, OUTFILE );
    [ st, error ] = readHeader( denum );
    [ scan, error ] = read_asc( denum, fn, dontSave );
    [ fn ] = asc_name( denum, year );
    [ fn ] = eph_name( denum, year );
    function [ NCOEFF ] = numCoeff( IPT )
      % calculate number of coefficients in record based on IPT matrix
      KMX = 0;
      KHI = 0;
      [~,nc] = size(IPT);
      for I = 1:nc
        if (IPT(1,I) > KMX) && (IPT(2,I) > 0)
          KMX = IPT(1,I);
          KHI = I;
        end
      end
      if KHI == 0
        NCOEFF = 0;
        return;
      end

      if (KHI == Ephem.TTState)
        ND = 1;
      elseif (KHI == Ephem.NutationsState)
        ND = 2;
      else
        ND = 3;
      end

      NCOEFF = IPT(1,KHI)+ND*IPT(2,KHI)*IPT(3,KHI)-1;
    end
    function [rtn] = createAllEphFiles( denum, regen )
      % download and convert all ASCII ephemeris files to binary
      if nargin < 2
        regen = false;
      end
      if isnumeric(denum)
        denum = sprintf('%d',denum);
      end
      index = {};
      for i=1:length(Ephem.de_numbers)
        if strcmp(Ephem.de_numbers{i},denum)
          index = i;
          break;
        end
      end
      if isempty(index)
        rtn = 1;
        return;
      end
      for yr=Ephem.minyears(index):Ephem.incyears(index):Ephem.maxyears(index)
        infile = Ephem.asc_name(denum,yr);
        outfile = Ephem.eph_name(denum,yr);
        if regen || exist(outfile,'file') ~= 2
          rtn = Ephem.asc2eph(denum,infile,outfile);
        end
      end
    end
    rtn = TESTZERODAY( denum, testCase );
    rtn = TESTPOINTS( denum, testCase );
  end
%   methods (Test)
%     function test102(testCase)
%       % test DE102
%       Ephem.TESTZERODAY(102, testCase);
%       Ephem.TESTPOINTS(102, testCase);
%     end
%     function test200(testCase)
%       % test DE200
%       %Ephem.TESTZERODAY(200, testCase);
%       Ephem.TESTPOINTS(200, testCase);
%     end
%     function test405(testCase)
%       % test DE405
%       Ephem.TESTZERODAY(405, testCase);
%       Ephem.TESTPOINTS(405, testCase);
%     end
%     function test430(testCase)
%       % test DE430
%       Ephem.TESTZERODAY(430, testCase);
%       Ephem.TESTPOINTS(430, testCase);
%     end
%     function test430t(testCase)
%       % test DE430t
%       Ephem.TESTZERODAY('430t', testCase);
%       Ephem.TESTPOINTS('430t', testCase);
%     end
%     function test431(testCase)
%       % test DE431
%       Ephem.TESTZERODAY(431, testCase);
%       Ephem.TESTPOINTS(431, testCase);
%     end
%   end

end

