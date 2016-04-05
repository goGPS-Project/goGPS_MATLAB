function [ target_pos, target_vel, target_acc, error ] = state( s, jed, target )
% position, velocity, and acceleration of target
% 
% short int state (double *jed, short int target,
% 
%                  double *target_pos, double *target_vel)
% 
% ------------------------------------------------------------------------
% 
%    PURPOSE:
%       This function reads and interpolates the JPL planetary
%       ephemeris file.
% 
%    REFERENCES:
%       Standish, E.M. and Newhall, X X (1988). "The JPL Export
%          Planetary Ephemeris"; JPL document dated 17 June 1988.
% 
%    INPUT
%    ARGUMENTS:
%       *jed (double)
%          2-element Julian date (TDB) at which interpolation is wanted.
%          Any combination of jed[0]+jed[1] which falls within the time
%          span on the file is a permissible epoch.  See Note 1 below.
%       target (short int)
%          The requested body to get data for from the ephemeris file.
%          The designation of the astronomical bodies is:
%                  1 = Mercury                    7 = Uranus
%                  2 = Venus                      8 = Neptune
%                  3 = Earth-Moon barycenter      9 = Pluto
%                  4 = Mars                      10 = geocentric Moon
%                  5 = Jupiter                   11 = Sun
%                  6 = Saturn                    12 = Nutations
%                                                13 = Lunar Librations
%                                                14 = Lunar Euler angle rates
%                                                15 = TT-TDB
% 
%    OUTPUT
%    ARGUMENTS:
%       *target_pos (double)
%          The barycentric position vector array of the requested object,
%          in AU.
%          (If target object is the Moon, then the vector is geocentric.)
%       *target_vel (double)
%          The barycentric velocity vector array of the requested object,
%          in AU/Day.
% 
%          Both vectors are referenced to the Earth mean equator and
%          equinox of epoch.
% 
%    RETURNED
%    VALUE:
%       (short int)
%          0...everything OK.
%          1...error reading ephemeris file.
%          2...epoch out of range.
% 
%    GLOBALS
%    USED:
%       KM                eph_manager.h
%       EPHFILE           eph_manager.h
%       IPT               eph_manager.h
%       BUFFER            eph_manager.h
%       NRL               eph_manager.h
%       RECORD_LENGTH     eph_manager.h
%       SS                eph_manager.h
%       JPLAU             eph_manager.h
% 
%    FUNCTIONS
%    CALLED:
%       split             eph_manager.h
%       fseek             stdio.h
%       fread             stdio.h
%       interpolate       eph_manager.h
%       ephem_close       eph_manager.h
% 
%    VER./DATE/
%    PROGRAMMER:
%       V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
%       V1.1/07-93/WTH (USNO/AA): Update to C standards.
%       V2.0/07-98/WTH (USNO/AA): Modify to make position and velocity
%                                 two distinct vector arrays.  Routine set
%                                 to compute one state per call.
%       V2.1/11-07/WKP (USNO/AA): Updated prolog.
%       V2.2/10-10/WKP (USNO/AA): Renamed function to lowercase to
%                                 comply with coding standards.
% 
%    NOTES:
%       1. For ease in programming, the user may put the entire epoch in
%          jed[0] and set jed[1] = 0. For maximum interpolation accuracy,
%          set jed[0] = the most recent midnight at or before
%          interpolation epoch, and set jed[1] = fractional part of a day
%          elapsed between jed[0] and epoch. As an alternative, it may
%          prove convenient to set jed[0] = some fixed epoch, such as
%          start of the integration and jed[1] = elapsed interval between
%          then and epoch.
% 
% ------------------------------------------------------------------------
   
  %   Check epoch.
  jd = zeros(1,4);
  
  if length(jed) == 2
    js = double(jed(1)) - 0.5;
    [jd(1),jd(2)] = Ephem.split (js);
    [sp1,sp2] = Ephem.split (double(jed(2)));
    jd(2) = jd(2) + sp2;
    [jd(3),jd(4)] = Ephem.split (jd(2));
    jd(1) = jd(1) + sp1 + 0.5 + jd(3);
  else
    js = double(jed) - 0.5;
    [jd(1),jd(2)] = Ephem.split (js);
    %jd(2) = jd(2) + 0.0;
    [jd(3),jd(4)] = Ephem.split (jd(2));
    jd(1) = jd(1) + 0.5 + jd(3);
  end

  %   Return error code if date is out of range.

  if (jd(1) < s.SS(1)) || ...
     ((jd(1) + jd(4)) > s.SS(2))
    if jd(1) < s.SS(1)
      if s.output > 0
        fprintf(s.output,'state: JD %f days too early %g < %g\n',...
          -jd(1)+s.SS(1),jd(1),s.SS(1));
      end
    end
    if (jd(1) + jd(4)) > s.SS(2)
      if s.output > 0
        fprintf(s.output,'state: JD %f days too late %g > %g\n',...
          jd(1) + jd(4)-s.SS(2),jd(1) + jd(4),s.SS(2));
      end
    end
    target_pos = Ephem.newVector();
    target_vel = Ephem.newVector();
    target_acc = Ephem.newVector();
    error = 2;
    return;
  end

  %   Calculate record number and relative time interval.

  nr = int32(floor ((jd(1) - s.SS(1)) / s.SS(3)));
  if jd(1) == s.SS(2)
    nr = nr - 1;
  end
  t = double (nr) * s.SS(3) + s.SS(1);
  t = ((jd(1) - t) + jd(4)) / s.SS(3);
  if t < 0 || t > 1
    if s.output > 0
      fprintf(s.output,'state: time out of range [0,1] %g\n',t);
    end
    %error = 99;
  end
   
  %   Set units based on value of the 'KM' flag.

  if s.KM == 1
    daysfac =  1.0 / 86400.0 / s.SS(3); % days/second
    aufac = 1.0; % km or radians
  else
    daysfac = 1.0 / s.SS(3); % days
    if target < Ephem.NutationsState
      aufac = 1.0 / s.JPLAU; % AU/km
    else
      aufac = 1.0; % radians
    end
  end

  %   Check and interpolate for requested body.

  [target_pos,target_vel,target_acc] = s.interpolate (nr, target, t);

  target_pos = target_pos * aufac; % km or AU
  target_vel = target_vel * aufac * daysfac; % km/s or AU/day
  target_acc = target_acc * aufac * daysfac * daysfac; % km/s^2 or AU/day^2
  error = 0;

end

