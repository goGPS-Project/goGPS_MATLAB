function [ position, velocity, acceleration, error ] = planet_ephemeris( s, tjd, target, center )
% This function accesses the JPL planetary ephemeris
% to give the position and velocity of the target object with respect to the
% center object.
% 
% short int planet_ephemeris (double tjd[2], short int target,
%                             short int center,
% 
%                             double *position, double *velocity)
%
% ------------------------------------------------------------------------
% 
%    PURPOSE:
%       This function accesses the JPL planetary ephemeris to give the
%       position and velocity of the target object with respect to the
%       center object.
% 
%    REFERENCES:
%       Standish, E.M. and Newhall, X X (1988). "The JPL Export
%          Planetary Ephemeris"; JPL document dated 17 June 1988.
% 
%    INPUT
%    ARGUMENTS:
%       tjd[2] (double)
%          Two-element array containing the Julian date, which may be
%          split any way (although the first element is usually the
%          "integer" part, and the second element is the "fractional"
%          part).  Julian date is in the TDB or "T_eph" time scale.
%       target (short int)
%          Number of 'target' point.
%       center (short int)
%          Number of 'center' (origin) point.
%          The numbering convention for 'target' and'center' is:
%             1 = Mercury           8 = Neptune
%             2 = Venus             9 = Pluto
%             3 = Earth            10 = Moon
%             4 = Mars             11 = Sun
%             5 = Jupiter          12 = Solar system bary.
%             6 = Saturn           13 = Earth-Moon bary.
%             7 = Uranus           14 = Nutations (long int. and obliq.)
%                                  15 = Librations
%                                  16 = Lunar Euler angle rates
%                                  17 = TT-TDB
%             (If nutations are desired, set 'target' = 14;
%              'center' should be zero.)
% 
%    OUTPUT
%    ARGUMENTS:
%       *position (double)
%          Position vector array of target relative to center, measured
%          in AU or km.
%          Nutations returned in longitude (1) and obliquity (2) in radians
%          Librations returned in radians
%       *velocity (double)
%          Velocity vector array of target relative to center, measured
%          in AU/day or km/sec.
%          Nutations rate returned in longitude (1) and obliquity (2) in
%          radians/day or radians/sec
%          Librations rate returned in radians/day or radians/sec
%       *acceleration (double)
%          acceleration vector array of target relative to center, measured
%          in AU/day^2 or km/sec^2.
%          Nutations rate returned in longitude (1) and obliquity (2) in
%          radians/day^2 or radians/sec^2
%          Librations rate returned in radians/day^2 or radians/sec^2
% 
%    RETURNED
%    VALUE:
%       (short int)
%          0  ...everything OK.
%          1,2...error returned from State.
% 
%    GLOBALS
%    USED:
%       EM_RATIO          eph_manager.h
% 
%    FUNCTIONS
%    CALLED:
%       state             eph_manager.h
% 
%    VER./DATE/
%    PROGRAMMER:
%       V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
%       V1.1/07-93/WTH (USNO/AA): Update to C standards.
%       V2.0/07-98/WTH (USNO/AA): Modified for ease of use and linearity.
%       V3.0/11-06/JAB (USNO/AA): Allowed for use of input 'split' Julian
%                                 date for higher precision.
%       V3.1/11-07/WKP (USNO/AA): Updated prolog and error codes.
%       V3.1/12-07/WKP (USNO/AA): Removed unreferenced variables.
%       V3.2/10-10/WKP (USNO/AA): Renamed function to lowercase to
%                                 comply with coding standards.
% 
%    NOTES:
%       None.
% 
% ------------------------------------------------------------------------
  persistent zeroVector;
  if isempty(zeroVector)
    zeroVector = Ephem.newVector();
  end
  position = zeroVector;
  velocity = zeroVector;
  acceleration = zeroVector;

  if s.SS(3) == 0
    if s.output > 0
      fprintf(s.output,'planet_ephemeris: No starting JD in DE%d\n',...
       s.de_number);
    end
    error = 99;
    return;
  end

  error = 0;

  % Initialize 'jed' for 'state' and set up component count.
  jed = tjd;

  % Check for target point = center point.

  if target == center
    return;
  end
  
  %   Check for instances of target or center being Earth or Moon,
  %   and for target or center being the Earth-Moon barycenter.

  if (target == Ephem.EarthMoonBarycenter) || (center == Ephem.EarthMoonBarycenter) || ...
     (target == Ephem.Moon) || (center == Ephem.Moon) || ...
     (target == Ephem.Earth) || (center == Ephem.Earth)
    [pos_earthmoon,vel_earthmoon,acc_earthmoon,error] = s.state (jed, Ephem.EarthMoonBarycenterState);
    if error
      return;
    end
  end

  if (target == Ephem.Earth) || (center == Ephem.Earth) || ...
     (target == Ephem.Moon) || (center == Ephem.Moon)
    [pos_moon,vel_moon,acc_moon,error] = s.state (jed, Ephem.MoonState);
    if error
      return;
    end
  end
  
  % Check for cases of Earth as target and Moon as center or vice versa.

  if (target == Ephem.Earth) && (center == Ephem.Moon)
    % selenocentric Earth
    position = -pos_moon;
    velocity = -vel_moon;
    acceleration = -acc_moon;
    return;
  elseif (target == Ephem.Moon) && (center == Ephem.Earth)
    % geocentric Moon
    position = pos_moon;
    velocity = vel_moon;
    acceleration = acc_moon;
    return;
  end

  %   Make call to State for target object.
  switch target
    case Ephem.EarthMoonBarycenter
      targetState = Ephem.EarthMoonBarycenterState;
    case Ephem.Nutations
      targetState = Ephem.NutationsState;
    case Ephem.Librations
      targetState = Ephem.LibrationsState;
    case Ephem.LunarEulerRates
      targetState = Ephem.LunarEulerState;
    case Ephem.TTmTDB
      targetState = Ephem.TTState;
    otherwise
      targetState = target;
  end

  if target <= 0 || target == Ephem.SolarSystemBarycenter
    target_pos = zeroVector;
    target_vel = zeroVector;
    target_acc = zeroVector;
  elseif target == Ephem.EarthMoonBarycenter
    target_pos = pos_earthmoon;
    target_vel = vel_earthmoon;
    target_acc = acc_earthmoon;
  elseif target == Ephem.Earth
    target_pos = pos_earthmoon - (pos_moon / (1.0 + s.EM_RATIO));
    target_vel = vel_earthmoon - (vel_moon / (1.0 + s.EM_RATIO));
    target_acc = acc_earthmoon - (acc_moon / (1.0 + s.EM_RATIO));
  elseif target == Ephem.Moon
    target_pos = pos_earthmoon + pos_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
    target_vel = vel_earthmoon + vel_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
    target_acc = acc_earthmoon + acc_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
  elseif target >= Ephem.Nutations
    [~,nc] = size(s.IPT);
    if targetState <= nc && s.IPT(2,targetState) > 0
      [target_pos,target_vel,target_acc,error] = s.state (jed, targetState);
      if target == Ephem.TTmTDB
        target_pos(2) = 0.0;
        target_pos(3) = 0.0;
        target_vel(2) = 0.0;
        target_vel(3) = 0.0;
        target_acc(2) = 0.0;
        target_acc(3) = 0.0;
      elseif target == Ephem.Nutations
        target_pos(3) = 0.0;
        target_vel(3) = 0.0;
        target_acc(3) = 0.0;
      end
    else
      error = 98;
    end
    if error
      return;
    end
  else
    [target_pos,target_vel,target_acc,error] = s.state (jed, targetState);
    if error
      return;
    end
  end

  % Make call to State for center object.
  switch center
    case Ephem.EarthMoonBarycenter
      centerState = Ephem.EarthMoonBarycenterState;
    case Ephem.Nutations
      centerState = Ephem.NutationsState;
    case Ephem.Librations
      centerState = Ephem.LibrationsState;
    case Ephem.LunarEulerRates
      centerState = Ephem.LunarEulerState;
    case Ephem.TTmTDB
      centerState = Ephem.TTState;
    otherwise
      centerState = center;
  end

  if center <= 0 || center == Ephem.SolarSystemBarycenter

    % If the requested center is the Solar System barycenter,
    % then don't bother with a second call to State.

    center_pos = zeroVector;
    center_vel = zeroVector;
    center_acc = zeroVector;

  elseif center == Ephem.EarthMoonBarycenter
    
    % Center is Earth-Moon barycenter, which was already computed above.
    
    center_pos = pos_earthmoon;
    center_vel = vel_earthmoon;
    center_acc = acc_earthmoon;
    
  elseif center == Ephem.Earth
    
    % Center is Earth, which was already computed above.
    
    center_pos = pos_earthmoon - (pos_moon / (1.0 + s.EM_RATIO));
    center_vel = vel_earthmoon - (vel_moon / (1.0 + s.EM_RATIO));
    center_acc = acc_earthmoon - (acc_moon / (1.0 + s.EM_RATIO));
   
  elseif center == Ephem.Moon
    
    center_pos = pos_earthmoon + pos_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
    center_vel = vel_earthmoon + vel_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
    center_acc = acc_earthmoon + acc_moon * (s.EM_RATIO/ (1.0 + s.EM_RATIO));
    
  elseif center >= Ephem.Nutations
    [~,nc] = size(s.IPT);
    if centerState <= nc && s.IPT(2,centerState) > 0
      [center_pos,center_vel,center_acc,error] = s.state (jed, centerState);
      if center == Ephem.TTmTDB
        center_pos(2) = 0.0;
        center_pos(3) = 0.0;
        center_vel(2) = 0.0;
        center_vel(3) = 0.0;
        center_acc(2) = 0.0;
        center_acc(3) = 0.0;
      elseif center == Ephem.Nutations
        center_pos(3) = 0.0;
        center_vel(3) = 0.0;
        center_acc(3) = 0.0;
      end
    else
      error = 98;
    end
    if error
      return;
    end
  else
    [center_pos,center_vel,center_acc,error] = s.state (jed, centerState);
    if error
      return;
    end
  end

  % Compute position and velocity vectors.

  position = target_pos - center_pos;
  velocity = target_vel - center_vel;
  acceleration = target_acc - center_acc;

end

