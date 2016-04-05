function [ position, velocity, acceleration ] = interpolate( s, nrl, target, t )
% This function differentiates and interpolates a set of Chebyshev coefficients to give position and velocity.
% 
% void interpolate (double *buf, double *t, long int ncf, long int na,
% 
%                   double *position, double *velocity)
% 
% ------------------------------------------------------------------------
% 
%    PURPOSE:
%       This function differentiates and interpolates a set of
%       Chebyshev coefficients to give position and velocity.
% 
%    REFERENCES:
%       Standish, E.M. and Newhall, X X (1988). "The JPL Export
%          Planetary Ephemeris"; JPL document dated 17 June 1988.
% 
%    INPUT
%    ARGUMENTS:
%       *buf (double)
%          Array of Chebyshev coefficients of position.
%       *t (double)
%          t[0] is fractional time interval covered by coefficients at
%          which interpolation is desired (0 <= t[0] <= 1).
%          t[1] is length of whole interval in input time units.
%       ncf (long int)
%          Number of coefficients per component.
%       na (long int)
%          Number of sets of coefficients in full array
%          (i.e., number of sub-intervals in full interval).
% 
%    OUTPUT
%    ARGUMENTS:
%       *position (double)
%          Position array of requested object.
%       *velocity (double)
%          Velocity array of requested object.
% 
%    RETURNED
%    VALUE:
%       None.
% 
%    GLOBALS
%    USED:
%       NP                eph_manager.h
%       NV                eph_manager.h
%       PC                eph_manager.h
%       VC                eph_manager.h
%       TWOT              eph_manager.h
% 
%    FUNCTIONS
%    CALLED:
%       fmod              math.h
% 
%    VER./DATE/
%    PROGRAMMER:
%       V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
%       V1.1/07-93/WTH (USNO/AA): Update to C standards.
%       V1.2/07-98/WTH (USNO/AA): Modify to make position and velocity
%                                 two distinct vector arrays.
%       V1.3/11-07/WKP (USNO/AA): Updated prolog.
%       V1.4/12-07/WKP (USNO/AA): Changed ncf and na arguments from short
%                                 int to long int.
%       V1.5/10-10/WKP (USNO/AA): Renamed function to lowercase to
%                                 comply with coding standards.
% 
%    NOTES:
%       None.
% 
% ------------------------------------------------------------------------
  
  persistent PC VC AC NP TWOT;
  if isempty(PC)
    PC = zeros(1,15); % ncf <= 15
    VC = zeros(1,15); % ncf <= 15
    AC = zeros(1,15); % ncf <= 15
    TWOT = 0.0;
    PC(1) = 1.0; % T_0(x) = 1
    PC(2) = 0.0; % T_1(x) = x
    VC(1) = 0.0; % T'_0(x) = 0
    VC(2) = 1.0; % T'_1(x) = 1
    AC(1) = 0.0; % T"_0(x) = 0
    AC(2) = 0.0; % T"_1(x) = 0
    AC(3) = 4.0; % T"_2(x) = 4
    NP = 2;
  end
  if target == Ephem.TTState
    NDIM = 1;
  elseif target == Ephem.NutationsState
    NDIM = 2;
  else
    NDIM = 3;
  end
  
  buf = s.IPT(1,target) + nrl*s.indexInc + s.indexOffset;
  ncf = s.IPT(2,target);
  na  = s.IPT(3,target);
  table = s.scan;
  
  position = zeros(1,NDIM); % row vector
  velocity = zeros(1,NDIM); % row vector
  acceleration = zeros(1,NDIM); % row vector
  if ncf <= 0
    if s.output > 0
      fprintf(s.output,'ERROR: interpolate DE%s. No coefficients available for target %d\n', ...
        s.de_number,target);
    end
    return;
  end
  %   Get correct sub-interval number for this set of coefficients and
  %   then get normalized Chebyshev time within that subinterval.
  dna = double(na);
  temp = dna * double(t);

  %   'tc' is the normalized Chebyshev time (-1 <= tc <= 1).

  tc = 2.0 * mod(temp,1.0) - 1.0;

  %   Check to see whether Chebyshev time has changed, and compute new
  %   polynomial values if it has.  (The element PC[1] is the value of
  %   t1[tc] and hence contains the value of 'tc' on the previous call.)

  if tc ~= PC(2)
    NP = 3;
    TWOT = tc + tc;
    PC(2) = tc; % T_1(x) = x
    PC(3) = TWOT * PC(2) - PC(1); % T_2(x) = 2x^2-1
    VC(3) = 2.0 * TWOT; % T'_2(x) = 4x
  end

  %   Be sure that at least 'ncf' polynomials have been evaluated and
  %   are stored in the array 'PC'.
  %   Chebyshev polynomial recurrence relationship PC[i+1] = T_i(tc)

  if NP < ncf
    for i = (NP+1):ncf
      % T_i-1(x) = 2x T_i-2(x) - T_i-3(x)
      PC(i) = TWOT * PC(i-1) - PC(i-2);
      % T'_i-1(x) = 2x T'_i-2(x) + 2 T_i-2(x) - T'_i-3(x)
      VC(i) = TWOT * VC(i-1) + PC(i-1) + PC(i-1) - VC(i-2);
      % T"_i-1(x) = 2x T"_i-2(x) + 4 T'_i-2(x) - T"_i-3(x)
      AC(i) = TWOT * AC(i-1) + 4*VC(i-1) - AC(i-2);
    end
    NP = ncf;
  end

  vfac = 2.0 * dna;
  %   Interpolate to get position for each component.

  k = buf + int32(floor(temp)) * (NDIM * ncf);
  if k + NDIM*ncf > length(table)
    if s.output > 0
      fprintf(s.output,'interpolate: past end of coefficients by %d/%d\n',...
        k + NDIM*ncf-length(table),length(table));
    end
    return;
  end
  % p = sum( a_i * T_i(x) )
  % v = sum( a_i * T'_i(x) )*2*n
  % a = sum( a_i * T"_i(x) )*4*n^2
  for i = 1:NDIM
    p = 0.0;
    v = 0.0;
    a = 0.0;
    for j = ncf:-1:1
      p = p + PC(j) * table(k + j);
      v = v + VC(j) * table(k + j);
      a = a + AC(j) * table(k + j);
    end
    position(i) = p;
    velocity(i) = v*vfac;
    acceleration(i) = a*vfac*vfac;
    k = k + ncf;
  end

end

