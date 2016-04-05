function pos2 = grvdef (tjd, loc, pos1, pobs)

% this function computes the total gravitational deflection of
% light for the observed object due to the major gravitating bodies
% in the solar system. this function is valid for an observed body
% within the solar system as well as for a star. see klioner
% (2003), astronomical journal 125, 1580-1597, section 6.

%      tjd    = tdb julian date of observation

%      loc    = code for location of observer, determining
%               whether the gravitational deflection due to the
%               earth itself is applied (in)
%               set loc=0 for no earth deflection (normally means
%                         observer is at geocenter)
%               set loc=1 to add in earth deflection (normally
%                         means observer is on or above surface
%                         of earth, including earth orbit)

%      pos1   = position vector of observed object, with respect to
%               origin at observer (or the geocenter), referred
%               to icrs axes, components in au (in)

%      pobs   = position vector of observer (or the geocenter),
%               with respect to origin at solar system barycenter,
%               referred to icrs axes, components in au (in)

%      pos2   = position vector of observed object, with respect to
%               origin at observer (or the geocenter), referred
%               to icrs axes, corrected for gravitational
%               deflection, components in au (out)

% the following list of names identifies which gravitating bodies
% (aside from the earth) are potentially used -- list is taken from
% klioner's table 1, the order based on area of sky affected (col 2)
     
% names = ['sun'; 'jup'; 'sat'; 'moo'; 'ven'; 'ura'; 'nep'];

% change value of nbody to include or exclude gravitating bodies
% (nbody = 0 means no deflection calculated, nbody = 1 means sun only,
% nbody = 2 means sun + jupiter, etc.)

nbody = 3;

% get c, the speed of light in au/day

c = astcon ('c(au/day)', 1.0d0);

% id numbers and reciprocal masses of gravitating bodies

id(1) = idss_novas('sun');

id(2) = idss_novas('jup');

id(3) = idss_novas('sat');

rmass(1) = astcon('mass_sun', 1.0d0);

rmass(2) = astcon('mass_jup', 1.0d0);

rmass(3) = astcon('mass_sat', 1.0d0);

ide = idss_novas ('earth');

rmasse = astcon ('mass_earth', 1.0d0);

% initialize output vector of observed object to equal input vector

for j = 1:3
    
    pos2(j) = pos1(j);
    
end

% option for no deflection

if (nbody <= 0)
    
    return
    
end

% compute light-time to observed object

tlt = sqrt (pos1(1)^2 + pos1(2)^2 + pos1(3)^2) / c;

% cycle through gravitating bodies

for i = 1:nbody
    
    if (id(i) == -9999)
        
        break
        
    end
    
    % get position of gravitating body wrt ss barycenter at time tjd
    
    [pbody, vbody, ierr] = solsys (tjd, id(i), 0);
    
    % get position of gravitating body wrt observer at time tjd
    
    [pbodyo, x] = geocen (pbody, pobs);
    
    % compute light-time from point on incoming light ray that
    % is closest to gravitating body
    
    dlt = dlight (pos2, pbodyo);
    
    % get position of gravitating body wrt ss barycenter at time
    % when incoming photons were closest to it
    
    tclose = tjd;
    
    if (dlt > 0.0d0)
        
        tclose = tjd - dlt;
        
    end
    
    if (tlt < dlt)
        
        tclose = tjd - tlt;
        
    end
    
    [pbody, vbody, ierr] = solsys (tclose, id(i), 0);
    
    % compute deflection due to gravitating body
    
    pos2 = grvd (pos2, pobs, pbody, rmass(i));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if observer is not at geocenter, add in deflection due to earth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (loc ~= 0)
    
    % get position of earth wrt ss barycenter at time tjd
    
    [pbody, vbody, ierr] = solsys (tjd, ide, 0);
    
    % compute deflection due to earth
    
    pos2 = grvd (pos2, pobs, pbody, rmasse);
    
end

