function [pos, tlight] = littim (tjd, idbody, pose, tlite)

% this function computes the position of a solar system body,
% as antedated for light-time.

%      tjd    = tdb julian date of observation (in)

%      idbody = id number of body, used in calls to solsys (in)

%      pose   = position vector of observer (or the geocenter),
%               with respect to origin at solar system barycenter,
%               referred to icrs axes, components in au (in)

%      tlite  = first approximation to light-time, in days (in)
%               (can be set to 0.0d0 if unknown)

%      pos    = position vector of body, with respect to origin at
%               observer (or the geocenter), referred to icrs axes,
%               components in au (out)

%      tlight = final light-time, in days (out)

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% set light-time convergence tolerance

tol = 1.0d-9;

t0 = 0.0d0;

t1 = tjd - t0;

t2 = t1 - tlite;

% iterate to obtain correct light-time (usually converges rapidly)

for iter = 1:10
    
    [pos1, vel1, ierr] = solsys (t2, str2num(idbody), 0);
    
    [pos, tlight] = geocen (pos1, pose);
    
    if (ierr ~= 0)
        
        fprintf ('\nplace: cannot obtain coordinates of object at jd %16.8f', t0 + t2);
        
        return
        
    end
    
    t3 = t1 - tlight;
    
    if (abs(t3 - t2) > tol)
        
        t2 = t3;
        
    else
        
        break
        
    end
    
end


