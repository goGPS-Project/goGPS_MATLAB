function pos2 = grvd (pos1, pobs, pbody, rmass)

% this function corrects position vector for the deflection
% of light in the gravitational field of an arbitrary body.  adapted
% from murray (1981) mon. notices royal ast. society 195, 639-648.
% see also formulae in section b of the astronomical almanac, or
% kaplan et al. (1989) astronomical journal 97, 1197-1210, section
% iii f.  this function is valid for an observed body within the
% solar system as well as for a star.

%      pos1   = position vector of observed object, with respect to
%               origin at observer (or the geocenter), components
%               in au (in)
%      pobs   = position vector of observer (or the geocenter),
%               with respect to origin at solar system barycenter,
%               components in au (in)
%      pbody  = position vector of gravitating body, with respect to
%               origin at solar system barycenter, components
%               in au (in)
%      rmass  = reciprocal mass of gravitating body in solar mass
%               units, that is, sun mass / body mass (in)
%      pos2   = position vector of observed object, with respect to
%               origin at observer (or the geocenter), corrected for
%               gravitational deflection, components in au (out)

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% get c, the speed of light in meters/second

c = astcon ('c', 1.0d0);

% get mau, the length of the au in meters

mau = astcon ('au', 1.0d0);

% get gs, the heliocentric gravitational constant

gs = astcon ('gs', 1.0d0);

% construct vector pq from gravitating body to observed object and
% construct vector pe from gravitating body to observer

for j = 1:3
    
    pq(j) = pobs(j) + pos1(j) - pbody(j);
    
    pe(j) = pobs(j) - pbody(j);
    
end

% compute vector magnitudes and unit vectors

pmag = sqrt(pos1(1)^2 + pos1(2)^2 + pos1(3)^2);

emag = sqrt(pe(1)^2 + pe(2)^2 + pe(3)^2);

qmag = sqrt(pq(1)^2 + pq(2)^2 + pq(3)^2);

for j = 1:3
    
    phat(j) = pos1(j) / pmag;
    
    ehat(j) =   pe(j) / emag;
    
    qhat(j) =   pq(j) / qmag;
    
end

% compute dot products of vectors

pdotq = phat(1) * qhat(1) + phat(2) * qhat(2) + phat(3) * qhat(3);

edotp = ehat(1) * phat(1) + ehat(2) * phat(2) + ehat(3) * phat(3);

qdote = qhat(1) * ehat(1) + qhat(2) * ehat(2) + qhat(3) * ehat(3);

% if gravitating body is observed object, or is on a straight line
% toward or away from observed object to within 1 arcsec,
% deflection is set to zero

if (abs (edotp) > 0.99999999999d0)
    
    for j = 1:3
        
        pos2(j) = pos1(j);
        
    end
    
    return
    
end

% compute scalar factors

fac1 = 2.0d0 * gs / (c * c * emag * mau * rmass);

fac2 = 1.0d0 + qdote;

% construct corrected position vector pos2

for j = 1:3
    
    p2j = phat(j) + fac1 * (pdotq * ehat(j) - edotp * qhat(j)) / fac2;
    
    pos2(j) = p2j * pmag;
    
end

