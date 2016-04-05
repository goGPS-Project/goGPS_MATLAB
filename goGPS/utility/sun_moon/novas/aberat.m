function pos2 = aberat (pos1, ve, tlight)

% this function corrects position vector for aberration of light.
% algorithm includes relativistic terms.  adapted from murray (1981)
% mon. notices royal ast. society 195, 639-648.

% input

%  pos1   = position vector of observed object, with reespect to
%           origin at observer (or the geocenter), components in au

%  ve     = velocity vector of observer (or the geocenter),
%           with respect to origin at solar system barycenter,
%           components in au/day (in)

%  tlight = light time from body to observer (or the geocenter) in days

% output

%  pos2 = position vector of observed object, with respect to
%         origin at observer (or the geocenter), corrected
%         for aberration, components in au

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% light-time for one astronomical unit in seconds, from de-405

ausec = 499.0047838061d0;

% speed of light in au/day

c = 86400.0d0 / ausec;

tl = tlight;

p1mag = tl * c;

if (tl == 0.0d0)

    p1mag = sqrt(pos1(1)^2 + pos1(2)^2 + pos1(3)^2);

    tl = p1mag / c;
end

vemag = sqrt(ve(1)^2 + ve(2)^2 + ve(3)^2);

beta = vemag / c;

rdotv = pos1(1) * ve(1) + pos1(2) * ve(2) + pos1(3) * ve(3);

cosd = rdotv / (p1mag * vemag);

gammai = sqrt(1.0d0 - beta^2);

p = beta * cosd;

q = (1.0d0 + p / (1.0d0 + gammai)) * tl;

r = 1.0d0 + p;

for j = 1:1:3
    
    pos2(j) = (gammai * pos1(j) + q * ve(j)) / r;
    
end

