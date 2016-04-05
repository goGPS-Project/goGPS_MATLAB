function [alimb, afrac] = limang (pos1, poso)

% this function determines the angle of an object above or below
% the earth's limb (horizon).  the geometric limb is computed,
% assuming the earth to be an airless sphere (no refraction or
% oblateness is included).  the observer can be on or above the
% earth.  for an observer on the surface of the earth, this
% function returns the approximate unrefracted altitude.

% input

%  pos1 = position vector of observed object, with respect to
%         origin at geocenter, components in au

%  poso = position vector of observer, with respect to origin
%         at geocenter, components in au

% output

%  alimb = angle of observed object above (+) or below (-) limb in degrees

%  afrac = nadir angle of observed object as a fraction of
%          apparent radius of limb

%          afrac < 1.0d0 means below the limb
%          afrac = 1.0d0 means on the limb
%          afrac > 1.0d0 means above the limb

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

halfpi = 0.5d0 * pi;

degcon = 180.0d0 / pi;

au = astcon('au', 1.0d0);

erad = astcon('erad', 1.0d0);

rade = erad / au;

disobj = sqrt(pos1(1)^2 + pos1(2)^2 + pos1(3)^2);

disobs = sqrt(poso(1)^2 + poso(2)^2 + poso(3)^2);

% compute apparent angular radius of earth's limb

if (disobs > rade)

    aprad = asin(rade / disobs);

else

    aprad = halfpi;

end

% compute zenith distance of earth's limb

zdlim = pi - aprad;

% compute zenith distance of observed object

coszd = (pos1(1) * poso(1) + pos1(2) * poso(2) ...
    + pos1(3) * poso(3)) / (disobj * disobs);

if (coszd <= -1.0d0)

    zdobj = pi;

elseif (coszd >= 1.0d0)

    zdobj = 0.0d0;

else

    zdobj = acos(coszd);

end

% angle of object wrt limb is difference in zenith distances

alimb = (zdlim - zdobj) * degcon;

% nadir angle of object as a fraction of angular radius of limb

afrac = (pi - zdobj) / aprad;

