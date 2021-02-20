function diflt = dlight (pos1, pe)

% this function returns the difference in light-time, for a star,
% between the barycenter of the solar system and the observer (or
% the geocenter).

% input

%  pos1  = position vector of star, with respect to origin at
%          solar system barycenter

%  pe    = position vector of observer (or the geocenter),
%          with respect to origin at solar system barycenter,
%          components in au

% output

%  diflt = difference in light time, in the sense star to
%          barycenter minus star to earth, in days

% -or-

% this function returns the light-time from the observer (or the
% geocenter) to a point on a light ray that is closest to a
% specific solar system body.

% input

%  pos1 = position vector toward observed object, with respect
%         to origin at observer (or the geocenter)

%  pe   = position vector of solar system body, with respect
%         to origin at observer (or the geocenter), components in au

% output

%  diflt = light time to point on line defined by pos1 that is
%          closest to solar system body (positive if light
%          passes body before hitting observer, i.e., if
%          pos1 is within 90 degrees of pe)

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% speed of light in au/day

c = astcon('c(au/day)', 1.0d0);

% from pos1, form unit vector u1 in direction of star or light source

dis = sqrt(pos1(1)^2 + pos1(2)^2 + pos1(3)^2);

for j = 1:1:3

    u1(j) = pos1(j) / dis;

end

% light-time returned is the projection of vector pe onto the unit
% vector u1 (formed from pos1), divided by the speed of light

diflt = (pe(1) * u1(1) + pe(2) * u1(2) + pe(3) * u1(3)) / c;


