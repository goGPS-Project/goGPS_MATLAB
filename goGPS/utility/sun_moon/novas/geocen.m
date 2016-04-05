function [pos2, tlight] = geocen (pos1, pe)

% this function moves the origin of coordinates from the
% barycenter of the solar system to the observer (or the
% geocenter).  i.e., this function accounts for parallax
% (annual + geocentric or just annual).

% input

%  pos1 = position vector of star or planet, with respect to
%         origin at solar system barycenter, components in au

%  pe   = position vector of observer (or the geocenter),
%         with respect to origin at solar system barycenter,
%         components in au

% output

%   pos2   = position vector of star or planet, with respect to
%            origin at observer (or the geocenter), components in au

%   tlight = light-time from star or planet to observer (or the
%            geocenter) in days

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% speed of light in au/day

c = astcon('c(au/day)', 1.0);

for j = 1:1:3

    pos2(j) = pos1(j) - pe(j);

end

tlight = sqrt(pos2(1)^2 + pos2(2)^2 + pos2(3)^2) / c;





