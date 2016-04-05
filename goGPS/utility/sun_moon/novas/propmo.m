function pos2 = propmo (tjd1, pos1, vel1, tjd2)

% this function applies proper motion, including foreshortening
% effects, to a star's position.

% input

%  tjd1 = tdb julian date of first epoch

%  pos1 = position vector of star at first epoch

%  vel1 = velocity vector of star at first epoch

%  tjd2 = tdb julian date of second epoch

% output

%  pos2 = position vector of star at second epoch

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:1:3
    
    pos2(j) = pos1(j) + vel1(j) * (tjd2 - tjd1);

end
