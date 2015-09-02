function svprint(r, v)

% print position and velocity vectors and magnitudes

% input

%  r = position vector (kilometers)
%  v = velocity vector (kilometers/second)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position magnitude (kilometers)

rmag = norm(r);

% velocity magnitude (kilometers)

vmag = norm(v);

% print state vector and magnitudes

fprintf ('\n        rx (km)                 ry (km)                rz (km)                rmag (km)');

fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e  \n', r(1), r(2), r(3), rmag);

fprintf ('\n        vx (kps)                vy (kps)               vz (kps)               vmag (kps)');

fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e  \n\n', v(1), v(2), v(3), vmag);

