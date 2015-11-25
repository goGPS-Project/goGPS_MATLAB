function tmatrix = j2000_icrs(itype)

% transformation matrix from the mean dynamical
% equator and equinox at j2000 to the
% International Celestial Reference System (ICRS)

% input

%  itype = type of transformation
%      1 = LLR + VLBI
%      2 = Chapront et al. LLR

% output

%  tmatrix = transformation matrix

% reference: Rotation Matrix from the Mean Dynamical
% Equator and Equinox at J2000.0 to the ICRS, Astronomy
% and Astrophysics, October 1, 2003.

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (itype == 1)
    % LLR and VLBI
    
    tmatrix(1, 1) = 0.99999999999999332;
    tmatrix(1, 2) = 7.2e-8;
    tmatrix(1, 3) = -8.056e-8;

    tmatrix(2, 1) = -7.2e-8;
    tmatrix(2, 2) = 0.99999999999999332;
    tmatrix(2, 3) = -3.306e-8;

    tmatrix(3, 1) = 8.056e-8;
    tmatrix(3, 2) = 3.306e-8;
    tmatrix(3, 3) = 0.999999999999996208;
end

if (itype == 2)
    % Chapront et al. LLR
    
    tmatrix(1, 1) = 0.9999999999999938;
    tmatrix(1, 2) = 7.1e-8;
    tmatrix(1, 3) = -8.6e-8;

    tmatrix(2, 1) = -7.1e-8;
    tmatrix(2, 2) = 0.9999999999999938;
    tmatrix(2, 3) = -2.6e-8;

    tmatrix(3, 1) = 8.6e-8;
    tmatrix(3, 2) = 2.6e-8;
    tmatrix(3, 3) = 0.99999999999999598;
end
