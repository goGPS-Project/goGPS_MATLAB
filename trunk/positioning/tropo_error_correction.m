function [corr] = tropo_error_correction(el, h)

% SYNTAX:
%   [corr] = tropo_error_correction(el, h);
%
% INPUT:
%   el = satellite elevation
%   h  = receiver ellipsoidal height
%
% OUTPUT:
%   corr = tropospheric error correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to tropospheric refraction.
%   Saastamoinen algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Laboratorio di Geomatica, Polo Regionale di Como,
%    Politecnico di Milano, Italy
%----------------------------------------------------------------------------------------------

global tropo_model

%Saastamoinen model requires positive ellipsoidal height
h(h < 0) = 0;

% If in ATM_model section in config file tropo is set to 0 it does not use tropospheric model
% [ATM_model]
% tropo=0
if (tropo_model == 0)
    corr = zeros(size(el));
else
    if (h < 5000)

	%conversion to radians
	el = abs(el) * pi/180;

	    %Standard atmosphere - Berg, 1948 (Bernese)
	%pressure [mbar]
	Pr = 1013.25;
	%temperature [K]
	Tr = 291.15;
	%numerical constants for the algorithm [-] [m] [mbar]
	Hr = 50.0;

	P = Pr * (1-0.0000226*h).^5.225;
	T = Tr - 0.0065*h;
	H = Hr * exp(-0.0006396*h);

	%----------------------------------------------------------------------

	%linear interpolation
	h_a = [0; 500; 1000; 1500; 2000; 2500; 3000; 4000; 5000];
	B_a = [1.156; 1.079; 1.006; 0.938; 0.874; 0.813; 0.757; 0.654; 0.563];

	t = zeros(length(T),1);
	B = zeros(length(T),1);

	for i = 1 : length(T)

	    d = h_a - h(i);
	    [dmin, j] = min(abs(d));
	    if (d(j) > 0)
		index = [j-1; j];
	    else
		index = [j; j+1];
	    end

	    t(i) = (h(i) - h_a(index(1))) ./ (h_a(index(2)) - h_a(index(1)));
	    B(i) = (1-t(i))*B_a(index(1)) + t(i)*B_a(index(2));
	end

	%----------------------------------------------------------------------

	e = 0.01 * H .* exp(-37.2465 + 0.213166*T - 0.000256908*T.^2);

	%tropospheric error
	corr = ((0.002277 ./ sin(el)) .* (P - (B ./ (tan(el)).^2)) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e);
    else
	corr = zeros(size(el));
    end
end