function mode = getmod

% this function returns the 'mode' value, which
% determines the method used for the computation of sidereal
% time and the terrestrial-to-celestial transformation, and the
% accuracy of nutation and related calculations.

%  mode = selection for method and accuracy (out)

%         mode = 0 means cio-based method, full accuracy
%         mode = 1 means cio-based method, reduced accuracy
%         mode = 2 means equinox-based method, full accuracy
%         mode = 3 means equinox-based method, reduced accuracy

global imode

mode = imode;