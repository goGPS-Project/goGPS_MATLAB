function setmod (mode)

% this function allows the user to specify the 'mode' value,
% which determines the method used for the computation of sidereal
% time and the terrestrial-to-celestial transformation, and the
% accuracy of nutation and related calculations.

%  mode = selection for method and accuracy (in)

%         set mode = 0 for cio-based method, full accuracy
%         set mode = 1 for cio-based method, reduced accuracy
%         set mode = 2 for equinox-based method, full accuracy
%         set mode = 3 for equinox-based method, reduced accuracy

% note: other entry points are provided to allow the method and
% accuracy to be specified in a more obvious way:

%  mode = 0 can be set by call ciotio and call hiacc
%  mode = 1 can be set by call ciotio and call loacc
%  mode = 2 can be set by call eqinox and call hiacc
%  mode = 3 can be set by call eqinox and call loacc

global imode lmode

lmode = imode;

imode = mode;





