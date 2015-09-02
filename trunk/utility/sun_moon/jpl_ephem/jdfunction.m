function fx = jdfunction (jdin)

% objective function for tdb2utc

% input

%  jdin = current value for UTC julian date

% output

%  fx = delta julian date

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdsaved

tai_utc = findleap(jdin);

fx = utc2tdb (jdin, tai_utc) - jdsaved;

end