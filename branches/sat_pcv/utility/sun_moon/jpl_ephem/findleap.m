function leapsecond = findleap(jdate)

% find number of leap seconds for utc julian date

% input

%  jdate = utc julian date

% input via global

%  jdateleap = array of utc julian dates
%  leapsec   = array of leap seconds

% output

%  leapsecond = number of leap seconds

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdateleap leapsec

ndata = length(jdateleap);

if (jdate <= jdateleap(1))
    % date is <= 1972
    
    leapsecond = leapsec(1);
elseif (jdate >= jdateleap(ndata))
    % date is >= end of current data
    
    leapsecond = leapsec(ndata);
else
    % find data within table
    
   for i = 1:1:ndata - 1
       if (jdate >= jdateleap(i) && jdate < jdateleap(i + 1))
           leapsecond = leapsec(i);
           
           break;
       end
   end
end

    