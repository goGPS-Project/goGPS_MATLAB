function [ems_data_available] = check_ems_availability(GPS_time_fc, GPS_time_ltc, GPS_time_ic, week_R, time_R)

% SYNTAX:
%   [ems_data_available] = check_ems_availability(GPS_time_fc, GPS_time_ltc, GPS_time_ic, week_R, time_R);
%
% INPUT:
%   GPS_time_fc  = EMS prc data GPS time
%   GPS_time_ltc = EMS long term corrections data GPS time
%   GPS_time_ic  = EMS ionospheric corrections data GPS time
%   week_R = reference vector of GPS week numbers
%   time_R = reference vector of GPS time of week
%
% OUTPUT:
%   ems_data_available = boolean flag for data availability check
%
% DESCRIPTION:
%   Function that checks the availability of EMS data 15 minutes before
%   and after the survey timespan.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
% Adapted by Eugenio Realini, 2013
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

%buffer
buf_min = 15; %minutes
buf_dtn = datenum([0 0 0 0 buf_min 0]);

%survey start and end timings (plus buffer)
ts = week_R(1)*604800   + round(time_R(1))   - buf_dtn;     
te = week_R(end)*604800 + round(time_R(end)) + buf_dtn; 

%EMS corrections timings
t_fc  = GPS_time_fc(:,1)*604800  + GPS_time_fc(:,2);
t_ltc = GPS_time_ltc(:,1)*604800 + GPS_time_ltc(:,2);
t_ic  = GPS_time_ic(:,1)*604800  + GPS_time_ic(:,2);

%first and last useful prc data
% start_ems02 = find(t_fc <= ts, 1, 'last');
% end_ems02   = find(t_fc >  te, 1, 'first');

[dist_start_ems02, start_ems02] = min(abs(t_fc-ts));
[dist_end_ems02,     end_ems02] = min(abs(t_fc-te));

%first and last useful long term corrections data
% start_ems25 = find(t_ltc <= ts, 1, 'last');
% end_ems25   = find(t_ltc >  te, 1, 'first');

[dist_start_ems25, start_ems25] = min(abs(t_ltc-ts));
[dist_end_ems25,     end_ems25] = min(abs(t_ltc-te));

%first and last useful ionospheric corrections data
% start_ems26 = find(t_ic <= ts, 1, 'last');
% end_ems26   = find(t_ic >  te, 1, 'first');

[dist_start_ems26, start_ems26] = min(abs(t_ic-ts));
[dist_end_ems26,     end_ems26] = min(abs(t_ic-te));

clear ts te t_fc t_ltc t_ic 

%boolean flag
ems_data_available = 1;

%seconds in 1 day
sday = 86400;

%check if EMS files contain all needed data
if ( isempty(start_ems02) | isempty(end_ems02) | ...
     isempty(start_ems26) | isempty(end_ems26) | ...
     isempty(start_ems25) | isempty(end_ems25) | ...
     dist_start_ems02/sday > 1 | dist_end_ems02/sday > 1 | ...
     dist_start_ems25/sday > 1 | dist_end_ems25/sday > 1 | ...
     dist_start_ems26/sday > 1 | dist_end_ems26/sday > 1 )
 
    ems_data_available = 0;
 
    fprintf('EMS data not sufficient: additional data are needed.\n')
end
