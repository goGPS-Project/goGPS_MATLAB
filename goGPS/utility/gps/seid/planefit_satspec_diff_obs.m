function [til_obs] = planefit_satspec_diff_obs(diff_obs, commontime, ipp_lon, ipp_lat, PRN, target_sta, diff_phase_flag)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------


%find idx to execute interpolation
not_nan_idx=~isnan(ipp_lon(target_sta,:));

timestamps = 1 : length(commontime);
timestamps = timestamps(not_nan_idx);
interp_dobs = diff_obs(PRN, :, target_sta);
n_sta = length(ipp_lon(:,1));

for i=1:length(timestamps)
    if timestamps(i)>length(interp_dobs);break;end
    %verify the stations used for interpolation
    %notnan stations
    plane_idx=~isnan(ipp_lon(:,timestamps(i)));

    %exclude L1 stations
    plane_idx(target_sta)=0;

    if any(plane_idx)==1
        %idx for all stations
        raw_idx=1:n_sta;

        %stations idx used for interpolation
        interp_idx=intersect(raw_idx,raw_idx(plane_idx));

        %execute interpolation
        [interp_dobs(timestamps(i)),~]=plane_fitting(ipp_lon(interp_idx,timestamps(i)),ipp_lat(interp_idx,timestamps(i)),squeeze(diff_obs(PRN,timestamps(i),interp_idx)),ipp_lon(target_sta,timestamps(i)),ipp_lat(target_sta,timestamps(i)));

    else
        %apply value of 1 epoch ago
        if i>2 && ~isnan(interp_dobs(timestamps(i-1)))==1
            interp_dobs(timestamps(i))=interp_dobs(timestamps(i-1));
        else
            interp_dobs(timestamps(i))=NaN;
        end
    end
end

til_obs=NaN(1,length(commontime));
if (diff_phase_flag)
    %compute sum of ~dobs
    notnan_id=find(not_nan_idx==1);
    %check if isnan_idx exceeds size of interp_dobs
    if any(notnan_id>length(interp_dobs))
        not_nan_idx(notnan_id(notnan_id>length(interp_dobs)))=0;
    end

    dobs_for_sum=interp_dobs(not_nan_idx);
    %  search NaN and replace to zero
    nan_find=isnan(dobs_for_sum);
    if (any(nan_find))==1
        dobs_for_sum(nan_find==1)=0;
    end

    sum_dobs=zeros(1,length(dobs_for_sum));
    for i=2:length(dobs_for_sum)
        sum_dobs(i)=sum(dobs_for_sum(1:i-1));
    end

    til_obs(not_nan_idx) = sum_dobs;
else
    til_obs(not_nan_idx) = interp_dobs(timestamps);
end
