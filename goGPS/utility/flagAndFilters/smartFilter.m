function data_out = smartFilter(data_in, min_arc)
% SYNTAX:
%   data_out = smartFilter(data_in, min_arc)
%
% INPUT:
%   data_in = dataset
%   min_arc = minimum arc length to be kept
%
% OUTPUT:
%   data_out = filtered dataset
%
% DESCRIPTION:
%   The software has been calibrated to filter pseudo-range data

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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

in_size = size(data_in);

% serialize as column array the input data
data_in = data_in(:);

% set to true to ignore small jumps in the detection of the intervals
reduce_jumps = false;
max_gap_size = 3;        % If I loss the satellite for max_gap_size try to fill the gap
med_win_size = 7;        % The signal to be spline filtered is the one previously filtered by a median filter
spline_size = 10;        % Size of the spline base (a cubic spline have a total width of 4 x spline_size epochs)
jump_detection_thr = 20; % Minimum jump to be recognized
jump_detection_win = 31; % Size of the moving windows to filter the signal and recognize jump
%--------------------------------------------------------------------------
% First step: find jumps (they are in proximity of NaN)

% find  intervals of NaNs
fI = getOutliers(isnan(data_in));

% trim the NaNs at the border of the dataset
% I cannot process them
rNaN  = [0; 0];
if (~isempty(fI) && fI(end) == numel(data_in))
    rNaN(2) = length(data_in) -  fI(end, 1) + 1;
    data_in = data_in(1:fI(end, 1)-1);
end
if (~isempty(fI) && fI(1) == 1)
    data_in = data_in((fI(1, 2)+1):end);
    rNaN(1) = fI(1, 2);
end
% To restore the full size of the set e.g. data_out = [nan(rNaN(1),1); data_out(:); nan(rNaN(2),1)];

ids = 1:length(data_in);
inan = isnan(data_in);

% find intervals of NaNs
fI = getOutliers(isnan(data_in));

% find the intervals of good data_in
vI = [[1; fI(:,2)+1] [fI(:,1)-1; numel(data_in)]];

% for each valid interval compute the median interpolant line
data_m = zeros(size(data_in));
s = 0;
for i = 1 : size(vI,1)
    tmp = data_in(vI(i,1):vI(i,2));
    % DEBUG: figure(i); clf;
    % DEBUG: plot(tmp - median(tmp), 'o-'); hold on;
    % DEBUG: plot(medfilt_mat(tmp - median(tmp), 3)); setLastLineWidth(2);
    % DEBUG: plot(medfilt_mat(medfilt_mat(tmp - median(tmp), 3), 13)); setLastLineWidth(2);
    jump_ids = find(abs(diff(medfilt_mat(tmp,min(jump_detection_win, length(tmp)-(mod(length(tmp)+1,2)))))) > jump_detection_thr);
    sub_segment = [[1; jump_ids + 1] [jump_ids; length(tmp)]];
    % DEBUG: plot(sub_segment(:,1), tmp(sub_segment(:,1)) - median(tmp),'o', sub_segment(:,2),tmp(sub_segment(:,2)) - median(tmp),'*');
    for s = 1 : size(sub_segment)
        if ((sub_segment(s,2) - sub_segment(s,1) + 1) >= min_arc)
            sub_tmp = tmp(sub_segment(s,1):sub_segment(s,2));
            data_m((vI(i,1) + sub_segment(s,1) - 1):(vI(i,1) + sub_segment(s,2) - 1)) = median(sub_tmp(not(isnan(sub_tmp))));
        else
            inan((vI(i,1) + sub_segment(s,1) - 1):(vI(i,1) + sub_segment(s,2) - 1)) = true;
        end
    end

    % DEBUG: [i std(data_in(vI(i,1):vI(i,2)) - median(data_in(vI(i,1):vI(i,2))))]
end
% DEBUG: dockAllFigures;

%--------------------------------------------------------------------------
% % reduce the number of jumps -> if a jump is not significant -> ignore it
%
% remove all the intervals under 0.1
if reduce_jumps
    thr = 20;
    data_m(inan) = 0;
    data_m = medfilt_mat(data_m, 21);
    jump_id = find(abs(diff(data_m)) > thr) + 1;
    segment = [[1; jump_id + 1] [jump_id; length(data_m)]];

    % for each valid interval compute the median interpolant line
    for i = 1 : length(segment)
        tmp = data_in(segment(i,1):segment(i,2));
        data_m(segment(i,1):segment(i,2)) = median(tmp(not(isnan(tmp))));
    end
    clear tmp;
end
data_m(inan) = nan;
%--------------------------------------------------------------------------
% Reduce the dataset using the median interpolant

% remove the reference
% DEBUG: figure(100); clf;
% DEBUG: plot(data_in - data_m,'o');
% DEBUG: hold on;

% Let's use a simple filler based on the median to fill small nans
tmp = data_in - data_m;
tmp(inan) = 0;
tmp = medfilt_mat(tmp, 2*max_gap_size+1);
inan_new = (tmp == 0) & (inan);
y1 = data_in - data_m;
y1(inan) = tmp(inan); % fill with median
clear tmp
y1 = medfilt_mat(y1, med_win_size); % improve median filter
y1(inan_new) = nan;
% DEBUG: plot(y1);

% % I cannot use this filler (it's in C++) :-(
% y1 = fill1D(data_in - data_m, fI, 2);
% plot(y1,'.-');

% smooth the dataset
% and restore the jumps

y1s = nan(size(data_in));
y1s(not(inan)) = splinerMat(ids(not(inan)), y1(not(inan)), spline_size, 0);
% DEBUG: plot(y1s); setLastLineWidth(2);
% DEBUG: ylim([-100 100]);

%--------------------------------------------------------------------------
% Restore the dataset adding back the median interpolant

% DEBUG: figure(200); clf;
% DEBUG: plot(data_in,'o');
% DEBUG: hold on;
% DEBUG: plot(data_m);
data_out = y1s + data_m;
% DEBUG: plot(data_out);
% DEBUG: setLastLineWidth(2);

data_out = [nan(rNaN(1),1); data_out(:); nan(rNaN(2),1)];
data_out = reshape(data_out, in_size(1), in_size(2));
