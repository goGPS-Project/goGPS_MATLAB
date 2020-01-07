function [pos_R_crd, flag_XR, pos_M_crd, flag_XM] = load_CRD(filename_sta, marker_R, marker_M)

% SYNTAX:
%   [pos_R_crd, flag_XR, pos_M_crd, flag_XM] = load_CRD(filename_sta, marker_R, marker_M);
%
% INPUT:
%   filename_sta = path to the stations coordinates file [string]
%   marker_R = rover marker name [string]
%   marker_M = master marker name [string]
%
% OUTPUT:
%   pos_R_crd = rover ECEF coordinates [3x1 vector]
%   flag_XR   = 0: unknown
%               1: approximated
%               2: fixed
%   pos_M_crd = master ECEF coordinates [3x1 vector]
%   flag_XM   = 0: unknown
%               1: approximated
%               2: fixed
%
% DESCRIPTION:
%   Tool for loading .CRD files: stations coordinates.

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

pos_R_crd = [];
flag_XR = 0;
pos_M_crd = [];
flag_XM = 0;

if (isempty(filename_sta))
    return
end

%parse file containing station coordinates
[markers, coords_X, coords_Y, coords_Z, flags] = textread(filename_sta,'%s%f%f%f%d'); %#ok<DTXTRD>

% master
%find the correct marker
marker_idx = find(strcmp(markers, marker_M),1,'last');

if ~isempty(marker_idx)
    %extract the corresponding coordinates
    XM = coords_X(marker_idx);
    YM = coords_Y(marker_idx);
    ZM = coords_Z(marker_idx);

    %set master station position
    pos_M_crd = [XM; YM; ZM];

    flag_XM = flags(marker_idx);
end

% rover
%find the correct marker
pos_R_crd = zeros(3,1,size(marker_R,3));
flag_XR = zeros(size(marker_R,3),1);
for m = 1 : size(marker_R,3)
    if (iscell(marker_R))
        target_marker = marker_R{1,1,m};
    else
        target_marker = marker_R;
    end
    marker_idx = find(strcmpi(markers, target_marker),1,'last');

    if ~isempty(marker_idx)
        %extract the corresponding coordinates
        XR = coords_X(marker_idx);
        YR = coords_Y(marker_idx);
        ZR = coords_Z(marker_idx);

        %set rover station position
        pos_R_crd(:,1,m) = [XR; YR; ZR];

        flag_XR(m) = flags(marker_idx);
    end
end
flag_XR = min(flag_XR);
