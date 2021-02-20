function [time, xyz, enu, std_enu, std_3d] = steCrdImporter (file_name)
% SYNTAX:
%   [time, xyz, enu, std_enu, std_3d] = steCRDimporter (file_name);
%
% EXAMPLE:
%   [time, xyz, enu, std_enu, std_3d] = steCRDimporter (file_name);
%
% INPUT:
%   file_name           [n x 1] observation time
%   y           [n x 1] observation value
%               [n x 2] observation value and variances (in this case the variances of the data are taken into account)
%   dxs         [1 x 1] spline base size
%   regFactor   [1 x 1] regularization factor on the first derivative
%   x_ext       [m x 1] points in which to compute the interpolation
%
% OUTPUT:
%   ySplined    [n x 1] interpolated observation (on the observation epochs)
%   xS          [o x 1] center of the (o) splines
%   wS          [o x 1] weights of the splines
%   y_ext       [m x 1] spline interpolated in x_ext positions
%
% DESCRIPTION:
%   Interpolate with cubic splines a given dataset
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
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

    log = Core.getLogger();

    fid = fopen(file_name,'r');

    if (fid < 0)
        log.addError(['Failed to open ', file_name]);
        time = [];
        xyz = [];
        enu = [];
        std_enu = [];
        std_3d = [];
    else
        txt = fread(fid,'*char')';
        log.addMessage(['Reading ', file_name]);
        fclose(fid);

        data = textscan(txt(392:end)',' %4d-%1d %4d-%03d %4d-%2d-%2d %1s %2d:%2d:%2d %2d:%2d:%2d %4s %6s %15f %15f %15f %15f %15f %15f %1s %8f %8f %8f %8f ');
        time = GPS_Time(datenum(double([data{:,[5:7]} data{:,[9:11]}])));
        xyz = [data{:,17:19}];
        enu = [data{:,20:22}];
        std_enu = [data{:,24:26}];
        std_3d = data{:,27};
    end
end
