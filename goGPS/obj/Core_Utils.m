
%   CLASS Core_Utils
% =========================================================================
%
% DESCRIPTION
%   Class to manages utilities
%
% EXAMPLE
%   % set of static utilities
%   Core_Utils.diffAndPred
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, Giulio Tagliaferro, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Core_Utils < handle
    properties (Constant)
        V_LIGHT = 299792458;                % Velocity of light in the void [m/s]
    end
    
    methods (Static)
        %--------------------------------------------------------------------------
        % FILTERS AND INTERPOLATORS
        %--------------------------------------------------------------------------
        
        function diff_data = diffAndPred(data, n_order, t_ref, method)
            % compute diff predicting epoch 0 of each arc
            % using interp 1 spline method
            %
            % SYNTAX
            %   Core_Utils.diffAndPred(data, n_order, t_ref, method)
            
            if nargin < 3 || isempty(t_ref)
                t_ref = 1 : size(data,1);
            end
            if nargin < 2 || isempty(n_order)
                n_order = 1;
            end
            if nargin < 4 || isempty(method)
                method = 'spline';
            end
            diff_data = nan(size(data));
            % Add n_order rows to data
            data = [repmat(data(1,:), n_order, 1); data];
            for s = 1 : size(data, 2)
                % Get the original good data for column s
                tmp = data(1 + n_order : end, s);
                id_ok = ~isnan(tmp);
                if sum(id_ok) > 2
                    lim = getOutliers(id_ok);
                    % Interpolate data beginning
                    % interpolate the "left" of the first element of an arc
                    % because diff "eat" the first value
                    %if (length(id_ok) > (n_order + 1)) && any(id_ok(1))
                    %    id_est = find(id_ok(lim(1,1):lim(1,2)));
                    %    data(1 : n_order, s) = interp1(t_ref(id_est), tmp(id_est), 1 - n_order : 0, 'spline', 'extrap');
                    %end
                    
                    lim_short = lim(lim(:,2) - lim(:,1) < 2 & lim(:,1) > 1, :);
                    % short arcs cannot be differenciated efficiently
                    for l = 1 : size(lim_short, 1)
                        data(lim_short(l, 1), s) = data(lim_short(l, 1)+1, s);
                    end
                    
                    % differenciate only limits larger than 2
                    lim = lim(lim(:,2) - lim(:,1) > 1, :);
                    for l = 1 : size(lim, 1)
                        id_data = lim(l, 1) : lim(l, 2);
                        id_est = 0 : (n_order - 1);
                        
                        % slower approach with interp1
                        % data(lim(l, 1) + id_est, s) = interp1(t_ref(id_data), tmp(id_data), lim(l, 1) - 1 - fliplr(id_est), 'spline', 'extrap');
                        
                        % faster approach skipping a lot of checks
                        % this is the internal implementation of interp1
                        if strcmp(method, 'zeros')
                            data(lim(l, 1) + id_est, s) = 0;
                        else
                            fun = griddedInterpolant(t_ref(id_data), tmp(id_data), method);
                            data(lim(l, 1) + id_est, s) = fun(lim(l, 1) - 1 - fliplr(id_est));
                        end
                        
                        diff_data(id_data, s) = diff(data(lim(l, 1) : (lim(l, 2) + n_order), s), n_order);
                        % restore data for the next interval
                        data(1 + n_order : end, s) = tmp;
                    end
                end
            end
            % diff_data = diff(data, n_order); % now it is done arc by arc
        end
        
        
        function y_out = interp1LS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 1D dataset
            %
            % SYNTAX
            %   y_out = interp1LS(x_in, y_in, degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                if size(y_in, 1) == 1
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~isnan(x_in(:));
                x_tmp = x_in(id_ok);
                y_tmp = y_tmp(id_ok);
                
                n_obs = numel(x_tmp);
                A = zeros(n_obs, degree + 1);
                A(:, 1) = ones(n_obs, 1);
                for d = 1 : degree
                    A(:, d + 1) = x_tmp .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = numel(x_out);
                    A2 = zeros(n_out, degree + 1);
                    A2(:, 1) = ones(n_out, 1);
                    for d = 1 : degree
                        A2(:, d + 1) = x_out .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    %y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = A2 * (A \ y_tmp(:));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 2));
                else
                    y_out(:,c) = A2 * ((A' * A + 1e-6 * eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')
            end
        end
        
        function [y_out] = interpPlaneLS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 3D dataset
            % Estimate a plane on x, y + polynomial on z
            % y_in is an array containing m set of data (1 set per column)
            % e.g. degree = 2
            %   f(x,y,z) =  a*x + b*y + c + d*z^2 + e*z
            %
            % SYNTAX
            %   y_out = interpPlane1LS(x_in, y_in, z_degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                % try to correct orientation of x
                if size(y_in, 1) == 1 && size(y_tmp, 1) ~= size(x_tmp, 1)
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~any(isnan(x_in), 2);
                x_tmp = x_in(id_ok, :);
                y_tmp = y_tmp(id_ok);
                
                n_obs = size(x_tmp, 1);
                A = zeros(n_obs, degree + 3);
                A(:, 1) = ones(n_obs, 1);
                A(:, 2) = x_tmp(:, 1); % y
                A(:, 3) = x_tmp(:, 2); % x
                for d = 1 : degree
                    A(:, d + 3) = x_tmp(:, 3) .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = size(x_out, 1);
                    A2 = zeros(n_out, degree + 3);
                    A2(:, 1) = ones(n_out, 1);
                    A2(:, 2) = x_out(:, 1); % y
                    A2(:, 3) = x_out(:, 2); % x
                    for d = 1 : degree
                        A2(:, d + 3) = x_out(:, 3) .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 3));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')
            end
        end
        
        function [y_out] = interp22nLS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 3D dataset
            % Estimate a plane on x, y + polynomial on z
            % y_in is an array containing m set of data (1 set per column)
            % e.g. degree = 2
            %   f(x,y,z) =  a*x + b*y + c + d*z^2 + e*z
            %
            % SYNTAX
            %   y_out = interpPlane1LS(x_in, y_in, z_degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                % try to correct orientation of x
                if size(y_in, 1) == 1 && size(y_tmp, 1) ~= size(x_tmp, 1)
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~any(isnan(x_in), 2);
                x_tmp = x_in(id_ok, :);
                y_tmp = y_tmp(id_ok);
                
                n_obs = size(x_tmp, 1);
                A = zeros(n_obs, degree + 6);
                A(:, 1) = ones(n_obs, 1);
                A(:, 2) = x_tmp(:, 1); % y
                A(:, 3) = x_tmp(:, 2); % x
                A(:, 4) = x_tmp(:, 1) .^ 2; % x^2
                A(:, 5) = x_tmp(:, 2) .^ 2; % y^2
                A(:, 6) = x_tmp(:, 1) .* x_tmp(:, 2); % x*y
                for d = 1 : degree
                    A(:, d + 6) = x_tmp(:, 3) .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp) &&  (sum(x_out(:) - x_tmp(:)) == 0)
                    A2 = A;
                else
                    n_out = size(x_out, 1);
                    A2 = zeros(n_out, degree + 6);
                    A2(:, 1) = ones(n_out, 1);
                    A2(:, 2) = x_out(:, 1); % y
                    A2(:, 3) = x_out(:, 2); % x
                    A2(:, 4) = x_out(:, 1) .^ 2; % x^2
                    A2(:, 5) = x_out(:, 2) .^ 2; % y^2
                    A2(:, 6) = x_out(:, 1) .* x_out(:, 2); % x*y
                    for d = 1 : degree
                        A2(:, d + 6) = x_out(:, 3) .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 3));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')
            end
        end
        
        function val = linInterpLatLonTime(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % Interpolate values froma data on a gepgraphical grid with multiple epoch
            % data structure:
            %        first dimension : dlat (+) south pole -> north pole
            %        second dimension : dlon (+) west -> east
            %        third dimension : dr (+) time usual direction
            %        NOTE: dlat, dlon,dt do not have to be positive
            %
            % INPUT:
            %      data - the data to be interpolate
            %      fist_lat - value of first lat value (max lat)
            %      dlat - px size lat
            %      first_lon - value of first lon value
            %      dlon - px size lon
            %      first_t - value of first time
            %      dt - px size time
            %      lat - lat at what we want to interpolate
            %      lon - lon at what we ant to interpolate
            %      gps_time - time at what we want to interpolate
            % NOTES 1 - all lat values should have same unit of measure
            %       2 - all lon values should have same unit of measure
            %       3 - all time values should have same unit of measure
            %       4 - the method will interpolate first in the dimesnion with less time
            % IMPORTANT : no double values at the borders should coexist: e.g. -180 180 or 0 360
            [nlat , nlon, nt] = size(data);
            n_in_lat = length(lat);
            n_in_lon = length(lon);
            n_in_t = length(t);
            assert(n_in_lat == n_in_lon);
            [ it, st, ilons, ilone, slon, ilat, slat] = Core_Utils.getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon,t);
            if n_in_lat > n_in_t % time first
                
                it = it*ones(size(ilat));
                % interpolate along time
                % [ 1 2  <= index of the cell at the smae time
                %   3 4]
                idx1 = sub2ind([nlat nlon nt], ilat, ilons, it);
                idx2 = sub2ind([nlat nlon nt], ilat, ilons, it+1);
                vallu = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt], ilat   , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat   , ilone , it+1);
                valru = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilons , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilons , it+1);
                valld =  data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilone , it+1);
                valrd =  data(idx1).*(1-st) + data(idx2).*st;
                
                %interpolate along long
                valu = vallu.*(1-slon) + valru.*slon;
                vald = valld.*(1-slon) + valrd.*slon;
                
                %interpolate along lat
                val = valu.*(1-slat) + vald.*slat;
                
            else %space first % NOTE: consider speed up in case only one time is present, unnecessary operations done
                % interpolate along lon
                valbu = permute(data(ilat   , ilons , it  ).*(1-slon) + data(ilat   , ilone , it  ).*slon,[3 1 2]);
                valau = permute(data(ilat   , ilons , min(it+1,size(data,3))).*(1-slon) + data(ilat   , ilone , min(it+1,size(data,3))).*slon,[3 1 2]);
                valbd = permute(data(ilat+1 , ilons , it  ).*(1-slon) + data(ilat+1 , ilone , it  ).*slon,[3 1 2]);
                valad = permute(data(ilat+1 , ilons , min(it+1,size(data,3))).*(1-slon) + data(ilat+1 , ilone , min(it+1,size(data,3))).*slon,[3 1 2]);
                
                %interpolate along lat
                valb = valbd.*(1-slat) + valbu.*slat;
                vala = valad.*(1-slat) + valau.*slat;
                
                %interpolate along time
                val = valb.*(1-st) + vala.*st;
            end
            
        end
        
        
        function [val] = cubicSpline(t)
            % Compute matrix entry for cubic spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),4);
            val(:,1) = (1 - t).^3/6;
            val(:,2) = ((2-t).^3 - 4*(1-t).^3)/6;
            val(:,3) = ((1+t).^3 - 4*(t).^3)/6;
            val(:,4) = (t).^3/6;
        end
        
        function [val] = linearSpline(t)
            % Compute matrix entry for linear spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),2);
            val(:,1) = 1 -t;
            val(:,2) = t;
        end
        
        function [idx, val] = hemisphereCubicSpline(n_az, n_el, az, el)
            % give the index of the hemisphere spline idx
            % first the equator then the first parallel then the second
            % parallel
            %
            % SYNTAX:
            %  [idx] = hemisphereSpline(n_az,n_el,az,el)
            el_step = 90/n_el;
            idx_el = repmat(ceil(el / el_step),1,4);
            idx_el(:,2) = idx_el(:,2) + 1;
            idx_el(:,3) = idx_el(:,3) + 2;
            idx_el(:,4) = idx_el(:,4) + 3;
            az_step = 360/n_el;
            idx_az = repmat(ceil(az / az_step),1,4);
            idx_az(:,2) = idx_az(:,2) + 1;
            idx_az(:,3) = idx_az(:,3) + 2;
            idx_az(:,4) = idx_az(:,4) + 3;
            idx_az(idx_az > n_az) = idx_az(idx_az > n_az) - n_az;
            idx = idx_az + (idx_el -1).*n_el;
            idx = [idx idx+1 idx+2 idx+3];
            
            t_el = rem(el/el_step);
            t_az = rem(az/az_step);
            val_el = Core_Utils.cubicSpline(t_el);
            val_az = Core_Utils.cubicSpline(t_az);
            
            val = [val_az.*repmat(val_el(:,1),1,4) val_az.*repmat(val_el(:,2),1,4) val_az.*repmat(val_el(:,3),1,4) val_az.*repmat(val_el(:,4),1,4)];
        end
        
        function x = despline(x,spline_base)
            % despline signal x
            %
            % SYNTAX:
            %   x = despline(x,<spline_base>)
            if nargin <2
                spline_base = round(size(x,1)/7);
            end
            for i = 1: size(x,2)
                x(:,2) = x(:,2) - splinerMat(1:length(x(:,2)),x(:,2),spline_base);
            end
        end
        
        function [az_grid, el_grid] = getPolarGrid(step_az, step_el)
            % Get a regularly spaced lnots coordinates for a polar grid
            %
            % INPUT
            %   step_az     step of the grid in azimuth   [deg]
            %   step_el     step of the grid in elevation [deg]
            %
            % OUTPUT 
            %   el_grid     grid centers in azimuth   [deg]
            %   el_grid     grid centers in elevation [deg]
            %
            % SYNTAX
            %   [az_grid, el_grid] = Core_Utils.getPolarGrid(step_az, step_el)
            [el_grid, az_grid] = getGrid([step_el step_az], 0, 90, -180, 180);
        end
            
        function id_ok = polarCleaner(az, el, data, step_deg)
            % Remove observations above 3 sigma 
            %
            % SYNTAX
            %   id_ok = Core_Utils.polarCleaner(az, el, data, step_deg)
            
            % az -180 : 180
            % el 0 : 90
            az_grid = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
            el_grid = flipud(((step_deg(end) / 2) : step_deg(end) : 90 - (step_deg(end) / 2))' .* (pi/180));
            n_az = numel(az_grid);
            n_el = numel(el_grid);
            
            % Find map indexes
            col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
            row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
            
            uid = row + (col-1) * numel(el_grid);
            id_ok = true(size(data));
            for b = unique(uid)'
                id_set = find(uid == b);
                dset = data(id_set);
                sigma = std(dset);
                id_ok(id_set(abs(dset) > 3 * sigma)) = false;
            end
        end
        
        function [data_map, n_map, az_grid_out, el_grid_out] = polarGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells, n_min)
            % Grid points on a regularly gridded semi sphere
            %
            % INPUT 
            %   az      azimuth
            %   el      elevation
            %
            % SYNTAX
            %   [data_map, n_map, az_grid, el_grid] = Core_Utils.polarGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells)
            
            % Define grid
            % az -180 : 180
            % el 0 : 90
            el_grid = flipud(((step_deg(end) / 2) : step_deg(end) : 90 - (step_deg(end) / 2))' .* (pi/180));
            flag_congurent_cells = nargin >= 6 && ~isempty(flag_congurent_cells) && flag_congurent_cells;
            if nargin < 7
                n_min = 0; % by default grid all the data
            end            
            if flag_congurent_cells
                step_az = 360 ./ round((360 / step_deg(1)) * cos(el_grid));
                az_grid = {};
                for i = 1 : numel(step_az)
                   az_grid{i} = ((-180 + (step_az(i) / 2)) : step_az(i) : (180 - step_az(i) / 2)) .* (pi/180);
                   n_az(i) = numel(az_grid{i});
                end                
            else
                az_grid = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                n_az = numel(az_grid);
            end
            n_el = numel(el_grid);
            
            % Find map indexes
            row = max(1, min(floor((pi/2 - el) / (step_deg(end) / 180 * pi)) + 1, length(el_grid)));
            if flag_congurent_cells
                col = max(1, min(floor((az + pi) ./ (step_az(row) / 180 * pi) ) + 1, n_az(i)));
            else
                col = max(1, min(floor((az + pi) / (step_deg(1) / 180 * pi) ) + 1, length(az_grid)));
            end
            
            % init maps
            [n_map, data_map] = deal(zeros(n_el, max(n_az)));
                        
            % fill maps
            for i = 1 : numel(data)
                n_map(row(i), col(i)) = n_map(row(i), col(i)) + 1;
                data_map(row(i), col(i)) = data_map(row(i), col(i)) + data(i);
            end
            data_map(n_map > 0) = data_map(n_map > 0) ./ (n_map(n_map > 0));
            
            % zeros cells with a minimum number of data < n_min
            n_map = nan2zero(n_map);
            data_map(n_map <= n_min) = 0;

            flag_debug = false;
            if flag_congurent_cells && flag_debug
                data_congruent = {};
                for i = 1 : numel(el_grid)
                    data_congruent{i} = data_map(i, 1 : numel(az_grid{i}));
                end
                Core_Utils.plotSphPatchGrid(el_grid, az_grid, data_congruent);
            end
            
            % distort map
            if (nargin >= 5) && ~isempty(step_deg_out)
                data_map_in = data_map;
                decl_n = ((pi/2 - el_grid)/(pi/2));
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg(1) / 2)) : step_deg(1) : (180 - step_deg(1) / 2)) .* (pi/180);
                el_grid_out = el_grid;
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                
                if flag_congurent_cells
                    % Interpolate elevation by elevation
                    [az, el] = deal(zeros(n_el, max(n_az)));
                    [az_mg, el_mg] = meshgrid(az_grid_out, el_grid_out);
                    for i = 1 : numel(el_grid)
                        az_tmp = [az_grid{i} nan(1, max(n_az) - n_az(i))];
                        az(i, :) = az_tmp;
                        el(i, :) = el_grid(i);
                        
                        if sum(n_map(i, :) > 0) < 2
                            data_map(i, :) = data_map_in(i,1);
                        else
                            az_tmp = az_tmp(n_map(i, :) > 0)';
                            az_tmp = [az_tmp-2*pi; az_tmp; az_tmp+2*pi];
                            data_tmp = data_map_in(i, n_map(i, :) > 0)';
                            data_tmp = [data_tmp; data_tmp; data_tmp];
                            data_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'linear');
                            az_tmp = [az_grid{i}'-2*pi; az_grid{i}'; az_grid{i}'+2*pi];
                            data_tmp = [n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'; n_map(i, 1 : n_az(i))'];
                            n_map(i, :) = interp1(az_tmp, data_tmp, az_mg(i, :)', 'nearest');
                        end
                    end
                    data_map(n_map == 0) = 0;
                    n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    data_map_in = data_map;
                    
                    % get polar coordinates
                    x = sin(az_grid_out) .* decl_n;
                    y = cos(az_grid_out) .* decl_n;
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                else
                    % get polar coordinates
                    x = sin(az_grid) .* decl_n;
                    y = cos(az_grid) .* decl_n;
                    
                    data_map_in(n_map == 0) = 0;
                    n_map(n_map == 0) = 0.1; % this is to cheat the next scatteredInterpolant
                    
                    funGridder = scatteredInterpolant(x(n_map > 0), y(n_map > 0), data_map_in(n_map > 0), 'linear' );
                    x = -1 : 0.005 : 1;
                    y = x;
                    [x_mg, y_mg] = meshgrid(x, y);
                    polar_data = nan(numel(x), numel(y));
                    id_ok = hypot(x_mg, y_mg) < 1;
                    polar_data(id_ok) = funGridder(x_mg(id_ok), y_mg(id_ok));
                    
                    % Prepare polar gridder
                    funGridder = scatteredInterpolant(x_mg(id_ok), y_mg(id_ok), polar_data(id_ok), 'linear');
                end
                
                
                % Define output grid
                az_grid_out = ((-180 + (step_deg_out(1) / 2)) : step_deg_out(1) : (180 - step_deg_out(1) / 2)) .* (pi/180);
                el_grid_out = flipud(((step_deg_out(end) / 2) : step_deg_out(end) : 90 - (step_deg_out(end) / 2))' .* (pi/180));
                n_az_out = numel(az_grid_out);
                n_el_out = numel(el_grid_out);
                data_map = zeros(n_el_out, n_az_out);
                
                % Get polar coordinates
                decl_n = ((pi/2 - el_grid_out)/(pi/2));
                x = sin(az_grid_out) .* decl_n;
                y = cos(az_grid_out) .* decl_n;
                
                data_map(:) = funGridder(x(:),y(:));
            else
                el_grid_out = el_grid;
                az_grid_out = az_grid;
            end
        end
        
        
        function [z, l, m] = getAllZernike(l_max, m_max, az, r)
            % Generate all the Zernike parameters combinations
            %
            % INPUT 
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   [z, l, m] = getAllZernike(l_max, az, r)
            %   [z, l, m] = getAllZernike(l_max, m_max, az, r)
            if nargin == 3
                r = az;
                az = m_max;
                m_max = l_max;
            end
            
            n_par = l_max * (l_max + 3) / 2 + 1;
            
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
            
            z = Core_Utils.getZernike(l, m, az, r);
        end
        
        function [n] = getAllNormCoeffZernike(l_max, m_max)
            % Generate all the Zernike normalization coefficient
            %
            % SINTAX
            %   [n] = getAllNormCoeffZernike(l_max, m_max)
            n_par = l_max * (l_max + 3) / 2 + 1;
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
            n = sqrt((1+(m~=0)).*(l+1)/pi);
        end
        
        function [z, l, m] = getAllZernikeNorm(l_max, m_max, az, r)
            % Generate all the Zernike parameters combinations
            %
            % INPUT
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   [z, l, m] = getAllZernike(l_max, az, el)
            %   [z, l, m] = getAllZernike(l_max, m_max, az, el)
            if nargin == 3
                r = az;
                az = m_max;
                m_max = l_max;
            end
            
            n_par = l_max * (l_max + 3) / 2 + 1;
            
            l = zeros(n_par, 1);
            m = zeros(n_par, 1);
            i = 0;
            for degree = 0 : l_max
                i = i(end) + (1 : degree + 1);
                l(i) = degree;
                m(i) = -degree : 2 : degree;
            end
            
            l(abs(m) > m_max) = [];
            m(abs(m) > m_max) = [];
            
            z = Core_Utils.getZernikeNorm(l, m, az, r);
        end
        
        function z = getZernikeNorm(l, m, az, r)
            % Get Zernike values for the polynomials
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   z = getZernike(l, m, az, el)
            
            z = zernfun(l, m, r(:), az(:), 'norm');
        end
        
        function z = getZernike(l, m, az, r)
            % Get Zernike values for the polynomials
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       radius [0..1]
            %
            % SINTAX
            %   z = getZernike(l, m, az, r)
                                    
            z = zernfun(l, m, r(:), az(:));
        end
                
        function [S] = reorthZernikeMask(lat, lon, el_thrsh, n_sample)
            % get an orthogonal basis of zernike function that maximize the
            % power in the area of interest
            %
            % SYNTAX:
            % [S] = Core_Utils.reorthZernikeMask(lat,lon,el_thrsh)
            if nargin < 4
                n_sample = 500000;
            end
            %             N = zeros(741,741);
            %             for j = 1 :10
            xy = (rand(n_sample*2,2)-0.5)*2;
            len= sqrt(xy(:,1).^2 + xy(:,2).^2);
            xy(len > 1,:) = [];
            len(len > 1) = [];
            r = len;
            theta = atan2(xy(:,2),xy(:,1));
            el = pi/2 - r*pi/2;
            cc = Constellation_Collector();
            xy(el < (el_thrsh-(2/180*pi)),:) = []; %remove under treshold
            [mask_north, mask_sud] = cc.getGPS.getPolarMask(lat,lon,5);
            
            xy_mask_north = [sin(mask_north(:,1)).*(1 - mask_north(:,2)/(pi/2)) cos(mask_north(:,1)).*(1 - mask_north(:,2)/(pi/2)) ];
            xy_mask_north = [xy_mask_north; xy_mask_north(1,:)];
            idx_inside = true(size(xy,1),1);
            for i = 1 : (size(xy_mask_north,1)-1)
                idx_inside = idx_inside & ((xy_mask_north(i+1,1) - xy_mask_north(i,1)) * (xy(:,2) - xy_mask_north(i,2)) - (xy_mask_north(i+1,2) - xy_mask_north(i,2)) * (xy(:,1) - xy_mask_north(i,1))) > 0;
            end
            xy(idx_inside,:) = [];
            
            len= sqrt(xy(:,1).^2 + xy(:,2).^2);
            r = len;
            theta = atan2(xy(:,2),xy(:,1));
            el = pi/2 - r*pi/2;
            
            [z1] = Core_Utils.getAllZernike(10, 10, theta, r);
            %z1 = [zernfun(1,1,r,theta,'norm')      zernfun(1,-1,r,theta,'norm') zernfun(5,5,r,theta,'norm')];%           2/sqrt(pi)
            %z1 = [r .* sin(theta)      r .* cos(theta)];
            %       1   -1    r * sin(theta)                 2/sqrt(pi)
            
            %             N = N + z1'*z1;
            %             j
            %             end
            N = z1'*z1;
            [U,D,V] = svd(N);
            d = diag(D);
            % get the flexum
            exclude_fisrt = round(length(d)/5);
            [power_smooth] = splinerMat(exclude_fisrt:length(d),d(exclude_fisrt:end),20);
            [~,flex] = max(Core_Utils.diffAndPred(power_smooth,3));
            
            stop = exclude_fisrt + flex;
            
            
            S = U(:,1:stop);
        end
        
        function [z_interp, l, m] = zSinthesysAll(l_max, m_max, az, r, z_par)
            % Get Zernike interpolation given the coefficients
            % of their polynomials
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   z_par   Zernike coefficients
            %
            % SINTAX
            %   [z_interp] = zSinthesysAll(l_max, m_max, az, el, z_par)
            z_interp = nan(size(az));
            
            id_ok = ~isnan(az);
            [A, l, m] = Core_Utils.getAllZernike(l_max, m_max, az, r);
            z_interp(id_ok) = A * z_par;
        end
        
        function [z_interp, l, m] = zSinthesys(l, m, az, r, z_par)
            % Get Zernike interpolation given the coefficients
            % of their polynomials
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   z_par   Zernike coefficients
            %
            % SINTAX
            %   [z_interp] = zSinthesys(l, m, az, el, z_par)
            z_interp = nan(size(az));
            
            id_ok = ~isnan(az);
            A = Core_Utils.getZernike(l, m, az, r);
            z_interp(id_ok) = A * z_par;
        end
        
        function [z_par, l, m, A] = zAnalisysAll(l_max, m_max, az, r, data, max_reg)
            % Get Zernike polynomials parameters 
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized
            %   max_reg maximum regularization
            %
            % SINTAX
            %   [z_par, l, m, A] = zAnalisysAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:));
            [A, l, m] = Core_Utils.getAllZernike(l_max, m_max, az(id_ok), r(id_ok));
            if nargin == 6 && ~isempty(max_reg) && max_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
            end
            z_par = (A'*A + diag(reg_fun .* ones(size(A, 2), 1))) \ (A' * data(id_ok));    
        end
        
         
        function [z_par, l, m, A] = zroAnalisysAll(l_max, m_max, az, r, data, S)
            % Get Zernike polynomials parameters form reothonormalize
            % function
            %
            % SINTAX
            %   [z_par, l, m, A] = zAnalisysAll(l_max, m_max, az, el, data, S)
            
            id_ok = ~isnan(data(:));
            [A, l, m] = Core_Utils.getAllZernike(l_max, m_max, az(id_ok), r(id_ok));
            N = A'*A;
            B = A'*data;
            N = S'*N*S;
            B = S'*B;
            z_par = N\B;
            z_par = S*z_par;
        end
        
        function [z_par, l, m, A] = zAnalisys(l, m, az, r, data, max_reg)
            % Get Zernike polynomials parameters 
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized
            %   max_reg maximum regularization
            %
            % SINTAX
            %   [z_par, l, m, A] = zAnalisysAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:));
            A = Core_Utils.getZernike(l, m, az(id_ok), r(id_ok));
            if nargin == 6 && ~isempty(flag_reg) && flag_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
            end
            z_par = (A' * A + reg_fun .* diag(ones(size(A, 2), 1))) \ (A' * data(id_ok));
        end
                
        function [filtered_data, z_par, l, m,  A] = zFilter(l_max, m_max, az, r, data, max_reg)
            % Get Zernike polynomials parameters 
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized
            %   max_reg maximum regularization
            %
            % SINTAX
            %   [z_par, A] = zFilter(az, el, data, <max_reg = 1>)
            
            if nargin == 6
                [z_par, l, m, A] = Core_Utils.zAnalisysAll(l_max, m_max, az, r, data, max_reg);
            else
                [z_par, l, m, A] = Core_Utils.zAnalisysAll(l_max, m_max, az, r, data);
            end
            filtered_data = A * z_par;
        end

        function flag_is_hold = isHold()
            % Return if there is an open figure on hold
            %
            % SYNTAX
            %   Core_Utils.isHold()
            
            flag_is_hold = ~isempty(findobj('Type', 'figure')) && ishold;
        end
        
        function fh = showZernike(l, m, z_par, el_min, funMapElevation)
            % Show 3D plot of Zernike polynomials 
            %
            % SINTAX
            %   fh = showZernike(l, m, z_par)
            
            %if ~Core_Utils.isHold()
                is_hold = false;
                fh = figure();
            %else
            %    is_hold = true;
            %    fh = gcf;
            %end
            %%% INTERNAL PARAMETER
            scale = 1;
            %%%

            x = -1 : 0.005 : 1;
            y = x;
            [X,Y] = meshgrid(x,y);
            [theta, r_prj] = cart2pol(X,Y); % This radius is the correct one for my polar projection             
            if nargin == 5
                r_prj = funMapElevation(r_prj * (pi/2)) / (pi/2);
            end
            r_zern = r_prj;
            if nargin >= 4 && ~isempty(el_min)
                r_max = 1 - (2 * el_min / pi);
                idx = r_prj <= r_max;
            else
                idx = r_prj <= 1;
            end
            z = nan(size(X));      
            z(idx) = zernfun(l, m, r_zern(idx), 2 * pi - theta(idx) + pi/2) * z_par;
            
            %h = scatter(X(~isnan(z)),Y(~isnan(z)),160,z(~isnan(z)),'filled');
            h = imagesc(x,y,z);
            h.AlphaData = ~isnan(z);
            ax = gca;
            ax.YDir = 'normal';
            plot_bg = ~is_hold;
            if plot_bg
                hold on
                %plot parallel
                az_l = [0:pi/200:2*pi];
                d_step = 15/180*pi;
                decl_s = ([0:d_step:pi/2]/(pi/2))*scale;
                for d = decl_s
                    x = cos(az_l).*d;
                    y = sin(az_l).*d;
                    plot(x,y,'color',[0.6 0.6 0.6]);                    
                    text(cos(80/180*pi)*d,sin(80/180*pi)*d,sprintf('%d',round(d*90)),'HorizontalAlignment','center', 'FontWeight', 'bold');
                end
                %plot meridian
                az_step = 30/180 *pi;
                az_s = [0:az_step:2*pi];
                decl_l = ([0 1])*scale;
                for a = az_s
                    x = cos(a).*decl_l;
                    y = sin(a).*decl_l;
                    plot(x,y,'color',[0.6 0.6 0.6]);
                    if abs(a-2*pi) > 0.0001
                        text(cos(a)*1.1,sin(a)*1.1,sprintf('%d', mod(round((2*pi - a + pi/2) / pi * 180), 360)), 'HorizontalAlignment','center', 'FontWeight', 'bold');
                    end
                end
                axis equal
                axis off
                set(gcf,'color','w');
                if ~is_hold
                    hold off
                end
                xlim([-1.15 1.15]); ylim([-1.15 1.15]);
                colormap(jet);
                colorbar;
            end
            fh = gcf; Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'light');            
        end
        
        function fh = showZernike3(l, m, z_par, el_min)
            % Show 3D plot of Zernike polynomials 
            %
            % SINTAX
            %   fh = showZernike3(l, m, z_par)
            
            % [x, y] = pol2cart(theta, r_synt);
            if nargin == 4
                r_max = 1 - (2 * el_min / pi);
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, r_max, 101));
            else
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, 1, 101));
            end
            % r_prj is the correct radius for my polar projection 
            r_zern = r_prj;
    
            z = nan(size(theta));
            z(:) = zernfun(l, m, r_zern(:), theta(:)) * z_par;
            %z = r_zern.*cos(theta);
            
            
            fh = figure();
            title('Zernike expansion')
            %polarplot3d(z, 'PlotType','surfn');
            polarplot3d(z,'PlotType','surfn','PolarGrid',{4 24}, 'PolarDirection', 'cw', 'TickSpacing',5,...
                   'RadLabels',4, 'RadLabelLocation',{180 'max'}, 'RadLabelColor','black', 'AxisLocation', 'mean');
            ar = get(gca,'DataAspectRatio');
            set(gca, 'DataAspectRatio', [1 1 ar(3)]);
            
            colormap(jet);
            material([ 0.4 0.9 0.55])
            l1 = light('position',[200 -300 400], 'color', [0.6 0.6 0.6]);
            l2 = light('position',[-600 600 900], 'color', [0.6 0.6 0.6]);
            l3 = light('position',[0 2 100], 'color', [0.6 0.6 0.6]);
            view(35, 45);
            Core_UI.beautifyFig(fh, 'dark');            
        end
        
        function fh = showZernike3StylePCV(l, m, z_par, el_min, limits)
            % Show 3D plot of Zernike polynomials 
            %
            % INPUT
            %   l       zernike degree [n x 1]
            %   m       zernike degree [n x 1]
            %   z_par   zernike degree [n x 1]
            %   el_min  cut-off angle  [rad]
            %   limits  min max (saturation) of the z map
            %
            % SINTAX
            %   fh = showZernike3(l, m, z_par, <el_min>, <limits>)
            
            % [x, y] = pol2cart(theta, r_synt);
            if nargin >= 4 && ~isempty(el_min)
                r_max = 1 - (2 * el_min / pi);
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, r_max, 101));
            else
                [theta, r_prj] = meshgrid(linspace(0, 2*pi, 361), linspace(0, 1, 101));
            end
            % r_prj is the correct radius for my polar projection 
            r_zern = r_prj;
    
            z = nan(size(theta));
            z(:) = zernfun(l, m, r_zern(:), theta(:)) * z_par;
            
            fh = figure();
            title('Zernike expansion')
            
            if nargin == 5 && ~isempty(limits)
                z = min(limits(2), max(z, limits(1)));
            end
            polarplot3d(z, 'PlotType', 'surf', 'RadialRange',[0 90] / 180 * pi, ...
                    'AxisLocation', 0, 'InterpMethod', 'cubic', ...
                    'PlotType', 'surfn', 'tickspacing', 15, ...
                    'GridColor', [0.7 0.7 0.7]);
                        
            colormap(flipud(Cmap.get('PuOr', 256)));
            
            axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38], ...
                'Xlim', 1.5 * [-90 90] / 180 * pi, 'Ylim', 1.5 * [-90 90] / 180 * pi, ...
                % 'XTick', [], 'YTick', [], 'Color', 'none', 'XColor', 'none', 'YColor', 'none'
                };
            ax = gca;
            set(ax, axprop{:});
            
            if nargin >= 5 && ~isempty(limits)
                caxis(limits);
                zlim(limits);
            end

            
            %Core_UI.beautifyFig(fh, 'dark');
        end
        
        function fh = polarZerMap(l_max, m_max, az, r, data)
            % Take scattered observation and plot a polar interpolation
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMap(l_max, m_max, az, r, data)
            [z_par, l, m] = Core_Utils.zAnalisysAll(l_max, m_max, az, r, data, 1);
            fh = Core_Utils.showZernike(l, m, z_par);
        end
        
        function fh = polarZerMapDual(l_max, m_max, az, r, data)
            % Take scattered observation and plot:
            %  - a scattered polar
            %  - a polar interpolation
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMapDual(l_max, m_max, az, el, data)
            fh = figure('Visible', 'off'); subplot(1,2,1); hold on;
            Core_Utils.polarZerMap(l_max, m_max, az, r, data);
            fh.Visible = false;
            cax = caxis();
            subplot(1,2,2); hold on; 
            polarScatter(az, (r * pi/2), 50, data, 'filled'); colorbar(); colormap(jet);
            cax_s = caxis();
            %cax = [min(cax(1), cax_s(1)) max(cax(2), cax_s(2))];
            cax = cax_s;
            caxis(cax);
            subplot(1,2,1);
            caxis(cax);
            fh.Visible = true;
            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
        end
        
        function fh = polarZerMapQuad(l_max, m_max, az, r, data)
            % Take scattered observation and plot:
            %  - a scattered polar
            %  - a polar interpolation
            %  - a scattered vs el
            %  - an interpolation vs el
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %   data    list of data to be analized            
            %
            % SYNTAX
            %   fh = polarZerMapQuad(l_max, m_max, az, el, data)
            if ~isempty(findobj('Type', 'figure')) && ~ishold
                fh = gcf;
            else
                fh = figure('Visible', 'on'); 
            end
            
            subplot(2,2,1); hold on;
            Core_Utils.polarZerMap(l_max, m_max, az, r, data);
            cax = caxis();
            
            subplot(2,2,2); hold on; 
            polarScatter(az, r * pi/2, 50, data, 'filled'); colorbar(); colormap(jet);
            cax_s = caxis();
            %cax = [min(cax(1), cax_s(1)) max(cax(2), cax_s(2))];
            cax = cax_s;
            caxis(cax);
            subplot(2,2,1);
            caxis(cax);
            
            subplot(2,2,3); 
            data_smooth = Core_Utils.zFilter(l_max, m_max, az, r, data); hold on;       
            plot((r * pi/2) / pi * 180, data_smooth, '.', 'Color', Core_UI.getColor(1));
            
            el_grid = 0 : 90;
            plot(el_grid, Core_Utils.interp1LS((r * pi/2) / pi * 180, data_smooth, 3, el_grid), '--', 'LineWidth', 2, 'Color', Core_UI.getColor(3));
            
            ylim(cax);
            xlim([0 90]);
            grid minor;

            subplot(2,2,4); hold on; 
            plot((r * pi/2) / pi * 180, data, '.', 'Color', Core_UI.getColor(2)); hold on;
            plot(el_grid, Core_Utils.interp1LS((r * pi/2) / pi * 180, data, 3, el_grid), '--', 'LineWidth', 2 , 'Color', Core_UI.getColor(5));
            ylim(cax);
            xlim([0 90]);
            grid minor;
            
            fh.Visible = true;
            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
        end
        
                
        function sphTest()
            %%  along track analysis            
            l_min = 0;
            l_max = 30;
            id = 1:20:length(data);
            [az el] = meshgrid(lambdaGrid, phiGrid);
            [ N_subSet, TN_subSet ] = analisiPointC_defineSYS (data(:), 1+0*data(:), el(:)./180*pi, az(:)./180*pi, el(:)*0+1, l_min, l_max, l_min, l_max, 1, 1, 0);
            tic;
            idM = 1:(l_max+1)^4;
            idM = reshape(idM,(l_max+1)^2,(l_max+1)^2);
            idM = tril(idM); idM = idM(idM>0);
            N = zeros((l_max+1)^2,(l_max+1)^2);
            N(idM) = N_subSet;
            N = tril(N) + tril(N)' + diag(diag(N));
            
            cLS = zeros(l_max+1);
            sLS = zeros(l_max+1);
            x = N\TN_subSet;
            i = 0;
            for az = 1 : l_max +1
                for m = 1 : az
                    i = i + 1;
                    cLS(az,m) = x(i);
                    if m > 1
                        i = i + 1;
                        sLS(az,m) = x(i);
                    end
                end
            end
            
            [ topo_map2 ] = sintesiGrid (phiGrid./180*pi, lambdaGrid./180*pi, cLS, sLS, 0, l_max, 0, l_max, 1, 1, 0, 0);
        end
        
        function plm = fplm(l, m, theta)
            % Computing Legendre polynomial
            %
            % INPUT
            %   l       degree
            %   m       order
            %   theta   theta angle [rad]
            %
            % SYNTAX
            %   plm = Core_Utils.fplm(l, m, decl)
            % 
            lMin = l(1);
            lMax = l(end);
            mMin = m(1);
            mMax = m(end);
            
            if (size(l,1)==1)
                l=l';
            end
            if (size(m,1)==1)
                m=m';
            end
            if (size(theta,1)==1)
                theta=theta';
            end
            
            %%
            r1 = zeros(mMax,1);
            
            % computing root 1 (bl)
            r1(1) = sqrt(3);
            i = (2:mMax);
            l = (1:lMax);
            
            r1(2:end) = sqrt((2*i+1)./(2*i));
            %%
            % Init P (result matrix)
            
            plm = zeros(mMax,length(m),length(theta));
            
            % Computing Pmm
            
            % Skip the calculous of the first l(1)-1 lines
            % Go from 0 to lMin-1
            Ptmp0 = ones(length(theta),1);
            for mfix = 1:mMin-1
                Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
                Ptmp0 = Ptmp1;
            end
            
            % for each m
            r2 = zeros(length(l)-1,1);
            r3 = zeros(length(l)-1,1);
            
            for mfix = mMin:mMax
                % Computing Pmm --------------------------------------------------
                Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
                Ptmp0 = Ptmp1;
                
                % Save in the results matrix the Pmm element
                plm(mfix, mfix-mMin+1, :) = Ptmp1(:);
                
                % Computing Plm --------------------------------------------------
                
                % get the row
                r = mfix+1;
                
                
                % computing root 2 (clm)
                r2 = sqrt(((2*l(r:end)+1).*(2*l(r:end)-1))./((l(r:end) + mfix).*(l(r:end) - mfix)));
                
                % computing root 3 (dlm)
                r3 = sqrt(((2*l(r:end)+1).*(l(r:end)+mfix-1).*(l(r:end)-mfix-1))./((2*l(r:end)-3).*(l(r:end)+mfix).*(l(r:end)-mfix)));
                
                Pl1 = Ptmp1;
                Pl2 = zeros(size(Ptmp1));
                for lfix = r:lMax
                    tmp = r2(lfix-r+1) * cos(theta).*Pl1 - r3(lfix-r+1).*Pl2;
                    plm(lfix,mfix-mMin+1, :) = tmp(:);
                    Pl2 = Pl1;
                    Pl1 = tmp(:);
                end
            end
            
            plm = plm(lMin:lMax,:,:);
        end

        %--------------------------------------------------------------------------
        % TRIGONOMETRIC manipulators
        %--------------------------------------------------------------------------

        function r_angle = deg2rad(d_angle)
            % Convert degrees to radians
            %
            % SYNTAX
            %   r_angle =Core_Utils.deg2rad(d_angle)
            r_angle = d_angle .* (pi/180);
        end
        
        function d_angle = rad2deg(r_angle)
            % Convert degrees to radians
            %
            % SYNTAX
            %   r_angle =Core_Utils.deg2rad(d_angle)
            d_angle = r_angle .* (180/pi);
        end
        
        %--------------------------------------------------------------------------
        % DATA manipulators
        %--------------------------------------------------------------------------
        
        function [t, data_set] = insertNan4Plots(t, data_set)
            % Insert a Nan in a regularly sampled dataset to make
            % plots interrupt continuous lines
            %
            % INPUT
            %   t      epoch of the data [matrix of column arrays]
            %   data   epoch of the data [matrix of column arrays]
            %
            % SYNTAX
            %   [t, data] = Core_Utils.insertNan4Plots(t, data)
                        
            t = t(:);
            if size(t, 1) ~= size(data_set, 1)
                % data should be a column array
                data_set = data_set';
            end
            n_set = size(data_set, 2);
            dt = diff(t);
            rate = median(dt);
            id_in = find(dt > 1.5 * rate);
            for x = numel(id_in) : -1 : 1
                t = [t(1 : id_in(x)); (t(id_in(x)) + 1.5 * rate); t((id_in(x)+1) : end)];
                data_set = [data_set(1 : id_in(x), :); nan(1, n_set); data_set((id_in(x)+1) : end, :)];
            end
        end
        
        function lh = plotSep(t, data, varargin)
            % Special wrapper to regular plot
            % Works on regularly sampled data
            % When there is a gap of data, it insert a nan value
            % to avoid linear interpolation of the data
            %
            % INPUT
            %   t           column array of epochs
            %   data        columns of data (could be a matrix)
            %   varagin     add other useful parameters of the plot
            %
            % SYNTAX
            %   lh = Core_Utils.plotSep(t, data, varagin);
            %
            % SEE ALSO
            %   plot
            [t, data] = Core_Utils.insertNan4Plots(t, data);
            if nargin <= 2
                lh = plot(t, data);
            else
                lh = plot(t, data, varargin{:});
            end
        end
        
        %--------------------------------------------------------------------------
        % OTHER FUNCTIONS
        %--------------------------------------------------------------------------
        
        function printEx(ex)
            % Print exception to screen
            %
            % SYNTAX:
            %  Core_Utils.printEx(ex)
            
            msg = sprintf(['goGPS encountered a problem, please open an issue on GitHub posting\n' ...
                    ' the following lines together with a copy of the full log' ]);
            fprintf('\n---------------------------------------------------------------------\n %s', msg);
            fprintf('\n---------------------------------------------------------------------\n MESSAGE: %s\n---------------------------------------------------------------------\n\n', ex.message);
            
            for i=1:numel(ex.stack)
                fprintf('  file: "%s"\n  line: %d\n  fun: %s\n\n', ex.stack(i).file, ex.stack(i).line, ex.stack(i).name);
            end
            % keyboard
        end
        
        function exportCurFig(out_path, mode)
            % Eport Current Figure
            %
            % EXAMPLE
            %   Core_Utilis.exportCurFig(fullfile('/Users/Andrea/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/GIMS/title.png'))
            if nargin == 2
                Core_Utils.exportFig(gcf, out_path, mode);
            else
                Core_Utils.exportFig(gcf, out_path);
            end
        end
        
        function exportFig(fh, out_path, mode)
            % Export Figure
            %
            % SYNTAX
            %   Core_Utils.exportFig(<fh>, out_path)
            %   Core_Utils.exportFig(out_path)
            %
            % EXAMPLE
            %   Core_Utilis.exportFig(fh, fullfile('/Users/Andrea/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/GIMS/title.png'))
            if nargin <= 1 
                out_path = fh;
                fh = gcf;
            end
            %fh.WindowStyle = 'normal'; export_fig(fh, out_path, '-transparent', '-r150'); fh.WindowStyle = 'docked';
            ws_bk = fh.WindowStyle;
            fh.WindowStyle = 'normal';
            if nargin == 3 && ~isempty(mode)
                Core_UI.beautifyFig(fh, mode);
            end
            col = fh.Color;
            Logger.getInstance.addMessage(sprintf('Exporting to "%s"', out_path));
            box = findall(fh, 'type', 'uicontainer');
            if isempty(box)
                export_fig(fh, out_path, '-transparent', '-r150');
            else
                % Special tricks in case of figure containing boxes
                
                % Use saveas instead of export_fig
                [~, ~, ext] = fileparts(out_path);
                if strcmp(ext, '.png')
                    bg_color = [1 255 1]; % Use green screen
                else
                    bg_color = [255 255 255];
                end
                % fallback
                bg_box = {};
                for b = 1: numel(box)
                    try
                        bg_box{b} = box(b).BackgroundColor;
                        box(b).BackgroundColor = bg_color/255;
                    catch ex
                    end
                end
                saveas(fh, out_path);
                for b = 1: numel(box)
                    try                        
                        box(b).BackgroundColor = bg_box{b};
                    catch ex
                    end
                end
                
                % saveas does not have transparency management
                if strcmp(ext, '.png')
                    % Tricks are just for PNG file type

                    % read the just saved image 
                    [im_out] = imread(out_path);
                    % convert it in hue saturation value
                    im_hsv = rgb2hsv(im_out);
                    bg_hsv = rgb2hsv(bg_color/255);
                    alpha_mask = (im_hsv(:,:,1) ~= bg_hsv(1));
                    color_mask = (~(im_out(:,:,1) == bg_color(1) & im_out(:,:,2) == bg_color(2) & im_out(:,:,3) == bg_color(3)));
                    saturation = im_hsv(:,:,2);
                    value = im_hsv(:,:,3);
                    hue = im_hsv(:,:,1);
                    saturation(alpha_mask ~= 1) = 0;
                    hue(alpha_mask ~= 1) = 1;
                    alpha_mask = double(alpha_mask);
                    
                    % get alpha value for the value component of the image
                    alpha_mask(logical(color_mask - alpha_mask)) = double(1-value(logical(color_mask - alpha_mask)));
                    
                    im_hsv(:,:,1) = hue;
                    im_hsv(:,:,2) = saturation;
                    im_hsv(:,:,3) = value;
                    im_out = hsv2rgb(im_hsv);
                    
                    imwrite(im_out, out_path, 'alpha', alpha_mask);
                end               
            end
                
            fh.WindowStyle = ws_bk;
            fh.Color = col;
        end
        
        function export2Plotly(fh, flag_offline, new_rate)
            % Export figure on plotly managing time axes
            %
            % SYNTAX
            %   export2Plotly(fig_handle, flag_offline, new_rate);
            %
            % NOTE 
            %   to use this function you must install plotly API:
            %       https://plot.ly/matlab/getting-started/#installation
            
            
            if nargin < 3
                new_rate = [];
            end
            
            if nargin < 2 || isempty(flag_offline)
                flag_offline = true;
            end
            
            if nargin > 1 && ~isempty(fh)
                figure(fh);
            else
                fh = gcf;
            end
            
            % Prepare times for plotly
            ax = findall(fh, 'type', 'axes');
            last = struct();
            last.mode = {};
            for i = 1 : numel(ax)
                last.mode{i} = ax(i).XTickLabelMode;
                ax(i).XTickLabelMode = 'auto';
            end
            
            line = findall(fh, 'type', 'line');
            last.flag_date = [];
            for i = 1 : numel(line)
                if ~isempty(new_rate)
                    time = line(i).XData * 86400;
                    rate = median(round(diff(time), 3));
                    [~, id_ok] = ismember((unique(round(time ./ new_rate)) .* new_rate), (round(time ./ (2*rate)) .* (2*rate)));
                    line(i).XData = line(i).XData(noZero(id_ok));
                    line(i).YData = line(i).YData(noZero(id_ok));
                end
                
                if all(line(i).XData > datenum('1970-01-01') & line(i).XData < datenum('2070-01-01'))
                    line(i).XData = convertDate(line(i).XData);
                    last.flag_date(i) = true;
                else
                    last.flag_date(i) = false;
                end
            end            
            
            % Export to plotly
            try
                fig_name = fh.UserData.fig_name;
            catch
                fig_name = 'temp_unknown';
            end
            plotlyfig = fig2plotly(gcf, 'filename', fig_name, 'offline', flag_offline);

            if any(last.flag_date)
                for i = 1 : numel(ax)
                    plotlyfig.layout.(sprintf('xaxis%d', i)).type = 'date';
                end
                plotlyfig.PlotOptions.FileOpt = 'overwrite';
                plotly(plotlyfig);
            end
            
            for i = 1 : numel(line)
                if last.flag_date(i)
                    funToMatTime = @(date) (date/(1000*60*60*24) + datenum(1969,12,31,19,00,00));
                    line(i).XData = funToMatTime(line(i).XData);
                end
            end
            
            % Restore time for matlab
            for i = 1 : numel(ax)
                ax(i).XTickLabelMode = last.mode{i};
                setTimeTicks(ax(i));
            end
            
            Core_UI.beautifyFig(fh);
        end
                
        function idx = findMO(find_list, to_find_el)
            % find the postion of the elements of to_find_el into find_list
            % find list should have unique elements
            idx = zeros(size(to_find_el));
            for i = 1: length(to_find_el)
                idx(i) = find(find_list == to_find_el(i),1);
            end
        end
        
        function [amb_idx, n_amb] = remEmptyAmbIdx(amb_idx, n_amb)
            % remove emtpy amb_idx
            %
            % SYNTAX:
            % amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx, <n_amb>)
            if nargin < 2
                n_amb = max(max(amb_idx));
            end
            i = 1;
            while (i <= n_amb)
                n_ep_amb = sum(sum(amb_idx == i));
                if n_ep_amb == 0
                    n_amb = n_amb - 1;
                    amb_idx(amb_idx > i) = amb_idx(amb_idx > i) - 1;
                else
                    i = i + 1;
                end
            end
        end
        
        function num = round_even(num)
            num = round((num-2)/2)*2+2;
        end
        
        function num = round_odd(num)
            num = round((num-1)/2)*2+1;
        end
        
        function r = xcorr(x)
            % compute cross correlation
            %
            % SYNTAX:
            % xcorr = Core_Utils.xcorr(x)
            %
            % NOTE:
            % thank you Amro https://stackoverflow.com/questions/3949324/calculate-autocorrelation-using-fft-in-matlab
            len = length(x);
            
            % autocorrelation
            nfft = 2^nextpow2(2*len-1);
            r = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );
            % rearrange and keep values corresponding to lags: -(len-1):+(len-1)
            r = [r(end-len+2:end) ; r(1:len)];
        end
        
        function [s,n] = getSemivariogram1D(x,mode)
            % compute 1 d semivariogram
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            if nargin < 2
                mode = 'mean';
            end
            max_lag = length(x)-1;
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            if strcmpi(mode,'mean')
                for l = 1 : max_lag
                    diffs = (x((l+1):end) - x(1:(end-l))).^2;
                    s(l) = mean(diffs,'omitnan')/2;
                    n(l) = sum(~isnan(diffs));
                end
            elseif strcmpi(mode,'fft')
                L = length(x);
                %                 Y = fft(x)
                %                 Ae = abs(Y);
                A = Core_Utils.getSpectrum(x);
                cova = ifft(([A; flipud(A(2:end))]).^2);
                s = var(x) - cova(1:101);
            else
                for l = 1 : max_lag
                    s(l) = median((x((l+1):end) - x(1:(end-l))).^2,'omitnan')/2;
                end
            end
        end
        
        function [s,n] = getSemivariogram1DCS(x,cs_lid)
            % compute 1 d semivariogram accoutnign for cycle slips
            %
            % SYNTAX:
            %     s = Core_Utils.getSemivariogram1D(x)
            cs = unique([1; find(cs_lid); length(x)]);
            max_lag = length(x)-1;
            
            s = nan(max_lag,1);
            n = zeros(max_lag,1);
            for c = 2 : length(cs)
                if sum(~isnan(x(cs(c-1):cs(c)))) > 1
                    [st,nt] = Core_Utils.getSemivariogram1D(x(cs(c-1):(cs(c)-1)));
                    ntt = n(1:length(st))+nt;
                    s(1:length(st)) = zero2nan( ( nan2zero(s(1:length(st))).*n(1:length(st)) + nan2zero(st).*nt )./ntt );
                    n(1:length(st)) = ntt;
                end
            end
        end
        
        function [a,b] = logLogLineEst(y,lims)
            % compute slope and intercept in log log plane
            %
            % SYNTAX:
            %     [a,b] = Core_Utils.logLogLineEst(y,lims)
            x = (1:length(y))';
            x = x(lims(1):lims(2));
            y = y(lims(1):lims(2));
            A = [log(x) ones(size(x))];
            y = log(y);
            est = A\y;
            a=est(1);
            b=est(2);
        end
        
        function [gh11,nh11] =variofl(x1,icode)
            % SOURCE : Marcotte, Denis. "Fast variogram computation with FFT." Computers & Geosciences 22.10 (1996): 1175-1186.
            %
            % function [gh11,nh11l=variofl(x1,icode);
            %
            % function to compute variograms or covariograms, in 1D or 2D
            % the data are on a (possibly incomplete) regular grid.
            % the program computes variograms in the frequency domain by
            % using 2D-FFT.
            %
            % input: x1: data matrix. Missing values are indicated by NaN
            %
            %icode: a code to indicate which function to compute
            %=l : variogram
            %
            % =2 : covariogram
            %
            %
            % gh11: variogram or covariogram depending on icode.
            % output:
            % nh11: number of pairs available
            %
            %
            % this program uses the functions FFT2, IFFTZ, FFTlSHIFT and CONJ which are
            % standard MATLAB functions.
            [n,p]=size(x1);
            nrows=2*n-1;
            ncols=2*p-1;
            % dimensions of data matrix
            % find the closest multiple of 8 to obtain a good compromise between
            % speed (a power of 2) and memory required
            approx = 8;
            nr2=ceil(nrows/approx)*approx;
            nc2=ceil(ncols/approx)*approx;
            % form an indicator matrix:
            %      l's for all data values
            %      O's for missing values
            %
            % in data matrix, replace missing values by 0;
            x1id=~isnan(x1);
            x1(~x1id)=zeros(sum(sum(-x1id)),1); % 1 for a data value; 0 for missing
            % missing replaced by 0
            fx1=fft2(x1,nr2,nc2); % fourier transform of x1
            if icode==1
                fx1_x1=fft2(x1.*x1,nr2,nc2);
            end
            clear x1;
            fx1id=fft2(x1id,nr2,nc2);
            clear x1id
            % fourier transform of x1*x1
            % fourier transform of the indicator matrix
            % compute number of pairs at all lags
            nh11=round(real(ifft2(conj(fx1id).*fx1id)));
            % compute the different structural functions according to icode
            if icode==1
                % variogram is computed
                gh11=real(ifft2(conj(fx1id).*fx1_x1+conj(fx1_x1).*fx1id-2*conj(fx1).*fx1));
                gh11=gh11./max(nh11,1)/2;
            else
                % covariogram is computed
                ml=real(ifft2(conj(fx1).*fx1id))./max(nh11,1);
                m2=real(ifft2(conj(fx1id).*fx1))./max(nh11,1);
                clear fx1id
                gh11=real(ifft2(conj(fx1).*fx1));
                gh11=gh11./max(nh11,1)-ml.*m2;
            end
            % compute tail mean
            % compute head mean
            clear fx1 fx1id fx1_fx1
            % reduce matrix to required size and shift so that the 0 lag appears at the center of each matrix
            nh11 = [nh11(1:n,1:p) nh11(1:n,nc2-p+2:nc2); nh11(nr2-n+2,1:p) nh11(nr2-n+2:nr2,nc2-p+2:nc2)];
            gh11 = [gh11(1:n,1:p) gh11(1:n,nc2-p+2:nc2); gh11(nr2-n+2,1:p) gh11(nr2-n+2:nr2,nc2-p+2:nc2)];
            
            %             gh11=fftshift(gh11);
            %             nh11=fftshift(nh11);
        end
        
        function s = remPeriod(s,p)
            % remove signal for a given period
            %
            % SYNTAX:
            %     s = Core_Utils.remPeriod(s,p)
            ph = (1:numel(s))'/p*2*pi;
            A = [cos(ph) sin(ph)];
            N = A'*A; % less stable but faster
            B = A'*s;
            x = N\B;
            s = s - A*x;
        end
        
        function num = code2Char2Num(str2)
            % Convert a 2 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str2(:,1:2) * [2^8 1]';
        end
        
        function str2 = num2Code2Char(num)
            % Convert a numeric value (float) of a 2 char string
            % SYNTAX
            %   str2 = Core_Utils.num2Code2ch(num)
            num = double(num);
            str2 = char(zeros(numel(num), 2));
            str2(:,1) = char(floor(num / 2^8));
            num = num - str2(:,1) * 2^8;
            str2(:,2) = char(num);
        end
        
        function num = code3Char2Num(str3)
            % Convert a 3 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str3(:,1:3) * [2^16 2^8 1]';
        end
        
        function str3 = num2Code3Char(num)
            % Convert a numeric value (float) of a 3 char string
            % SYNTAX
            %   str3 = Core_Utils.num2Code3ch(num)
            num = double(num);
            str3 = char(zeros(numel(num), 3));
            str3(:,1) = char(floor(num / 2^16));
            num = num - str3(:,1) * 2^16;
            str3(:,2) = char(floor(num / 2^8));
            num = num - str3(:,2) * 2^8;
            str3(:,3) = char(num);
        end
        
        function num = code4Char2Num(str4)
            % Convert a 4 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code4ch2Num(str4);
            
            num = str4(:,1:4) * [2^24 2^16 2^8 1]';
        end
        
        function str4 = num2Code4Char(num)
            % Convert a numeric value (float) of a 4 char string
            % SYNTAX
            %   str4 = Core_Utils.num2Code4Char(num)
            str4 = char(zeros(numel(num), 4));
            str4(:,1) = char(floor(num / 2^24));
            num = num - str4(:,1) * 2^24;
            str4(:,2) = char(floor(num / 2^16));
            num = num - str4(:,2) * 2^16;
            str4(:,3) = char(floor(num / 2^8));
            num = num - str4(:,3) * 2^8;
            str4(:,4) = char(num);
        end
        
        function str4 = unique4ch(str4)
            % Perform unique on an array of 4 char codes
            %
            % SYNTAX
            %   str4 = Core_Utilis.unique4ch(str4)
            str4 = Core_Utils.num2Code4Char(unique(Core_Utils.code4Char2Num(str4)));
        end
        
        function str3 = unique3ch(str3)
            % Perform unique on an array of 3 char codes
            %
            % SYNTAX
            %   str3 = Core_Utilis.unique3ch(str3)
            str3 = Core_Utils.num2Code3Char(unique(Core_Utils.code3Char2Num(str3)));
        end
        
        function str2 = unique2ch(str2)
            % Perform unique on an array of 2 char codes
            %
            % SYNTAX
            %   str2 = Core_Utilis.unique2ch(str3)
            str2 = Core_Utils.num2Code2Char(unique(Core_Utils.code2Char2Num(str2)));
        end
        
        function f_status_lst = aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, date_list, out_dir)
            % Try to download files using aria2C
            %
            % INPUT
            %   file_name_list      list of file_names to download (remote path)   [cell]
            %   f_ext_lst           extension of compression ('' is valid)         [cell]
            %   f_status_lst        bool array of files existing                   [bool]
            %   date_list           GPS_Time of days of interest                   [GPS_Time]
            %   out_dir             path to out folder                             [char]
            %
            % SYNTAX
            %   f_status_lst = Core_Utils.aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, <date_list>, <out_dir>)
            %
            
            log = Core.getLogger();
            fnp = File_Name_Processor();
            rm = Remote_Resource_Manager.getInstance;
            
            if ispc()
                aria2c_path = '.\utility\thirdParty\aria2-extra\aria2_win\aria2c.exe';
            elseif ismac()
                aria2c_path = '/usr/local/bin/aria2c';
            else % is linux
                aria2c_path = '/usr/local/bin/aria2c';
                if ~exist(aria2c_path, 'file')
                    aria2c_path = '/usr/bin/aria2c';
                end
            end
            
            if ~exist(aria2c_path, 'file')
                Core.getLogger.addError('aria2c is not working, is it installed?');
                f_status_lst = [];
                return
            end

            % Get file list path for all the files present in the list
            idf = find(~f_status_lst);
            fnl = file_name_lst(idf);
            fel = f_ext_lst(idf); % file extension list
            odl = {}; % out_dir_list
            ffp = {}; % final file path
            for i = 1 : numel(fnl)
                file_name = fnl{i};
                server = regexp(file_name,'(?<=\?{)\w*(?=})','match', 'once'); % saerch for ?{server_name} in paths
                if ~isempty(server)
                    file_name = strrep(file_name,['?{' server '}'],'');
                    [s_ip, port] = rm.getServerIp(server);
                    switch port
                        case '21'
                            fnl{i} = ['ftp://' s_ip ':' port file_name fel{i}];
                        otherwise
                            fnl{i} = ['http://' s_ip ':' port file_name fel{i}];
                    end
                else
                    fnl{i} = [file_name fel{i}];
                end
                if nargin < 5
                    out_dir = Core.getState.getFileDir(file_name);
                end
                if nargin >= 4 && ~isempty(date_list)
                    out_dir = fnp.dateKeyRep(out_dir, date_list.getEpoch(date_list.length - idf(i) + 1));
                end
                odl{i} = out_dir;
                [~, name, ext] = fileparts(fnl{i});
                if strcmp(ext,'.Z') || strcmp(ext,'.gz')
                    ffp{i} = fullfile(out_dir, name);
                else
                    ffp{i} = fullfile(out_dir, [name ext]);
                end
            end
            
            % if I have at least one file to download
            if numel(odl) > 0
                i = 0;
                file_name = fullfile('.', 'reserved', 'tmpAriaDownload.lst');
                if (exist(['.' filesep 'reserved' filesep], 'dir') == 0)
                    mkdir(['.' filesep 'reserved']);
                end
                fid = fopen(file_name, 'Wb');
                if fid < 0
                    log.addWarning(['Writing on "' file_name '" is not possible' char(10) 'aria2 could not work']);
                else
                    str = '';
                    old_od = odl{1};
                    while i <= numel(fnl)
                        i = i + 1;
                        if (i < numel(fnl)) && strcmp(odl{i}, old_od)
                            str = sprintf('%s%s\n', str, fnl{i});
                        else
                            
                            if i <= numel(fnl)
                                if i == numel(fnl) && strcmp(odl{i}, old_od)
                                    str = sprintf('%s%s\n', str, fnl{i});
                                end
                                
                                % call aria
                                % check for .Z or .gz compressed files too
                                fwrite(fid, str, '*char');
                                fclose(fid);
                                if ~isempty(str)
                                    if ~exist(old_od, 'file')
                                        mkdir(old_od);
                                    end
                                    log.addMessage(sprintf('Executing \n  aria2c -c -i %s -d %s\n  File download list:', file_name, old_od));
                                    log.addMessage(log.indent(sprintf('%s', str)));
                                    try
                                        if ispc()
                                            dos(sprintf('"%s" -j 20 -c -i %s -d %s >nul 2>&1', aria2c_path, file_name, old_od)); % suppress output
                                            % dos(sprintf('"%s" -j 20 -c -i %s -d %s', aria2c_path, file_name, old_od)); % do not suppress output
                                        else
                                            dos(sprintf('%s -j 20 -c -i %s -d %s &> /dev/null', aria2c_path, file_name, old_od));  % suppress output
                                            %dos(sprintf('%s -j 20 -c -i %s -d %s', aria2c_path, file_name, old_od));  % do not suppress output
                                        end
                                    catch
                                        this.log.addError('aria2c is not working, is it installed?');
                                    end
                                end
                                % Check for zero byte files (download errors)
                                for f = 1 : numel(fnl)
                                    [~, out_file_name, out_file_ext] = fileparts(fnl{f});
                                    out_file_path = [old_od, filesep, out_file_name, out_file_ext];
                                    if exist(out_file_path, 'file') == 2
                                        file_info = dir(out_file_path);
                                        if file_info.bytes == 0
                                            % the file is empty
                                            delete(out_file_path)
                                            log.addError(sprintf('%s download failed\nThe file is probably missing', [out_file_name, out_file_ext]));
                                        end
                                    end
                                end
                                % open file list for the next set
                                fid = fopen(file_name, 'Wb');
                                str = '';
                            end
                        end
                        if i <= numel(fnl)
                            decrement = 0;
                            if ~(strcmp(odl{i}, old_od))
                                decrement = 1; % this last file have not been downloaded!
                            end
                            old_od = odl{i};
                            i = i - decrement;
                        end
                    end
                    fclose(fid);
                    delete(file_name);
                end
                
                % once the file have been downloaded decompress and test presence
                for i = 1 : numel(fnl)
                    if strcmp(fel{i},'.Z') || strcmp(fel{i},'.gz')
                        if (isunix())
                            system(['gzip -d -f ' ffp{i} fel{i} '&> /dev/null']);
                        else
                            try
                                [status, result] = system(['.\utility\thirdParty\7z1602-extra\7za.exe -y x '  ffp{i} fel{i} ' -o'  odl{i} ]); %#ok<ASGLU>
                                if (status == 0)
                                    status = true;
                                end
                                delete([ffp{i} fel{i}]);
                            catch
                                this.log.addError(sprintf('Please decompress the %s file before trying to use it in goGPS!!!', fname));
                                status = false;
                            end
                        end
                    end
                end
                
                % update f_status_lst
                for i = 1 : numel(ffp)
                    f_status_lst(idf(i)) = exist(ffp{i}, 'file') == 2;
                end
            end
            
        end
        
        function [status] = downloadHttpTxtResUncompress(filename, out_dir)
            log = Core.getLogger();
            fnp = File_Name_Processor();
            try
                options = weboptions;
                options.ContentType = 'text';
                options.Timeout = 15;
                [remote_location, filename, ext] = fileparts(filename);
                filename = [filename ext];
                log.addMessage(log.indent(sprintf('downloading %s ...',filename)));
                compressed_name = '';
                status = true;
                if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                    mkdir(out_dir);
                end
                try
                    txt = websave(fullfile(out_dir, filename), ['http://' remote_location '/' filename]);
                catch ex
                    if instr(ex.message, '404')
                        try
                            compressed_name = [filename, '.gz'];
                            txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name]);
                        catch ex
                            if instr(ex.message, '404')
                                try
                                    compressed_name = [filename, '.Z'];
                                    txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name]);
                                catch
                                    status = false;
                                end
                            end
                        end
                    end
                end
                if status
                    status = false; %#ok<NASGU>
                    if ~isempty(compressed_name)
                        compressed_name = fnp.checkPath(fullfile(out_dir, compressed_name));
                        if (isunix())
                            system(['gzip -d -f ' compressed_name '&> /dev/null &']);
                        else
                            try
                                [status, result] = system(['.\utility\thirdParty\7z1602-extra\7za.exe -y x '  compressed_name ' -o'  out_dir ]); %#ok<ASGLU>
                                if (status == 0)
                                    status = true;
                                end
                                delete(compressed_name);
                            catch
                                this.log.addError(sprintf('Please decompress the %s file before trying to use it in goGPS!!!', compressed_name));
                                status = false;
                            end
                        end
                    end
                    status = true;
                    log.addMessage(' Done');
                end
            catch
                status = false;
            end
        end
        
        function [status, ext] = checkHttpTxtRes(filename)
            ext = '';
            if isunix()
                if ismac()
                    % Sometimes mac users does not have wget'
                    rem_check_cmd = 'curl --head ';
                else
                    % curl seems to have some problems with matlb libraries
                    % under some Linux installations, switching to wget --spyder
                    rem_check_cmd = 'wget --spider ';
                end
                
                [resp, txt] = system([rem_check_cmd filename]);
                if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) %#ok<STREMP>
                    status = true;
                else
                    [resp, txt] = system([rem_check_cmd filename '.gz']);
                    if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) %#ok<STREMP>
                        ext = '.gz';
                        status = true;
                    else
                        [resp, txt] = system([rem_check_cmd filename '.Z']);
                        if ~isempty(strfind(txt,' 200 OK')) || ~isempty(strfind(txt,' 302 Found')) %#ok<STREMP>
                            ext = '.Z';
                            status = true;
                        else
                            status = false;
                        end
                    end
                end
            else
                log = Logger.getInstance;
                log.addWarning('HTTP check is implemeted only for Unix systems')
                status = true; % !!! to be implemented
            end
        end
        
        function station_list = getStationList(dir_path_list, file_ext, flag_recursive)
            % Get the list of stations present in a folder (with keys substituted)
            %
            % SYNTAX
            %   station_list = Core_Utilis.getStationList(dir_path)
            
            if nargin == 1
                file_ext = '.';
            else
                file_ext = ['[' file_ext ']'];
            end
            
            if nargin == 3 && flag_recursive
                dir_path_list = strsplit(genpath(dir_path_list), ':');
            else
                dir_path_list = {dir_path_list};
            end
            
            recursive_dirs = {};
            for d = 1 : numel(dir_path_list)
                dir_path = dir_path_list{d};
                if ~isempty(dir_path)
                    try
                        % Calling dos is faster than dir with large directories
                        if isunix
                            [~, d] = dos(['ls ' dir_path]);
                            dir_list = strsplit(d);
                            dir_list = dir_list(1:end-1);
                        else
                            [~, d] = dos(['dir ' dir_path]);
                            dir_list = strsplit(d);
                            dir_list = dir_list(1:end-1);
                        end
                    catch
                        dir_list = dir(dir_path);
                        dir_list = {dir_list.name};
                    end
                    for i = 1 : numel(dir_list)
                        dir_list{i} = [dir_path(length(dir_path_list{1})+2:end) filesep dir_list{i}];
                    end
                end
                recursive_dirs = [recursive_dirs dir_list];
            end
            dir_list = recursive_dirs;
            %%
            % search for station files STAT${DOY}${S}${QQ}.${YY}
            file_list = {};
            for d = 1 : numel(dir_list)
                tmp = strsplit(dir_list{d}, filesep);
                file_name = tmp{end};
                file_name_len = numel(file_name);
                rin2_start = regexp(dir_list{d}, ['.{4}[0-9]{3}.{1}[0-9]{2}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once');
                rin3_start = regexp(dir_list{d}, '\_[0-9]{4}[0-9]{3}[0-9]{4}\_', 'once');
                if (file_name_len == 14) && ~isempty(rin2_start)
                    file_list = [file_list; {[dir_list{d}(1:rin2_start) '${DOY}${S}${QQ}.${YY}' dir_list{d}(end)]}];
                elseif ~isempty(rin3_start)
                    file_list = [file_list; {[dir_list{d}(1:rin3_start) '${YYYY}${DOY}' dir_list{d}(rin3_start + 8 : end)]}]; %#ok<AGROW>
                end
            end
            
            % search for station files STAT${DOY}${S}.${YY}
            for d = 1 : numel(dir_list)
                tmp = strsplit(dir_list{d}, filesep);
                file_name = tmp{end};
                file_name_len = numel(file_name);
                rin2_start = regexp(dir_list{d}, ['[0-9]{3}.{1}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once');
                if (file_name_len == 12) && ~isempty(rin2_start)
                    file_list = [file_list; {[dir_list{d}(1:rin2_start-1) '${DOY}${S}.${YY}' dir_list{d}(end)]}]; %#ok<AGROW>
                    %file_list = [file_list; dir_list{d}(1:4)];
                end
            end
            station_list = unique(file_list);
        end
        
        function data = injectData(data1, data2, idx1, idx2, data_size)
            % isert data2 into data1 at the position definied by idx1 and idx2
            % idx1 - 1 is the last element of data1 to be putted before data2 (0  if none)
            % idx2 + 1 is the first element of data1 to be put after data2 (data1 length if none)
            %
            % SYNTAX
            %   data = Core_Utils.injectData(data1, data2, idx1, idx2)
            if nargin == 5
                % check data size
                [m, n] = size(data2);
                if not(m == data_size(1) && n == data_size(2))
                    data2 = nan(data_size(1), data_size(2));
                end
            end
            if ~isempty(data2)
                if size(data1, 1) < idx2
                    data1((end + 1) : idx2, :) = nan;
                end
                data = [data1(1 : idx1 - 1, :); data2; data1(idx2 + 1 : end, :)];
            else
                data = data1;
            end
        end
        
        function data = injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time1, time2, id_stop, id_start, interpolate)
            % inject smoothed data
            %
            % INPUT:
            %   data_lft : all data left
            %   data_right : all data right
            %   idx_smt1 : which data of data_lft are to be smoothed
            %   idx_smt2 : which data of data_rgt are to be smoothed
            %   time_1: time of the left data to be smoothed
            %   time_2: time of the right data to be smoothed
            %   id_stop: first epoch of time_1 that should not be kept
            %   id_start: first epoch of time_2 to keep
            %
            % SYNTAX:
            %   data = Core_Utils.injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time_1, time_2, id_start)
            if nargin < 9
                interpolate = true;
            end
            
            % Get the part of data to interp
            data_tosmt_lft = data_lft(idx_smt1);
            data_tosmt_rgt = data_rgt(idx_smt2);
            
            % we use mat time, is easier and we do not need extreme precision
            time1 = time1.getMatlabTime();
            time2 = time2.getMatlabTime();
            
            [idx1, idx2, time_tot] = Core_Utils.intersectOrderedDouble(time1, time2, median([diff(time1); diff(time2)])/4); % 1/4 the rate tolerance
            
            mix_len = min(0.007, abs((time2(1) - time1(end)))/20); % <= empirically found
            w2 = 1 ./ (1 + exp(-((time_tot - mean(time_tot)) / mix_len))); % todo: scale to ensure [0 1]
            w1 = 1 - w2;
            n_out = size(time_tot);
            data1 = nan(n_out);
            data2 = nan(n_out);
            data1(idx1) = data_tosmt_lft;
            data2(idx2) = data_tosmt_rgt;
            %id_start = idx1(id_start);
            %id_ko = ((isnan(data1) & (1 : n_out)' < id_start) |(isnan(data2) & (1 : n_out)' >= id_start)) & ~(isnan(data1) & isnan(data2)); %?? should be used -> yes beacuse time is injected deleting overlapping times
            id_keep = unique([idx1(1:(id_stop-1)); idx2(id_start:end)]);
            % Interpolate missing data
            if interpolate
                is_nan = find(isnan(data1));
                extr_lft = is_nan(time_tot(is_nan) <= min(time1));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time1));
                data1(extr_lft) = data_tosmt_lft(1);
                data1(extr_rgh) = data_tosmt_lft(end);
                if any(~isnan(data_tosmt_lft)) &&  any(isnan(data1))
                    data1 = simpleFill1D(data1, isnan(data1), 'linear', time_tot);
                end
                is_nan = find(isnan(data2));
                extr_lft = is_nan(time_tot(is_nan) <= min(time2));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time2));
                data2(extr_lft) = data_tosmt_rgt(1);
                data2(extr_rgh) = data_tosmt_rgt(end);
                if any(~isnan(data_tosmt_rgt)) && any(isnan(data2))
                    data2 = simpleFill1D(data2, isnan(data2), 'linear', time_tot);
                end
            else
                data1(isnan(data1)) = data2(isnan(data1));
                data2(isnan(data2)) = data1(isnan(data2));
            end
            
            % Do not Merge nans
            last_ok = find(~isnan(data1), 1, 'last');
            if isempty(last_ok)
                last_ok = 0;
            end
            w1((last_ok + 1) : end) = 0;
            data1((last_ok + 1) : end) = 0;
            w2((last_ok + 1) : end) = 1;
            
            % Do not Merge nans
            first_ok = find(~isnan(data2), 1, 'first');
            if isempty(last_ok)
                first_ok = numel(data2) + 1;
            end
            w1(1 : (first_ok - 1)) = 1;
            w2(1 : (first_ok - 1)) = 0;
            data2(1 : (first_ok - 1)) = 0;
            
            data = w1.*data1 + w2.*data2;
            %data(id_ko) = [];
            data = data(id_keep);
            data = [data_lft(~idx_smt1); data; data_rgt(~idx_smt2)];
        end
        
        function [idx1, idx2, double_tot] = intersectOrderedDouble(double_1, double_2, threshold)
            % given two ordered double give the index of the two vector in the joint vector considering the threshold
            %
            % SYNTAX
            % [idx1, idx2] = Core_Utils.intersectOrderedDouble(double_1, double_2, threshold)
            l1 = length(double_1);
            l2 = length(double_2);
            idx1 = zeros(l1,1);
            idx2 = zeros(l2,1);
            i = 1;
            j = 1;
            tot = 1;
            while  i <=  l1 && j <= l2
                if abs(double_1(i) - double_2(j)) < threshold
                    idx1(i) = tot;
                    idx2(j) = tot;
                    i = i + 1;
                    j = j + 1;
                    tot = tot + 1;
                elseif double_1(i) < double_2(j)
                    idx1(i) = tot;
                    i = i + 1;
                    tot = tot + 1;
                else
                    idx2(j) = tot;
                    j = j + 1;
                    tot = tot +1;
                end
            end
            if j > l2 && i <= l1
                idx_end = (i : l1) -l1 + tot;
                idx1(i : l1) = idx_end;
            elseif i > l1 && j <= l2
                idx_end = (j : l2) -l2 + tot;
                idx2(j : l2) = idx_end;
            end
            double_tot = zeros(max(max(idx1), max(idx2)), 1);
            double_tot(idx1) = double_1;
            double_tot(idx2) = double_2;
        end
        
        function [ids] = findAinB(cellA,cellB)
            % find the index of cella in cellb
            %
            % SYNTAX
            %     [ids] = Core_Utils.findAinB(cellA,cellB)
            if ~iscell(cellA)
                cellA = {cellA};
            end
            lB = length(cellB);
            lA = length(cellA);
            ids = zeros(lA,1);
            for i = 1 : lA
                not_found = true;
                j = 1;
                while j <= lB && not_found
                    if strcmp(num2str(cellA{i}),num2str(cellB{j}))
                        ids(i) = j;
                        not_found = false;
                    end
                    j = j+1;
                end
            end
            %ids(ids==0) = [];
        end
        
        function [wl_cyle_out, frac_bias] = getFracBias(wl_cycle, weigth)
            % get the common frac bias between cycles
            % NOTE: very coarse/nobrain/empirical solution - > a simpler one should be found
            %
            % SYNTAX
            %   [wl_cyle, frac_bias] = Core_Utils.getFracBias(wl_cycle)
            if nargin < 2
                weigth = ones(size(wl_cycle));
            end
            frac_bias = zeros(1,3);
            wl_cycle_frac = zeros(size(wl_cycle,1),3);
            % get receiver wsb
            wl_cycle_frac(:,1) = zero2nan(wl_cycle) - floor(zero2nan(wl_cycle));
            wl_cycle_frac(:,2) = zero2nan(wl_cycle) - round(zero2nan(wl_cycle));
            wl_cycle_frac(:,3) = zero2nan(wl_cycle) - ceil(zero2nan(wl_cycle));
            frac_bias(1) = median(wl_cycle_frac(:,1),'omitnan');
            frac_bias(2) = median(wl_cycle_frac(:,2),'omitnan');
            frac_bias(3) = median(wl_cycle_frac(:,3),'omitnan');
            wl_cyle_var = mean(abs(wl_cycle_frac-repmat(frac_bias,size(wl_cycle_frac,1),1)),'omitnan');
            [~,idx] = min(wl_cyle_var);
            frac_bias = frac_bias(idx);
            wl_cycle_frac = wl_cycle_frac(:,idx);
            a = 0;
            idx_rw = ones(size(wl_cycle_frac));
            while sum(zero2nan(wl_cycle_frac - frac_bias) < -0.5) > 0  || a < 4
                wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) = wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) + 1;
                idx_rw(:) = 1;
                e = abs((wl_cycle_frac-frac_bias)) > 0.2;
                idx_rw(e) = 1./((wl_cycle_frac(e)-frac_bias)/0.2).^2; %idx_reweight
                idx_rw = idx_rw  .* weigth;
                idx_rw = idx_rw / sum(idx_rw);
                frac_bias = sum((wl_cycle_frac).*idx_rw);
                a = a+1;
            end
            wl_cyle_out = wl_cycle - frac_bias;
        end
        
        function [response] = timeIntersect(time1_st, time1_end, time2_st, time2_end)
            % check whether times of time bound 1 intersect with times of time bound 2
            %
            % SYNTAX
            % [response] = Core_Utils.timeIntersect(time1_st,time1_en, time2_st, time2_en)
            response = time1_st <= time2_end & time1_end >= time2_st;
            %(time2_st <= time1_st & time2_end >= time1_st) | (time2_st <= time1_en & time2_end >= time1_en);
        end
                
        function [ it, st, ilons, ilone, slon, ilat, slat] = getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % get  interpolating index
            [nlat , nlon, nt] = size(data);
            
            
            lon(lon < first_lon) = lon(lon < first_lon) + nlon * dlon; %% to account for earth circularity
            % find indexes and interpolating length
            % time
            it = max(min(floor((t - first_t)/ dt)+1,nt-1),1);
            st = max(min(t - first_t - (it-1)*dt, dt), 0) / dt;
            st = serialize(st);
            
            % lat
            ilat = max(min(floor((lat - first_lat)/ dlat)+1,nlat-1),1);
            slat = min(max(lat - first_lat - (ilat-1)*dlat, dlat), 0) / dlat;
            
            % lon
            ilons = max(min(floor((lon - first_lon)/ dlon)+1,nlon),1);
            ilone = ilons +1;
            ilone(ilone > nlon) = 1;
            slon = max(min(lon - first_lon- (ilons-1)*dlon, dlon), 0) / dlon;
        end
        
        function response = permutedEqual(str1, str2)
            % check if the two variables are permuted version of the same sequence
            %
            % SYNTAX:
            %     response =Core_Utils.pertutedEqual(var1, var2)
            if length(str1) ~= length(str2)
                response = false;
            else
                ll = length(str1);
                found_all = true;
                i = 1;
                j = 1;
                while i <= ll
                    found = false;
                    while j <= ll && ~found
                        found = found || str1(i) == str2(j);
                        j = j+1;
                    end
                    found_all = found_all && found;
                    i = i+1;
                end
                response = found_all;
            end
            
        end
        
        function A = remBFromA(A,B)
            % reomve the lement of B from A
            for i = 1 : length(B)
                A(A==B(i)) = [];
            end
        end
        
        function [iA] = inverseByPartsDiag(A,idx1,idx2)
            % inverse by partitioning of a diagonal block matrix
            %
            % SYNTAX:
            %  [iA] = inverseByPartsDiag(A,idx1,idx2)
            sz1 = sum(idx1);
            sz2 = sum(idx2);
            A11 = A(idx1, idx1);
            A22 = A(idx2, idx2);
            A21 = A(idx2, idx1);
            A12 = A(idx1, idx2);
            iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            r1 = A12*iA22*A21;
            r2 = A21*iA11*A12;
            A21 = -iA22*A21;
            A12 = -iA11*A12;
            A11 = A11 - r1;
            A22 = A22 - r2;
            iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            iA = sparse(size(A));
            iA(idx1,idx1) = iA11;
            iA(idx1,idx2) = A12*iA22;
            iA(idx2,idx1) = A21*iA11;
            iA(idx2,idx2) = iA22;
        end
        
        function [iA] = inverseByRecPartsDiag(A,indices)
            % inverse by recursive partitioning of a diagonal block matrix 
            % 
            % SYNTAX:
            %    [iA] = inverseByPartsDiags(A,indices)
            n_blk = length(indices);
            
            % slit the indece in two
            part1 = 1:ceil(n_blk/2);
            part2 = (ceil(n_blk/2)+1):n_blk;
            if length(part1) > 1
                composite1 = true;
            else
                composite1 = false;
                sz1 = sum(indices{part1});
            end
            if length(part2) > 1
                composite2 = true;
            else
                composite2 = false;
                sz2 = sum(indices{part2});
            end
            
            idx1 = false(size(indices{1}));
            idx2 = false(size(indices{1}));

            for i = 1 : length(part1)
                idx1 = idx1 | indices{part1(i)};
            end
            for i = 1 : length(part2)
                idx2 = idx2 | indices{part2(i)};
            end
            
            % actual inversion
            A11 = A(idx1, idx1);
            A22 = A(idx2, idx2);
            A21 = A(idx2, idx1);
            A12 = A(idx1, idx2);
            if composite1
                new_indeces1 = {};
                for i = 1 : length(part1)
                    idx_tmp  = indices{part1(i)};
                    new_indeces1{i} = idx_tmp(idx1);
                end
                iA11 = Core_Utils.inverseByRecPartsDiag(A11,new_indeces1);
            else
                iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            end
            if composite2
                new_indeces2 = {};
                for i = 1 : length(part2)
                    idx_tmp  = indices{part2(i)};
                    new_indeces2{i} = idx_tmp(idx2);
                end
                iA22 = Core_Utils.inverseByRecPartsDiag(A22,new_indeces2);
            else
                iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            end
            r1 = A12*iA22*A21;
            r2 = A21*iA11*A12;
            A21 = -iA22*A21;
            A12 = -iA11*A12;
            A11 = A11 - r1;
            A22 = A22 - r2;
            if composite1
                iA11 = Core_Utils.inverseByRecPartsDiag(A11,new_indeces1);
            else
                iA11 = spdiags(1./diag(A11), 0, sz1, sz1);
            end
            if composite2
                iA22 = Core_Utils.inverseByRecPartsDiag(A22,new_indeces2);
            else
                iA22 = spdiags(1./diag(A22), 0, sz2, sz2);
            end
            iA = sparse(size(A));
            iA(idx1,idx1) = iA11;
            iA(idx1,idx2) = A12*iA22;
            iA(idx2,idx1) = A21*iA11;
            iA(idx2,idx2) = iA22;
            
        end
        
        function [fb, frac_b_mat]= estimateFracBias(obs_cy, cycle_slip)
            % estimate the common factional bias to all the obesravtions
            %
            % SYNTAX:
            %    fb = Network.estimateFracBias(obs_cy, cycle_slip)
            amb_idx = Core_Utils.getAmbIdx(cycle_slip, obs_cy);
            frac_cy = obs_cy;
            n_arcs = max(amb_idx(:,end));
            frac_b = zeros(n_arcs,1);
            num_ep = zeros(n_arcs,1);
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb),'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = frac_cy(idx_amb) - amb;
                frac_b(i) = median(frac_cy(idx_amb),'omitnan'); %strongMean(frac_cy(idx_amb),1, 0.95,2.5);
                num_ep(i) = sum(sum(idx_amb));
                frac_b_mat(idx_amb) = frac_b(i);
            end
            frac_b_mat = frac_b_mat+0*obs_cy;
            fb = Core_Utils.circularModedRobustMean(frac_b_mat(:));
            frac_cy = nan(size(obs_cy));
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb) - fb,'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = obs_cy(idx_amb) - amb - fb;
                frac_b_mat(idx_amb) = amb;
            end
            fb = fb + strongMean(frac_cy(:),1, 0.90,2.5);
        end
        
        function fr_cy = circularMean(cycle, obs_weigth)
            % estimate the mean for data over 0 -1
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularMean(cycle)
            
            % elimnate the nan
            if nargin < 2
                obs_weigth = ones(size(cycle));
            end
            idx_nan = isnan(cycle);
            cycle(idx_nan) = [];
            obs_weigth(idx_nan) = [];
            cycle = cycle*2*pi;
            % compute the unit vectors, take the mean and recompute the
            % average
            
            unit = [cos(cycle) sin(cycle)];
            obs_weigth = obs_weigth./ sum(obs_weigth);
            mean_vec = [mean(unit(:,1) .* obs_weigth) mean(unit(:,2) .* obs_weigth)];
            fr_cy = atan2(mean_vec(1), mean_vec(2))/(2*pi);
            
            
        end
        
        function plotSphPatchGrid(el_grid, az_grid, data)
            % Plot a 3D spherical cap 
            % ad hoc function for plotting congruent cells
            % see. Generating Fuhrmann et al. 2013 statistically robust multipath stacking maps using congruent cells
            %
            % used for generateMultipath
            %
            % INPUT
            %   el_grid     array of elevation       [1 x n]
            %   az_grid     cell of array of azimuth [1 x n] of [1 x m(e)]            
            %   data        cell of data elevation x elevation in azimuth [1 x n] of [1 x m(e)]
            %   
            % SYNTAX: 
            %   plotSphPatchGrid(el_grid, az_grid, data)
            
            figure;
            hold on
            el_grid = flipud(el_grid);
            az_grid = fliplr(az_grid);
            data = fliplr(data);
            d_el  = diff(el_grid);
            el_grid = el_grid - [d_el; d_el(end)]/2;
            el_grid =  [el_grid; el_grid(end)+d_el(end)];
            for e = 1 : (length(el_grid)-1)
                data{e} = fliplr(data{e});
                d_az  = diff(az_grid{e}');
                if ~isempty(d_az)
                    az_gridw = az_grid{e}'- [d_az; d_az(end)]/2;
                    az_gridw =  [az_gridw; az_gridw(end)+d_az(end)];
                    for a = 1 : (length(az_gridw)-1)
                        plotSphPatch(el_grid(e:e+1)',az_gridw(a:a+1)', data{e}(a));
                    end
                end
            end
            
            function plotSphPatch(lats, lons, color)
                if lons(2) < lons(1)
                    lons(2) = lons(2)+2*pi;
                end
                dlons = (1 + (1 -cos(mean(lats)))*20)/180*pi;
                lonst = (lons(1):dlons:lons(2))';
                if lonst(end) ~= lons(2)
                    lonst = [lonst; lons(2)];
                end
                lons = lonst;
                lats = pi/2 -lats;
                xyz =[ [cos(lons)*sin(lats(1)) sin(lons)*sin(lats(1)) ones(size(lons))*cos(lats(1))];
                    flipud([cos(lons)*sin(lats(2)) sin(lons)*sin(lats(2)) ones(size(lons))*cos(lats(2))])];
                patch(xyz(:,1),xyz(:,2),xyz(:,3),color);
            end
        end
        
        function fr_cy = circularModedRobustMean(cycle)
            % estimate a roubust mean mean for data over 0 -1
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularModedRobustMean(cycle)
            
            % elimnate the nan
            cycle(isnan(cycle)) = [];
            mode_cy = mode(cycle);
            
            % center around the mode and then take the string(robust) mean
            idx_inf = (cycle - mode_cy) < -0.5;
            idx_sup = (cycle - mode_cy) > 0.5;
            cycle(idx_inf) = cycle(idx_inf) +0.5;
            cycle(idx_sup) = cycle(idx_sup) -0.5;
            
            fr_cy = strongMean(cycle,1,0.95,2.5);
            
        end
        
        function amb_idx = getAmbIdx(cycle_slip , obs)
            % get matrix of same dimesion of the observation showing the ambiguity index of the obsarvation
            %
            % SYNTAX:
            % this.getAmbIdx()
            
            amb_idx = ones(size(cycle_slip), 'uint16');
            n_epochs = size(amb_idx,1);
            n_stream = size(amb_idx,2);
            for s = 1:n_stream
                if s > 1
                    amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                end
                cs = find(cycle_slip(:, s) > 0)';
                for c = cs
                    amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                end
            end
            amb_idx = amb_idx .* uint16(obs ~= 0);
            amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx);
        end
        
        function createEmptyProject(base_dir, prj_name, prj_type)
            % create empty config file
            %
            % SYNTAX
            %    createEmptyProject(base_dir, prj_name)
            %    createEmptyProject(prj_name)
            
            fnp = File_Name_Processor();
            
            if nargin == 3
                state = Main_Settings('', fnp.checkPath([base_dir filesep prj_name]));
            else
                state = Main_Settings('');
            end
            
            if nargin == 1
                prj_name = base_dir;
                base_dir = fnp.getFullDirPath([state.getHomeDir filesep '..']);
            end
            
            log = Core.getLogger();
            log.addMarkedMessage(sprintf('Creating a new project "%s" into %s', prj_name, [base_dir filesep prj_name]));
            
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name]));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'config']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'out']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'out/log']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'RINEX']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station/CRD']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station/ocean']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station/MET']));
            state.setPrjHome(fnp.checkPath([base_dir filesep prj_name]));
            state.prj_name = prj_name;
            state.setOutDir('./out');
            state.setObsDir('./RINEX');
            state.setMetDir('./station/MET');
            state.setCrdDir('station/CRD');
            state.setMPDir('antenna/MP');
            if ~exist(state.getMPDir, 'dir')
                try
                    mkdir(state.getMPDir);
                catch ex
                end
            end
                
            crd_file = state.getCrdFile;
            if ~exist(crd_file, 'file')
                try
                    % try to create an empty ocean loading file
                    fid = fopen(crd_file, 'w');
                    fwrite(fid, '# Insert here the coordinates of the receivers');
                    fclose(fid);
                catch ex
                end
            end
            
            state.setOceanDir('station/ocean');
            ocean_file = fullfile(state.getHomeDir, state.getOceanFile);
            if ~exist(ocean_file, 'file')
                try
                    % try to create an empty ocean loading file
                    fid = fopen(ocean_file, 'w');
                    fwrite(fid, '$$ Insert here the ocean loading displacements');
                    fclose(fid);
                catch ex
                end
            end
            
            config_path = fnp.checkPath([base_dir filesep prj_name filesep 'config' filesep 'config.ini']);
            if nargin == 3
                switch prj_type
                    case 1 % PPP for tropo
                        state.setToTropoPPP();
                    case 2 % NET no iono, no tropo
                        state.setToShortNET();
                    case 3 % NET no iono
                        state.setToMediumNET();
                    case 4 % NET iono-free
                        state.setToLongNET();
                end
            end
            state.check();
            state.save(config_path);
            Core.getCurrentSettings.import(state);
        end
        
        function y = fillNan1D(y,x)
            % fill the nan into the y
            %
            % SYNTAX
            % y = Core_Utils.fillNan1D(y,<x>)
            
            if nargin < x
                x = 1: length(y);
            end
            idx_nan = isnan(y);
            if sum(~idx_nan) > 2
                int_data = interp1(x(~idx_nan),y(~idx_nan),x(idx_nan));
                y(idx_nan) = int_data;
            end
        end
        
        function [Amp,Phase,f] = getSpectrum(y,smpl_rate)
            % compute the spectrum with fft
            %
            % SYNTAX:
            %  [Amp,Phase,f] = Core_Utils.getSpectrum(y,smpl_rate);
            if nargin < 2
                smpl_rate = 1;
            end
            Y = fft(y);
            
            L = length(y);
            Fs = 1 /smpl_rate;
            f = Fs*(0:(L/2))/L;
            Y = Y(1:(L/2 +1));
            Amp = abs(Y/L)*2;
            Phase = angle(Y);
            
        end
        
        function [val] = spline(t,order)
            % Compute matrix entry for spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            switch order
                case 1
                    val = Core_Utils.linearSpline(t);
                case 2
                    %%% tBD
                case 3
                    val = Core_Utils.cubicSpline(t);
            end
            
        end
        
        function [id_ko] = snoopGatt(ssat_err, thr, thr_propagate)
            % mark the residual after one thresdol till their mov max  reenter the
            % second threshold
            %
            % SYNTAX:
            %    [w] = Core_Utils.snoopGatt(res, thr, thr_propagate)
            id_ko = false(size(ssat_err));
            for s = 1 : size(id_ko, 2)
                id_ko(:,s) = (movmax(abs(ssat_err(:,s)), 20) > thr_propagate) & flagExpand(abs(ssat_err(:,s)) > thr, 100);
            end
        end
                
        function [x,inv_diag] = fastInvDiag(N,B,mode)
            % solve the linear system and compute the diagonal entry of the inverse of N square matrix. This is
            % much faster than computing the whole inverse in case of very sparse
            % matrix. Its is done jointly with the resolution of the system in order to
            % keep the decomposition of the matrix.
            % SYNTAX
            %  inv_diag = Core_Utils.fastInvDiag(N,B,<mode>)
            
            if nargin < 3
                mode = 'ldl';
            end
            n_p = size(N,1);
            if strcmpi(mode,'ldl')
                % rememeber : inv(N) ==  P*iL'*inv(D)*iL*P'
                [L,D,P] = ldl(N);
                y = L\B;
                x = (D*L')\y;
                x= P*x;
                iL = inv(L);
                inv_diag = P'*sum(iL.*repmat(1./diag(D),1,n_p).*iL)';
            elseif  strcmpi(mode,'chol')
                [R] = chol(N);
                R = chol(N);
                y = R'\B;
                x=R\y;
                iR = inv(R);
                inv_diag = sum(iR.^2,2);
            end
        end
        
        function [f_nnan] = firstNoNan(mat)
            % find first no nan value of the colums
            %
            % SYNTAX
            %    [f_nnan] = Core_Utils.firstNoNan(mat)
            f_nnan = zeros(1,size(mat,2));
            for i = 1 : length(f_nnan)
                idx = find(~isnan(mat(:,i)),1,'first');
                if ~isempty(idx)
                    f_nnan(i) = mat(idx,i);
                else
                    f_nnan(i) = 0;
                end
            end
        end
        
        function [lid] = ordinal2logical(id,n_el)
            % ordinal 2 logical index
            %
            % SYNTAX
            %    [lid] = Core_Utils.ordinal2logical(id,n_el)
            lid = false(n_el,1);
            lid(id) = true;
        end
        
        
        function x = solveLDL(L,D,b)
            % solve system where normal matrix has beeen ldl decomposed
            % NOTE : A = L*D*L' (sparse matrices)
            %
            % SYNTAX
            %    x = Core_Utils.solveLDL(L,D,b)
            y = L\b;
            y = y./diag(D);
            x = L'\y;
        end
        
        function [dtm, lat, lon, georef, info] = getDTM(nwse, res)
            % Get the dtm of an area delimited by geographical coordinates nwse
            %
            % INPUT
            %   nwse(1)     North [deg: -90:90]
            %   nwse(2)     West  [deg: -180:180]
            %   nwse(3)     South [deg: -90:90]
            %   nwse(4)     East  [deg: -180:180]
            %   res         resolution ('high' / 'low') - default low
            %
            % If the DTM is not found in the DTM folder of the project
            % download it from -> use http://www.marine-geo.org/services/
            %
            % SYNTAX
            %  [dtm, lat, lon, georef, info] = Core_Utils.getDTM(nwse, res)
            
            if nargin == 1
                res = 'low';
            end
            
            dtm_path = File_Name_Processor.getFullDirPath(fullfile(Core.getInstallDir(), '../data/reference/DTM/'));
            local_dtm_path = fullfile(Core.getState.getHomeDir, 'reference' , 'DTM');
            if ~(exist(dtm_path, 'dir') == 7)
                dtm_path = local_dtm_path;
            end
            dtm_name = sprintf('dtm_N%06dW%07d_S%06dE%07d_%s.tiff', round(1e2*nwse(1)), round(1e2*nwse(2)), round(1e2*nwse(3)), round(1e2*nwse(4)), res);
            if ~exist(dtm_path, 'dir')
                try
                    mkdir(dtm_path);
                catch ex
                    Logger.getInstance.addWarning(sprintf('Could not create %s\nUsing local folder\nException: %s', dtm_path, ex.message))
                    % no write permission
                    dtm_path = fullfile('reference' , 'DTM');
                    dtm_name = sprintf('dtm_N%06dW%07d_S%06dE%07d_%s.tiff', round(1e2*nwse(1)), round(1e2*nwse(2)), round(1e2*nwse(3)), round(1e2*nwse(4)), res);
                    if ~exist(dtm_path, 'dir')
                        try
                            mkdir(dtm_path);
                        catch ex
                            Logger.getInstance.addError(sprintf('Could not create %s DTM cannot be saved locally\nException: %s\n', dtm_path, ex.message))
                        end
                    end
                end
            end
            if ~exist(fullfile(dtm_path, dtm_name), 'file')
                if exist(fullfile(local_dtm_path, dtm_name), 'file')
                    dtm_path = local_dtm_path;
                else
                    if ispc()
                        aria2c_path = '.\utility\thirdParty\aria2-extra\aria2_win\aria2c.exe';
                    elseif ismac()
                        aria2c_path = '/usr/local/bin/aria2c';
                    else % is linux
                        aria2c_path = '/usr/local/bin/aria2c';
                        if ~exist(aria2c_path, 'file')
                            aria2c_path = '/usr/bin/aria2c';
                        end
                    end
                    aria_call = sprintf('%s "%snorth=%f&west=%f&south=%f&east=%f%s%s" --dir="%s" --out="%s"', aria2c_path, 'https://www.gmrt.org/services/GridServer.php?', nwse(1), nwse(2), nwse(3), nwse(4) , '&layer=topo&format=geotiff&resolution=', res, dtm_path, dtm_name);
                    Logger.getInstance.addMarkedMessage(['Executing: "' aria_call '"']);
                    dos(aria_call)
                end
            end
            
            % Read DTM
            try
                [dtm, georef, lat, lon, info] = geotiffReader(fullfile(dtm_path, dtm_name));
            catch ex
                Logger.getInstance.addError(sprintf('Aria failed to download the required DTM file\nException: %s', ex.message));
                dtm = zeros(2,2);
                georef = [];
                lat = nwse([3 1]);
                lon = nwse([2 4]);
                info = [];
            end
        end
        
        function playAlert()
            % Warning at the end of a job! (Play a sound)
            load handel;
            sound(y(1:16000),Fs);
        end
    end
end
