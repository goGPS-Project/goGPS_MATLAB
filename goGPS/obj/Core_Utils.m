
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
%    |___/                    v 1.0 beta 4 ION
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
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 2));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
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
               
        
        function [z, l, m] = getAllZernike(l_max, m_max, az, el)
            % Generate all the Zernike parameters combinations
            %
            % SINTAX
            %   z = getAllZernike(l_max, az, el)
            %   z = getAllZernike(l_max, m_max, az, el)
            if nargin == 3
                el = az;
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
            
            z = Core_Utils.getZernike(l, m, az, el);
        end
            
        function z = getZernike(l, m, az, el)
            % Get Zernike values for the polynomials
            %
            % SINTAX
            %   z = getZernike(l, m, az, el)
            
            r = 1 - (2 * el(:) / pi);
            %r = cos(el(:));
            theta = az(:);
            z = zernfun(l, m, r, theta);
        end
        
        function [z_interp, l, m] = zAnalisysAll(l_max, m_max, az, el, z_par)
            % Get Zernike interpolation given the coefficients
            % of their polynomials
            %
            % SINTAX
            %   [z_interp] = zAnalisysAll(l_max, m_max, az, el, data)
            z_interp = nan(size(az));
            
            id_ok = ~isnan(az);
            [A, l, m] = Core_Utils.getAllZernike(l_max, m_max, l_max, m_max, az, el);
            z_interp(id_ok) = A * z_par;
        end
        
        function [z_interp, l, m] = zAnalisys(l, m, az, el, z_par)
            % Get Zernike interpolation given the coefficients
            % of their polynomials
            %
            % SINTAX
            %   [z_interp] = zAnalisys(l, m, az, el, data)
            z_interp = nan(size(az));
            
            id_ok = ~isnan(az);
            A = Core_Utils.getZernike(l, m, l_max, m_max, az, el);
            z_interp(id_ok) = A * z_par;
        end
        
        function [z_par, l, m, A] = zSinthesysAll(l_max, m_max, az, el, data, max_reg)
            % Get Zernike polynomials parameters 
            %
            % SINTAX
            %   [z_par, l, m, A] = zSinthesysAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:));
            [A, l, m] = Core_Utils.getAllZernike(l_max, m_max, az(id_ok), el(id_ok));
            if nargin == 6 && ~isempty(max_reg) && max_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 0;
            end
            z_par = (A'*A + reg_fun .* diag(ones(size(A, 2), 1))) \ A' * data(id_ok);
        end
        
        function [z_par, l, m, A] = zSinthesys(l, m, az, el, data, max_reg)
            % Get Zernike polynomials parameters 
            %
            % SINTAX
            %   [z_par, l, m, A] = zSinthesysAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:));
            A = Core_Utils.getZernike(l, m, az(id_ok), el(id_ok));
            if nargin == 6 && ~isempty(flag_reg) && flag_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 0;
            end
            z_par = (A' * A + reg_fun .* diag(ones(size(A, 2), 1))) \ A' * data(id_ok);
        end
                
        function [filtered_data, z_par, l, m,  A] = zFilter(l_max, m_max, az, el, data, max_reg)           
            % Get Zernike polynomials parameters 
            %
            % SINTAX
            %   [z_par, A] = zFilter(az, el, data, <max_reg = 1>)
            
            if nargin == 6
                [z_par, l, m, A] = Core_Utils.zSinthesysAll(l_max, m_max, az, el, data, max_reg);
            else
                [z_par, l, m, A] = Core_Utils.zSinthesysAll(l_max, m_max, az, el, data);
            end
            filtered_data = A * z_par;
        end

        function fh = showZerniche(l, m, z_par, el_min)
            % Get Zernike polynomials parameters 
            %
            % SINTAX
            %   fh = showZerniche(l, m, z_par)
            
            % [x, y] = pol2cart(theta, r_synt);
            if nargin == 4
                r_max = 1 - (2 * el_min / pi);
                [theta, r_synt] = meshgrid(linspace(0, 2*pi, 361), linspace(0, r_max, 101));
            else
                [theta, r_synt] = meshgrid(linspace(0, 2*pi, 361), linspace(0, 1, 101));
            end
                
            z = nan(size(theta));
            z(:) = zernfun(l, m, r_synt(:), theta(:)) * z_par;
            
            fh = figure(102);
            title('Zerniche expansion')
            %polarplot3d(z, 'PlotType','surfn');
            polarplot3d(z,'PlotType','surfn','PolarGrid',{4 24},'TickSpacing',8,...
                   'AngularRange',[0 360]*pi/180,'RadialRange', [0 1],...
                   'RadLabels',4,'RadLabelLocation',{180 'max'},'RadLabelColor','black', 'AxisLocation', 'mean');
            ar = get(gca,'DataAspectRatio');
            set(gca, 'DataAspectRatio', [1 1 ar(3)]);
            
            colormap(jet);
            material([ 0.4 0.9 0.55])
            l1=light('position',[200 -300 400], 'color', [0.6 0.6 0.6]);
            l2=light('position',[-600 600 900], 'color', [0.6 0.6 0.6]);
            
            Core_UI.beautifyFig(fh, 'dark');            
        end

        %--------------------------------------------------------------------------
        % OTHER FUNCTIONS
        %--------------------------------------------------------------------------
        
        function printEx(ex)
            % Print exception to screen
            %
            % SYNTAX:
            %  Core_Utils.printEx(ex)
            
            fprintf('\n---------------------------------------------------------------------\n MESSAGE: %s\n---------------------------------------------------------------------\n\n', ex.message);
            for i=1:numel(ex.stack)
                fprintf('  file: "%s"\n  line: %d\n  fun: %s\n\n', ex.stack(i).file, ex.stack(i).line, ex.stack(i).name);
            end
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
            export_fig(fh, out_path, '-transparent', '-r150');
            fh.WindowStyle = ws_bk;
            fh.Color = col;
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
            %   str3 = Core_Utils.num2Code3ch(num)
            str2 = char(zeros(numel(num), 2));
            str2(:,1) = char(floor(num / 2^8));
            num = num - str2(:,2) * 2^8;
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
            %   f_status_lst        bool array of files to download                [bool]
            %   date_list           GPS_Time of days of interest                   [GPS_Time]
            %   out_dir             path to out folder                             [char]
            %
            % SYNTAX
            %   f_status_lst = Core_Utils.aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst)
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
                fid = fopen(file_name, 'w');
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
                                fid = fopen(file_name, 'w');
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
                if strfind(txt,' 200 OK')
                    status = true;
                else
                    [resp, txt] = system([rem_check_cmd filename '.gz']);
                    if strfind(txt,' 200 OK')
                        ext = '.gz';
                        status = true;
                    else
                        [resp, txt] = system([rem_check_cmd filename '.Z']);
                        if strfind(txt,' 200 OK')
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
            station_list = unique(file_list)
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
            % Merge
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
            
            dtm_path = fullfile(Core.getState.getHomeDir, 'reference' , 'DTM');
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
    end
end
