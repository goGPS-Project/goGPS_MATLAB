%   CLASS Coordinates
% =========================================================================
%
% DESCRIPTION
%   Class to manage coordinates
%
% EXAMPLE
%   pos = Coordinates();
%
% CONSTRUCTOR SYNTAX
%   t = Coordinates.fromXYZ(XYZ);
%
% FOR A LIST OF CONSTANTs and METHODS use doc Coordinates
%
% COMMENTS
% The class stores arrays of time, not just a single element,
% it has been designed this way because MATLAB works best on arrays
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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


classdef Coordinates < Exportable & handle
    
    properties (Constant, GetAccess = public)
        ELL_A = GPS_SS.ELL_A;       % Ellipsoid semi-major axis
        ELL_E = GPS_SS.ELL_E;       % Ellipsoid eccentricity
        E2 = Coordinates.ELL_E^2;   % Square of the ellipsoidal eccentricity
        
        DEG2RAD = pi/180;           % Convert degree to radius
        RAD2DEG = 180/pi;           % Convert radius to degree
    end
    
    properties (SetAccess = public, GetAccess = public) % set permission have been changed from private to public (Giulio)
        xyz = [];                   % Coordinates are stored in meters in as cartesian XYZ ECEF
        precision = 0.0001;         % 3D limit [m] to check the equivalence among coordinates
        Cxx = [];
    end
        
    % =========================================================================
    %    
    % =========================================================================
    
    methods
        function this = Coordinates()
            % Constructor
            %
            % SYNTAX
            %   t = Coordinates();            
        end
                        
        function copyFrom(this, pos)
            % Copy from an object of the same type
            %
            % SYNTAX
            %   this.copyFrom(pos)
            this.xyz = pos.xyz;
            this.Cxx = pos.Cxx;
        end
        
        function copy = getCopy(this)
            % Get a copy of this
            %
            % SYNTAX
            %   copy = getCopy(this)
            copy = Coordinates();
            copy.copyFrom(this);
        end
        
        function this = append(this, pos)
            % Append a Coordinates object into the this
            %
            % SYNTAX
            %   this = append(this, pos)
            
            this.xyz = [this.xyz; pos.xyz];
            if ~isempty(this.Cxx) &&  ~isempty(pos.Cxx)
                this.Cxx = cat(this.Cxx, pos.Cxx, 3);
            end
        end
        
        function rem(this, idx)
            % Remove coordinates into the this
            %
            % SYNTAX
            %   this = rem(this, idx)
            this.xyz(idx,:) = [];
            if ~isempty(this.Cxx) 
                this.Cxx(:,:,idx) = [];
            end
        end
    end
        
    % =========================================================================
    %    GETTERS
    % =========================================================================
    
    methods
        function coo = getElement(this, id_el)
            % get a copy of the coordinates relative to a specific id (or epoch)
            %
            % UNFINISHED FUNCITON: 
            %   covariances and precisions are not extracted
            %
            % SYNTAX
            %   coo = getElement(this, id_el)
            try
                coo = Coordinates.fromXYZ(this.xyz(id_el, :));
            catch
                Core.getLogger.addWarning('Coordinates get Element out of bound');
                coo = Coordinates();
            end           
        end
                
        function [xyz, y, z] = getXYZ(this)
            % Get Coordinates as cartesian Earth Centered Earth Fixed Coordinates
            %
            % OUTPUT
            %   xyz      = coordinates     [m]
            %
            % SYNTAX 
            %   [xyz] = this.getXYZ()
            %   [x, y, z] = this.getXYZ()
           
            if nargout == 3
                xyz = this.xyz(:, 1);
                y =   this.xyz(:, 2);
                z =   this.xyz(:, 3);
            else
                xyz = this.xyz;
            end
        end
        
        function [east, north, utm_zone] = getUTM(this)
            % Get Coordinates as UTM coordinates
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX 
            %   [east, north, utm_zone] = this.getUTM();
            
            [lat, lon] = this.getGeodetic();
            [east, north, utm_zone] = this.geod2plan(lat, lon);
        end
                       
        function [east, north, up, utm_zone] = getENU(this, theta)
            % Get Coordinates as UTM ENUs coordinates
            %
            % INPUT
            %   theta   planar rotation angle [degree -360:360]
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   up       = up                [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX 
            %   [east, north, up, utm_zone] = this.getENU(<theta>);
            %   [enu, utm_zone] = this.getENU(<theta>);
            
            [lat, lon, up] = this.getGeodetic();
            [east, north, utm_zone] = this.geod2plan(lat, lon);
            
            if nargin == 2 && ~isempty(theta)
                tmp = [east(:) north(:)] * [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
                east = tmp(:,1);
                north = tmp(:,2)';
            end            
            if nargout < 3
                east = [east(:) north(:) up(:)];
                north = utm_zone;
            end            
        end        
        
        function [lat, lon, h_ellips, h_ortho] = getGeodetic(this)
            % Get Coordinates as Geodetic coordinates
            %
            % OUTPUT
            %   lat      = latitude                      [rad]
            %   lon      = longitude                     [rad]
            %   h_ellips = ellipsoidal height            [m]
            %   h_ortho  = orthometric height            [m]
            %
            % SYNTAX 
            %   [lat, lon, h_ellips, h_ortho] = this.getGeodetic()
           
            if nargout > 2
                [lat, lon, h_ellips] = Coordinates.cart2geod(this.xyz);
                if nargout == 4
                    h_ortho = h_ellips - this.getOrthometricCorrFromLatLon(lat, lon);
                end
            else
                [lat, lon] = Coordinates.cart2geod(this.xyz);
            end
        end
        
        function [lat_geoc, lon, h] = getGeocentric(this)
            % Get Coordinates as Geodetic coordinates
            %
            % OUTPUT
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %
            % SYNTAX 
            %   [lat, lon, h] = this.getGeocentric();
           
            [lat_geoc, lon, h] = Coordinates.cart2geoc(this.xyz);
        end
        
        function ondu = getOrthometricCorrection(this)
            % Get Orthometric correction from the geoid loaded in Core
            %
            % OUTPUT
            %   ondu     = geoid ondulation [m]
            %
            % SYNTAX 
            %   ondu = getOrthometricCorrection(this)
            
            [lat, lon] = this.getGeodetic();
            ondu = this.getOrthometricCorrFromLatLon(lat, lon);
        end
        
        function loc = getLocal(this, ref_pos)
            % Get Coordinates as Local coordinates with respect to ref_pos
            %
            % OUTPUT
            %   loc     = xyz local coordinates
            %
            % SYNTAX 
            %   loc = this.getLocal(ref_pos)
           
            xyz_ref = ref_pos.getXYZ;
            baseline = ref_pos.getXYZ - this.getXYZ;
            loc = Coordinates.cart2loca(xyz_ref, baseline);
        end
        
        function status = isEmpty(this)
            % Return the status of emptyness of the object
            %
            % SYNTAX
            %   status this.isEmpty();
            status = isempty(this.xyz);
        end  
        
        function dist = ellDistanceTo(this,coo)
            % return the distance on the ellipsoid between the object
            % coordinates and the coordinate coo using vincenty's formula
            % (https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)
            %
            % SYNTAX
            %     dist = this.distanceTo(coo)
            
            % reduced latitude
            a = this.ELL_A;
            b = sqrt((1 - this.E2)*a^2) ;
            f = (a - b)/a;
            [lat1, lon1] = this.getGeodetic();
            [lat2, lon2] = coo.getGeodetic();
            U1 = atan((1-f)*tan(lat1));
            U2 = atan((1-f)*tan(lat2));
            L = lon2 - lon1;
            lambda_old = 10000;
            lambda= L;
            
            while abs(lambda - lambda_old) > 1e-12
                sin_sigma = sqrt((cos(U2) * sin(lambda))^2 + (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda))^2);
                cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
                sigma = atan2(sin_sigma,cos_sigma);
                sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin_sigma;
                cos2_alpha = 1 - sin_alpha^2;
                cos_2sigmam = cos_sigma - 2*sin(U1)*sin(U2)/cos2_alpha;
                C = f/16*cos2_alpha*(4 + f*(4 -3*cos2_alpha));
                lambda_old = lambda;
                lambda = L + (1 - C)*f*sin_alpha*(sigma + C*sin_sigma*(cos_2sigmam + C*cos_sigma*(-1 + 2 * cos_2sigmam^2)));
            end
            u2 = cos2_alpha*(a^2 - b^2)/b^2;
            A = 1 + u2/16384*(4096 + u2*(-768 + u2*(320-175*u2)));
            B = u2/1024*(256 +u2*(-128 +u2*(74-47*u2)));
            Dsigma = B*sin_sigma*(cos_2sigmam + 1/4*B*(cos_sigma*(-1 + 2*cos_2sigmam^2) -B/6*cos_2sigmam*(-3 + 4 * sin_sigma^2)*(- 3 + 4 * cos_2sigmam^2)));
            dist = b*A*(sigma -Dsigma);
            
        end
    end
    
    % =========================================================================
    %    SETTERS
    % =========================================================================
    
    methods
        function setPosXYZ(this, xyz, y, z)
            % Set the Coordinates
            %
            % SYNTAX
            %   this.setPosXYZ(xyz)
            %   this.setPosXYZ(x, y, z)
            
            if nargin == 4
                xyz = [xyz(:) y(:) z(:)];
            end
            this.xyz = xyz;
        end
        
        function empty(this)
            % Empty the object
            %
            % SYNTAX
            %   this.empty();
            this.xyz = [];
        end
    end
    
    % =========================================================================
    %    STATIC CONSTRUCTOR
    % =========================================================================
    methods (Access = 'public', Static)
        function this = fromXYZ(xyz, y, z)
            % Set the Coordinates from XYZ coordinates
            %
            % SYNTAX
            %   this = Coordinates.fromXYZ(xyz)
            %   this = Coordinates.fromXYZ(x, y, z)
            
            this = Coordinates;
            if nargin > 2
                xyz = [x(:) y(:) z(:)];
            end
            this.setPosXYZ(xyz);
        end
        
        function this = fromStringXYZ(xyz_string)
            % Set the Coordinates from XYZ coordinates (String)
            %
            % SYNTAX
            %   this = Coordinates.fromStringXYZ(xyz)
            
            this = Coordinates;
            xyz = sscanf(xyz_string, '%f%f%f')';
            if numel(xyz) ~= 3
                xyz = [0 0 0];
            end            
            this.setPosXYZ(xyz);
        end
        
        function this = fromGeodetic(lat, lon, h_ellips, h_ortho)
            % Set the Coordinates from Geodetic coordinates
            %
            % INPUT
            %   lat, lon [radians]
            %
            % SYNTAX
            %   this = Coordinates.fromGeodetic(phi, lam, h_ellips);
            %   this = Coordinates.fromGeodetic(phi, lam, [], h_ortho);
            
            this = Coordinates;
            if nargin == 4
                h_ellips = h_ortho + this.getOrthometricCorrFromLatLon(lat, lon);
            end
            
            N = this.ELL_A ./ sqrt(1 - this.E2 * sin(lat).^2);
            
            x = (N + h_ellips) .* cos(lon) .* cos(lat);
            y = (N + h_ellips) .* sin(lon) .* cos(lat);
            z = (N * (1 - this.E2) + h_ellips) .* sin(lat);

            this.setPosXYZ(x, y, z);
        end                    
    end
    
    % =========================================================================
    %    SHOWs
    % =========================================================================
    
    methods (Access = 'public')
        function fh = showCoordinatesENU(coo_list, one_plot, time)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   this.showCoordinatesENU(coo_list, <one_plot>, time);
            if nargin < 2 || isempty(one_plot)
                one_plot = false;
            end
            
            set_time = false;
            log = Core.getLogger();
            for i = 1 : numel(coo_list)
                pos = coo_list(i);
                if ~pos.isEmpty
                    enu = pos.getENU;
                    if size(enu, 1) > 1
                        log.addMessage('Plotting coordinates');
                        
                        fh = figure; fh.Name = sprintf('%03d: PosENU', fh.Number); fh.NumberTitle = 'off';
                        color_order = handle(gca).ColorOrder;
                        
                        enu0 = median(enu, 1, 'omitnan');
                        
                        if nargin == 3 && ~isempty(time)
                            if isa(time, 'GPS_Time')
                                set_time = true;
                                if numel(time) == numel(coo_list)
                                    t = time(i).getMatlabTime;
                                else
                                    t = time.getMatlabTime;
                                end
                            else
                                t = time;
                            end
                        else                            
                            t = 1 : size(enu, 1);
                        end
                        
                        if ~one_plot, subplot(3,1,1); end
                        plot(t, (1e2 * (enu(:,1) - enu0(1))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        title(sprintf('ENU Coordinate\nstd = %.2f mm, %.2f mm, %.2fmm\\fontsize{5} \n', ...
                            std( 1e3 * (enu(:,1) - enu0(1)) ,'omitnan'), ...
                            std( 1e3 * (enu(:,2) - enu0(2)) ,'omitnan'), ...
                            std( 1e3 * (enu(:,3) - enu0(3)) ,'omitnan')))
                        ax(3) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('East [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if ~one_plot, subplot(3,1,2); end
                        plot(t, (1e2 * (enu(:,2) - enu0(2))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('North [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if ~one_plot, subplot(3,1,3); end
                        plot(t, (1e2 * (enu(:,3) - enu0(3))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('Up [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if one_plot
                            h = ylabel('ENU [cm]'); h.FontWeight = 'bold';
                        else                           
                            linkaxes(ax, 'x');
                        end
                        grid on;
                        
                    else
                        log.addMessage('Plotting a single point static Coordinates is not yet supported');
                    end
                end
            end
        end
        
        function fh = showCoordinatesXYZ(coo_list, one_plot, time)
            % Plot X Y Z coordinates
            %
            % SYNTAX
            %   this.showCoordinatesXYZ(coo_list, <one_plot>, <time>);
            if nargin == 1
                one_plot = false;
            end
                 
            set_time = false;
            log = Core.getLogger();
            for i = 1 : numel(coo_list)
                coo = coo_list(i);
                if ~coo.isEmpty
                    xyz = coo.getXYZ;
                    
                    if size(xyz, 1) > 1
                        log.addMessage('Plotting coordinates');
                        
                        fh = figure; fh.Name = sprintf('%03d: PosXYZ', fh.Number); fh.NumberTitle = 'off';
                        color_order = handle(gca).ColorOrder;
                        
                        xyz0 = median(xyz, 1, 'omitnan');
                        
                        if nargin == 3 && ~isempty(time)
                            if isa(time, 'GPS_Time')
                                set_time = true;
                                if numel(time) == numel(coo_list)
                                    t = time(i).getMatlabTime;
                                else
                                    t = time.getMatlabTime;
                                end
                            else
                                t = time;
                            end
                        else
                            t = 1 : size(enu, 1);
                        end
                        x = 1e2 * bsxfun(@minus, zero2nan(xyz(:,1)), xyz0(1));
                        y = 1e2 * bsxfun(@minus, zero2nan(xyz(:,2)), xyz0(2));
                        z = 1e2 * bsxfun(@minus, zero2nan(xyz(:,3)), xyz0(3));
                        
                        if ~one_plot, subplot(3,1,1); end
                        plot(t, x, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                        ax(3) = gca(); xlim([t(1) t(end)]);
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('X [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if ~one_plot, subplot(3,1,2); end
                        plot(t, y, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca(); xlim([t(1) t(end)]);
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('Y [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if ~one_plot, subplot(3,1,3); end
                        plot(t, z, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca(); xlim([t(1) t(end)]);
                        if set_time, setTimeTicks(4,'dd/mm/yyyy HH:MM'); end
                        h = ylabel('Z [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if one_plot
                            h = ylabel('XYZ [cm]'); h.FontWeight = 'bold';
                        end
                        linkaxes(ax, 'x');
                    else
                        log.addMessage('Plotting a single point static Coordinates is not yet supported');
                    end
                end
            end
        end
    end
    
    % =========================================================================
    %    OPERATIONS
    % =========================================================================
    
    methods (Access = 'public')                        
        function res = eq(gt_1, gt_2)
            %%% DESCRIPTION: check if two coordinates are equal
            d = sqrt(sum((gt_1.xyz - gt_2.xyz).^2, 2));
            res = d < gt_1.precision;
        end        
    end
        
    
    methods (Access = 'public', Static)
        function ondu = getOrthometricCorrFromLatLon(lat, lon)
            % Get Orthometric correction from the geoid loaded in Core
            %
            % INPUT
            %   lat, lon    [radians]
            % SYNTAX:
            %   ondu = getOrthometricCorrFromLatLon(lat, lon);
            
            geoid = Core.getRefGeoid();
            
            ondu = zeros(numel(lon), 1);
            for i = 1 : numel(lon)
                ondu(i) = grid_bilin_interp(lon(i) * Coordinates.RAD2DEG, lat(i) * Coordinates.RAD2DEG, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
            end            
        end
        
        function [lat, lon, h, lat_geoc] = cart2geod(xyz, y, z)
            % Get Geodetic coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(xyz)
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(x, y, z)
            
            if nargin == 1
                z = xyz(:, 3);
                y = xyz(:, 2);
                xyz = xyz(:, 1);
            end
            
            a = Coordinates.ELL_A;
            e = Coordinates.ELL_E;
            
            % radius computation
            r = sqrt(xyz.^2 + y.^2 + z.^2);
            
            % longitude
            lon = atan2(y,xyz);
            
            % geocentric latitude
            lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));
            
            % Coordinates transformation
            psi = atan(tan(lat_geoc)/sqrt(1-e^2));
            
            lat = atan((r.*sin(lat_geoc) + e^2*a/sqrt(1-e^2) * (sin(psi)).^3) ./ ...
                (r.*cos(lat_geoc) - e^2*a * (cos(psi)).^3));
            
            if nargout > 2
                N = a ./ sqrt(1 - e^2 * sin(lat).^2);
                
                % height
                h = r .* cos(lat_geoc)./cos(lat) - N;
            end
        end
        
        function [lat_geoc, lon, r] = cart2geoc(xyz, y, z)
            % Get Geocentric sphericals coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(xyz)
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(x, y, z)
            
            if nargin == 1
                z = xyz(:, 3);
                y = xyz(:, 2);
                xyz = xyz(:, 1);
            end
            
            a = Coordinates.ELL_A;
            e = Coordinates.ELL_E;
            
            % radius computation
            r = sqrt(xyz.^2 + y.^2 + z.^2);
            
            % longitude
            lon = atan2(y,xyz);
            
            % geocentric latitude
            lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));            
        end
        
        function [east, north, utm_zone] = geod2plan(lat, lon)
            % Conversion from geodetic coordinates to planimetric coordinates (UTM WGS84).
            %
            % SYNTAX:
            %   [east, north, utm_zone] = geod2plan(lat, lon);
            %
            % INPUT:
            %   lat = latitude [rad]
            %   lon = longitude [rad]
            %
            % OUTPUT:
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %       
            
            %number of input points
            n = size(lat);
            n = n(1,1);
            
            if (n > 0)
                %pre-allocation
                M = zeros(n,1);
                north_south = zeros(n,1);
                utm_zone(n,:) = '60 X';
                north = zeros(n,1);
                east = zeros(n,1);
                
                % conversion algorithm
                % UTM parameters
                EsMc = 500000; % false East
                
                % WGS84 ellipsoid parameters
                % semi-major equatorial axis [m]
                ell_a = Coordinates.ELL_A;
                
                % squared eccentricity (a^2-b^2)/a^2
                e2 = Coordinates.E2;
                
                % contraction factor
                contr = 0.9996;
                
                ell_a = ell_a * contr;
                ecc4 = e2 * e2;
                ecc6 = ecc4 * e2;
                ecc8 = ecc6 * e2;
                
                k0 = ell_a * (e2 / 4 + ecc4 * 3 / 64 + ecc6 * 5 / 256 + ecc8 * 175 / 16384);
                k = ell_a - k0;
                k1 = ell_a * (ecc4 * 13 / 96 + ecc6 * 59 / 384 + ecc8 * 1307 / 8192);
                k2 = ell_a * (ecc6 * 61 / 484 + ecc8 * 609 / 2048);
                k3 = ell_a * (ecc8 * 49561 / 322560);
                c1 = (e2 * 5 - ecc4) / 6;
                c2 = (ecc4 * 104 - ecc6 * 45) / 120;
                c3 = ecc6 * 1237 / 1260;
                
                % Sines, cosines and latitude powers
                Elix = [lat lon];
                latsessadec = lat ./ pi .* 180;
                lonsessadec = lon ./ pi .* 180;
                
                fiSin(:,1) = sin(Elix(:,1));
                fiCos(:,1) = cos(Elix(:,1));
                fiSin2(:,1) = fiSin .* fiSin;
                fiSin4(:,1) = fiSin2 .* fiSin2;
                fiSin6(:,1) = fiSin4 .* fiSin2;
                
                % UTM zone finding
                for i = 1 : n
                    
                    M(i, 1) = fix((180 + lonsessadec(i, 1)) / 6) + 1;
                    
                    if latsessadec(i, 1) >= 0
                        north_south(i, 1) = 1; %1 north, 0 south
                    else
                        north_south(i, 1) = 0;
                    end
                    
                    if     (latsessadec(i, 1) < -72), letter = 'C';
                    elseif (latsessadec(i, 1) < -64), letter = 'D';
                    elseif (latsessadec(i, 1) < -56), letter = 'E';
                    elseif (latsessadec(i, 1) < -48), letter = 'F';
                    elseif (latsessadec(i, 1) < -40), letter = 'G';
                    elseif (latsessadec(i, 1) < -32), letter = 'H';
                    elseif (latsessadec(i, 1) < -24), letter = 'J';
                    elseif (latsessadec(i, 1) < -16), letter = 'K';
                    elseif (latsessadec(i, 1) <  -8), letter = 'L';
                    elseif (latsessadec(i, 1) <   0), letter = 'M';
                    elseif (latsessadec(i, 1) <   8), letter = 'N';
                    elseif (latsessadec(i, 1) <  16), letter = 'P';
                    elseif (latsessadec(i, 1) <  24), letter = 'Q';
                    elseif (latsessadec(i, 1) <  32), letter = 'R';
                    elseif (latsessadec(i, 1) <  40), letter = 'S';
                    elseif (latsessadec(i, 1) <  48), letter = 'T';
                    elseif (latsessadec(i, 1) <  56), letter = 'U';
                    elseif (latsessadec(i, 1) <  64), letter = 'V';
                    elseif (latsessadec(i, 1) <  72), letter = 'W';
                    else,                             letter = 'X';
                    end
                    
                    utm_zone(i,:) = sprintf('%02d %c', nan2zero(M(i, 1)), letter);
                end
                
                for i = 1 : n
                    % Longitude of the central meridian
                    LonMeridianoCentrale = -177 + 6 * (M(i, 1) - 1);
                    
                    % Distance of the point from the central meridian
                    % la_sd --> distance in decimal degrees
                    % la    --> distance in radians
                    la_sd = lonsessadec(i, 1) - LonMeridianoCentrale;
                    la = la_sd / 180 * pi;
                    
                    if la == 0
                        laSin = 0;
                        laCos = 1;
                        laCot = 0;
                    else
                        laSin = sin(la);
                        laCos = cos(la);
                        laCot = laCos / laSin ;
                    end
                    
                    % longitude with respect to central meridian
                    laCot2 = laCot * laCot;
                    
                    % psi
                    psi = Elix(i, 1) - e2 * fiSin(i, 1) * fiCos(i, 1) *(1 + c1 * fiSin2(i, 1) + c2 * fiSin4(i, 1) + c3 * fiSin6(i, 1));
                    psiSin = sin(psi);
                    psiCos = cos(psi);
                    psiTan = psiSin / psiCos ;
                    psiSin2 = psiSin * psiSin ;
                    
                    % omega
                    ome = atan(psiTan / laCos);
                    
                    % sigma
                    if laSin ~= 0
                        sigSin =laSin * psiCos;
                        sig = asin(sigSin);
                    else
                        sigSin = 0;
                        sig = 0;
                    end
                    
                    sigSin2 = sigSin * sigSin;
                    sigSin4 = sigSin2 * sigSin2;
                    sigCos2 = 1 - sigSin2;
                    sigCos4 = sigCos2 * sigCos2;
                    
                    % chi
                    chi = sig / 2 + pi / 4;
                    chiTan = tan(chi);
                    chiLog = log(chiTan);
                    
                    % constants
                    aa = psiSin * psiCos * laCos * (1 + sigSin2) / sigCos4;
                    bb = sigSin * (sigCos2 - 2 * psiSin2) / sigCos4;
                    
                    if laCot ~= 0
                        a1 = (psiSin2 - sigSin4 * laCot2)/ sigCos4;
                        b1 = 2 * sigSin2 * psiSin * laCot / sigCos4;
                    else
                        a1 = psiSin2 / sigCos4 ;
                        b1 = 0;
                    end
                    a2 = a1 * a1 - b1 * b1 ;
                    b2 = 2 * a1 * b1 ;
                    a3 = a1 * a2 - b1 * b2 ;
                    b3 = a1 * b2 + b1 * a2 ;
                    rr = k0 - a1 * k1 + a2 * k2 - a3 * k3;
                    tt = b1 * k1 - b2 * k2 + b3 * k3;
                    
                    % X/Y coordinates
                    xx = k * ome + aa * rr + bb * tt;
                    yy = k * chiLog + bb * rr - aa * tt;
                    
                    % North and East
                    if north_south(i, 1) == 1
                        NoEq = 0;
                    else
                        NoEq = 10000000;
                    end
                    
                    north(i, 1) = NoEq + xx;
                    east(i, 1) = EsMc + yy;
                end                
            else
                utm_zone = [];
                north = [];
                east = [];
            end
        end
        
        function [loc,  rot_mat] = cart2local(xyz_ref, xyz_baseline)
            % cart2loca: from geocentric cartesian baselines (DX) to local coordinates in X0
            %
            % SYNTAX
            %   [loc,  rot_mat] = Coordinates.cart2local(xyz_ref, baseline)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
            loc = (rot_mat * xyz_baseline')';
        end
        
        function [loc,  rot_mat] = cart2loca(xyz_ref, xyz_baseline)
            [loc,  rot_mat] = Coordinates.cart2local(xyz_ref, xyz_baseline); % maybe the function was used with the wrong name
        end
        
        function [xyz_baseline, rot_mat] = loca2cart(xyz_ref, loc)
            % loca2cart: from local coordinates in X0 to geocentric cartesian baselines (DX)
            %
            % SYNTAX
            %   [baseline, rot_mat] = Coordinates.loca2cart(xyz_ref, loca)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]';
            xyz_baseline = (rot_mat * loc)';
        end                
    end
    
    % =========================================================================
    %    TESTS
    % =========================================================================
    
    methods (Static, Access = 'public')
        function test()
            % Testing function, tests some basic transformations
            %
            % SYNTAX
            %   test()
            log = Core.getLogger();
            
            log.addMessage('Testing Class Coordinates');
            tic;
            pos_diff = 0;

            if pos0 == pos1;
                log.addStatusOk('Passed');
            else
                log.addWarning(sprintf('Difference greater than 0.2 ms: %e',t_diff));
            end
            toc
        end
    end
    
end
