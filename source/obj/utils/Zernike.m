%   CLASS Core
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage a useful parameters of goGPS
%   This singleton class collects multiple objects containing various
%   parameters
%
% EXAMPLE
%   core = Core();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core
            
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro...
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

classdef Zernike < handle
    properties 
        mode = 0;  % -1 no discratization FileExchange
                   %  0 no discratization Recursive
                   %  1 R discratized Recursive
                   %  2 R and Az discratized and cached Recursive
        
        cut_off = 0; % cut_off level
        
        el_mf = 0; %  0 cos(el)
                   %  1 cos(el).^2
                   %  2 sin(pi/2 * cos(el).^2)
                   %  3 sin(pi/2 * cos(el))

       z_cache; % cache for mode 2
    end
    
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function this = Zernike(mode)
            % Initialisation of the variables
            if (nargin == 1) && ~isempty(mode)
                this.mode = mode;
            end
        end
    end
    
    % Getters and setters
    methods (Static)
        function this = getInstance()
            % Concrete implementation.  See Singleton superclass.
            persistent unique_instance_zernike__
            if isempty(unique_instance_zernike__)
                this = Zernike();
                unique_instance_zernike__ = this;
            else
                this = unique_instance_zernike__;
            end
        end
        
        function z_cache = getCache()
            % Get for the persistent object the current cache
            %
            % SYNTAX:
            %   z_cache = Zernike.getCache();
            
            zer = Zernike.getInstance;
            z_cache = zer.z_cache;
        end
        
        function setCache(z_cache)
            % Set for the persistent object the current cache
            %
            % SYNTAX:
            %   Zernike.setCache(z_cache);
            
            zer = Zernike.getInstance;
            zer.z_cache = z_cache;
        end
        
        function mode = getMode()
            % Get for the persistent object the discretization level
            %
            % SYNTAX:
            %   z_cache = Zernike.getMode();
            
            zer = Zernike.getInstance;
            mode = zer.mode;
        end
        
        function setMode(mode)
            % Set for the persistent object the discretization level
            %
            % SYNTAX:
            %   Zernike.setMode(mode);
            
            zer = Zernike.getInstance;
            zer.mode = mode;
        end
        
        function el_mf = getModeMF()
            % Get for the persistent object the mapping function type
            %
            % SYNTAX:
            %   z_cache = Zernike.getModeMF();
            
            zer = Zernike.getInstance;
            el_mf = zer.el_mf;
        end
        
        function setModeMF(el_mf)
            % Set for the persistent object the discretization level
            %
            % SYNTAX:
            %   Zernike.setModeMF(mode);
            
            zer = Zernike.getInstance;
            zer.el_mf = el_mf;
        end
        
        function cut_off = getCutOff()
            % Get for the persistent object the mapping function type
            %
            % SYNTAX:
            %   cut_off = Zernike.getCutOff();
            
            zer = Zernike.getInstance;
            cut_off = zer.cut_off;
        end
        
        function setCutOff(cut_off)
            % Set for the persistent object the discretization level
            %
            % SYNTAX:
            %   Zernike.setCutOff(cut_off);
            
            zer = Zernike.getInstance;
            zer.cut_off = cut_off;
        end
        
        function el2radius = getElFun(type, cut_off)
            if nargin == 2 && ~isempty(cut_off)
                % First approach:
                co = cut_off;
            else
                % Second approach
                co = Zernike.getCutOff();
            end
            
            if co == 0
                % First approach:
                remCutOff = @(el) (el - co) / (pi/2 - co) * (pi / 2);
            else
                % Second approach
                remCutOff = @(el) el;
            end
            
            if nargin < 1 || isempty(type)
                type = Zernike.getModeMF;
            end
            
            % mapping functions (Default 1)
            switch type
                case -1, el2radius = @(el) (pi/2 - remCutOff(el)) / (pi/2);
                case 0, el2radius = @(el) cos(remCutOff(el));
                case 1, el2radius = @(el) cos(remCutOff(el)).^2;
                case 2, el2radius = @(el) sin(pi/2*cos(remCutOff(el)).^2);
                case 3, el2radius = @(el) sin(pi/2*cos(remCutOff(el)));
                case 4, omf = @(el) 1 - sin(el);
                    el2radius = @(el) omf(pi/2 * (1-cos(pi/2*(1-cos(el)))));
            end
        end
        
    end
    
    % Methods exposed utilities
    methods (Static)
        function [z_interp, l, m] = synthesisAll(l_max, m_max, az, r, z_par)
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
            %   [z_interp] = synthesisAll(l_max, m_max, az, r, z_par)
            z_interp = nan(size(az));
            
            id_ok = ~isnan(az);
            [A, l, m] = Zernike.getAll(l_max, m_max, az, r);
            z_interp(id_ok) = A * z_par;
        end
        
        function [z_interp, l, m] = synthesis(l, m, az, r, z_par)
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
            %   [z_interp] = Zernike.synthesis(l, m, az, el, z_par)
            z_interp = nan(size(az));
            
            idx = find(~isnan(az));
            
            max_ram = 1 * 1024^3; % 1GB
            nobs_max = (max_ram / (numel(z_par) * 8));
            lim = [unique(floor([1 : nobs_max : numel(idx)-1]))' unique(floor([(nobs_max + 1) : nobs_max : numel(idx) numel(idx)]))'];
            
            if nobs_max < numel(idx)
                % If the synthesis is large perform it in blocks
                w_bar = Go_Wait_Bar.getInstance(size(lim, 1), 'Zernike Synthesis');
                w_bar.createNewBar
                for i = 1 : size(lim, 1)                    
                    A = Zernike.get(l, m, az(idx(lim(i,1):lim(i,2))), r(idx(lim(i,1):lim(i,2))));
                    z_interp(idx(lim(i,1):lim(i,2))) = A * z_par;
                    w_bar.goTime();
                end
                w_bar.close();                                
            else
                % If the synthesis is small perform it in a unique block
                A = Zernike.get(l, m, az(:), r(:));
                z_interp(idx) = A * z_par;
            end
        end
        
        function [z_grid, az_grid, el_grid] = synthesisGrid(l, m, z_par, grid_step, el_fun_type)
            % Get Zernike interpolation given the coefficients
            % of their polynomials
            %
            % INPUT
            %   l           list of degrees
            %   m           list of orders
            %   z_par       Zernike coefficients
            %   grid_step   step of the grid [deg]
            %   el_fun_type elevation 2 radius mapping function (see Zernike.getElFun)
            %
            % SINTAX
            %   [z_interp] = Zernike.synthesisGrid(l, m, az, z_par, grid_step, 1)
            az_grid = ((-180 + (grid_step(1) / 2)) : grid_step(1) : (180 - grid_step(1) / 2)) .* (pi/180);
            el_grid = flipud(((grid_step(end) / 2) : grid_step(end) : 90 - (grid_step(end) / 2))' .* (pi/180));
               
            [az_mgrid, el_mgrid] = meshgrid(az_grid, el_grid);

            if nargin < 5 || isempty(el_fun_type)
                el_fun_type = Zernike.getModeMF();
            end
            el2radius = Zernike.getElFun(el_fun_type);
            z_grid = zeros(size(az_mgrid));
            id_ok = el_mgrid >= (Zernike.getCutOff() / 180 * pi);
            z_grid(id_ok) = Zernike.synthesis(l, m, az_mgrid(id_ok), el2radius(el_mgrid(id_ok)), z_par);            
        end
                        
        function [filtered_data, z_par, l, m,  A] = filter(l_max, m_max, az, r, data, max_reg)
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
            %   [z_par, A] = filter(az, el, data, <max_reg = 1>)
            
            if nargin == 6
                [z_par, l, m, A] = Zernike.analysisAll(l_max, m_max, az, r, data, max_reg);
            else
                [z_par, l, m, A] = Zernike.analysisAll(l_max, m_max, az, r, data);
            end
            filtered_data = A * z_par;
        end
        
        
        function [z_par, l, m, A] = analysisAll(l_max, m_max, az, r, data, max_reg)
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
            %   [z_par, l, m, A] = analysisAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:)) & ~isnan(az) & ~isnan(r);
            [A, l, m] = Zernike.getAll(l_max, m_max, az(id_ok), r(id_ok));
            if nargin == 6 && ~isempty(max_reg) && max_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
            end
            z_par = (A'*A + diag(reg_fun .* ones(size(A, 2), 1))) \ (A' * data(id_ok));    
        end
        
        function [z_par, l, m, A] = zroAnalysisAll(l_max, m_max, az, r, data, S)
            % Get Zernike polynomials parameters form reothonormalize
            % function
            %
            % SINTAX
            %   [z_par, l, m, A] = analysisAll(l_max, m_max, az, el, data, S)
            
            id_ok = ~isnan(data(:));
            [A, l, m] = Zernike.getAll(l_max, m_max, az(id_ok), r(id_ok));
            N = A'*A;
            B = A'*data;
            N = S'*N*S;
            B = S'*B;
            z_par = N\B;
            z_par = S*z_par;
        end
        
        function [z_par, l, m, A] = zAnalysis(l, m, az, r, data, max_reg)
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
            %   [z_par, l, m, A] = analysisAll(l_max, m_max, az, el, data, <max_reg = 1>)
            
            id_ok = ~isnan(data(:));
            A = Zernike.get(l, m, az(id_ok), r(id_ok));
            if nargin == 6 && ~isempty(flag_reg) && flag_reg
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
            else
                reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
            end
            z_par = (A' * A + reg_fun .* diag(ones(size(A, 2), 1))) \ (A' * data(id_ok));
        end

        function [z_par, l, m] = analysisAllBlock(l_max, m_max, az, r, data, max_reg)
            % Get Zernike polynomials parameters using maximum 1GB a block for A matrix
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
            %   [z_par, l, m] = analysisAllBlock(l_max, m_max, az, el, data, <max_reg = 1>)
            
            max_ram = 512 * 1024^2; % max 1 GB
            
            id_ok = ~isnan(data(:));
            n_par_max = (l_max+1)*(m_max+2) / 2;
            n_obs = sum(id_ok);
            n_step = ceil(8 * n_obs * n_par_max / max_ram);
            
            % Split LS preparation in n_step
            if n_step > 1
                id_ok = find(id_ok);
                id_start = 1;
                w_bar = Go_Wait_Bar.getInstance(n_step, 'Block-wise Zernike analysis');
                w_bar.createNewBar
                for i = 1 : n_step
                    %fprintf('%s%3d / %3d\n', char(8 * ones(1, 10, 'uint8')), i, n_step);  drawnow;
                    id_end = min(n_obs, ceil(n_obs * i / n_step));
                    id_tmp = id_ok(id_start : id_end);
                    [A, l, m] = Zernike.getAll(l_max, m_max, az(id_tmp), r(id_tmp));
                    if i == 1
                        % init par
                        N = zeros(numel(l));
                        TN = zeros(numel(l), 1);
                    end
                    N = N + A'*A;
                    TN = TN + A' * data(id_tmp);
                    
                    id_start = id_end + 1;
                    w_bar.goTime();
                end
                if nargin == 6 && ~isempty(max_reg) && max_reg
                    reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
                else
                    reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
                end
                z_par = (N + diag(reg_fun .* ones(size(A, 2), 1))) \ TN;
                w_bar.close();
            else
                [A, l, m] = Zernike.getAll(l_max, m_max, az(id_ok), r(id_ok));
                if nargin == 6 && ~isempty(max_reg) && max_reg
                    reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5) * max_reg;
                else
                    reg_fun = 2 * ((1./(1 + exp(-l))) - 0.5);
                end
                z_par = (A'*A + diag(reg_fun .* ones(size(A, 2), 1))) \ (A' * data(id_ok));
            end
        end
        
        
        function [l, m] = getAllDegrees(l_max, m_max)
            % Get alll the valid degrees till l_max, m_max
            % 
            % SYNTAX
            %   [l, m] = Zernike.getAllDegrees(l_max, m_max)
            
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
        end
        
            
        function [z, l, m] = getAll(l_max, m_max, az, r)
            % Generate all the Zernike parameters combinations
            %
            % INPUT
            %   l_max   maximum degree
            %   m_max   maximum order
            %   az      list of azimuth angles   [rad]
            %   r       disc radius              [0 1]
            %
            % SINTAX
            %   [z, l, m] = Zernike.getAll(l_max, az, r)
            %   [z, l, m] = Zernike.getAll(l_max, m_max, az, r)
            if nargin == 3
                r = az;
                az = m_max;
                m_max = l_max;
            end
            
            [l, m] = Zernike.getAllDegrees(l_max, m_max);
            
            switch Zernike.getMode
                case -1,   z = zernfun(l, m, r(:), az(:));
                case 0,    z = Zernike.getRecursive(l, m, r(:), az(:));
                case 1,    z = Zernike.getDiscr(l, m, r(:), az(:));
                otherwise, z = Zernike.getDiscr2(l, m, r(:), az(:));
            end
                    
        end        
        
        function z = get(l, m, az, r)
            % Get Zernike values for the polynomials
            %
            % INPUT 
            %   l       list of degrees
            %   m       list of orders
            %   az      list of azimuth angles   [rad]
            %   r       radius [0..1]
            %
            % SINTAX
            %   z = get(l, m, az, r)
                                    
            switch Zernike.getMode
                case -1,   z = Zernike.zernfun(l, m, r(:), az(:));
                case 0,    z = Zernike.getRecursive(l, m, r(:), az(:));
                case 1,    z = Zernike.getDiscr(l, m, r(:), az(:));
                otherwise, z = Zernike.getDiscr2(l, m, r(:), az(:));
            end
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
            
            x = -1 : 0.0025 : 1;
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
            %z(idx) = zernfun(l, m, r_zern(idx), 2 * pi - theta(idx) + pi/2) * z_par;
            max_ram = 1 * 1024^3; % 1GB
            idx = find(idx);
            nobs_max = (max_ram / (numel(z_par) * 8));
            lim = [unique(floor([1 : nobs_max : numel(idx)-1]))' unique(floor([(nobs_max + 1) : nobs_max : numel(idx) numel(idx)]))'];
            for i = 1 : (size(lim, 1))
                i
                z(idx(lim(i,1):lim(i,2))) = Zernike.get(l, m, 2 * pi - theta(idx(lim(i,1):lim(i,2))) + pi/2, r_zern(idx(lim(i,1):lim(i,2)))) * z_par;
            end
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
    end
    
    % Actual Zernike polynomials implementation
    methods (Static)
        function z = zernfun(n,m,r,theta,nflag)
            %ZERNFUN Zernike functions of order N and frequency M on the unit circle.
            %   Z = ZERNFUN(N,M,R,THETA) returns the Zernike functions of order N
            %   and angular frequency M, evaluated at positions (R,THETA) on the
            %   unit circle.  N is a vector of positive integers (including 0), and
            %   M is a vector with the same number of elements as N.  Each element
            %   k of M must be a positive integer, with possible values M(k) = -N(k)
            %   to +N(k) in steps of 2.  R is a vector of numbers between 0 and 1,
            %   and THETA is a vector of angles.  R and THETA must have the same
            %   length.  The output Z is a matrix with one column for every (N,M)
            %   pair, and one row for every (R,THETA) pair.
            %
            %   Z = ZERNFUN(N,M,R,THETA,'norm') returns the normalized Zernike
            %   functions.  The normalization factor sqrt((2-delta(m,0))*(n+1)/pi),
            %   with delta(m,0) the Kronecker delta, is chosen so that the integral
            %   of (r * [Znm(r,theta)]^2) over the unit circle (from r=0 to r=1,
            %   and theta=0 to theta=2*pi) is unity.  For the non-normalized
            %   polynomials, max(Znm(r=1,theta))=1 for all [n,m].
            %
            %   The Zernike functions are an orthogonal basis on the unit circle.
            %   They are used in disciplines such as astronomy, optics, and
            %   optometry to describe functions on a circular domain.
            %
            %   The following table lists the first 15 Zernike functions.
            %
            %       n    m    Zernike function             Normalization
            %       ----------------------------------------------------
            %       0    0    1                              1/sqrt(pi)
            %       1    1    r * cos(theta)                 2/sqrt(pi)
            %       1   -1    r * sin(theta)                 2/sqrt(pi)
            %       2    2    r^2 * cos(2*theta)             sqrt(6/pi)
            %       2    0    (2*r^2 - 1)                    sqrt(3/pi)
            %       2   -2    r^2 * sin(2*theta)             sqrt(6/pi)
            %       3    3    r^3 * cos(3*theta)             sqrt(8/pi)
            %       3    1    (3*r^3 - 2*r) * cos(theta)     sqrt(8/pi)
            %       3   -1    (3*r^3 - 2*r) * sin(theta)     sqrt(8/pi)
            %       3   -3    r^3 * sin(3*theta)             sqrt(8/pi)
            %       4    4    r^4 * cos(4*theta)             sqrt(10/pi)
            %       4    2    (4*r^4 - 3*r^2) * cos(2*theta) sqrt(10/pi)
            %       4    0    6*r^4 - 6*r^2 + 1              sqrt(5/pi)
            %       4   -2    (4*r^4 - 3*r^2) * sin(2*theta) sqrt(10/pi)
            %       4   -4    r^4 * sin(4*theta)             sqrt(10/pi)
            %       ----------------------------------------------------
            %
            %   Example 1:
            %
            %       % Display the Zernike function Z(n=5,m=1)
            %       x = -1:0.01:1;
            %       [X,Y] = meshgrid(x,x);
            %       [theta,r] = cart2pol(X,Y);
            %       idx = r<=1;
            %       z = nan(size(X));
            %       z(idx) = zernfun(5,1,r(idx),theta(idx));
            %       figure
            %       pcolor(x,x,z), shading interp
            %       axis square, colorbar
            %       title('Zernike function Z_5^1(r,\theta)')
            %
            %   Example 2:
            %
            %       % Display the first 10 Zernike functions
            %       x = -1:0.01:1;
            %       [X,Y] = meshgrid(x,x);
            %       [theta,r] = cart2pol(X,Y);
            %       idx = r<=1;
            %       z = nan(size(X));
            %       n = [0  1  1  2  2  2  3  3  3  3];
            %       m = [0 -1  1 -2  0  2 -3 -1  1  3];
            %       Nplot = [4 10 12 16 18 20 22 24 26 28];
            %       y = zernfun(n,m,r(idx),theta(idx));
            %       figure('Units','normalized')
            %       for k = 1:10
            %           z(idx) = y(:,k);
            %           subplot(4,7,Nplot(k))
            %           pcolor(x,x,z), shading interp
            %           set(gca,'XTick',[],'YTick',[])
            %           axis square
            %           title(['Z_{' num2str(n(k)) '}^{' num2str(m(k)) '}'])
            %       end
            %
            %   See also ZERNPOL, ZERNFUN2.
            
            %   Paul Fricker  2/28/2012 original version at https://www.mathworks.com/company/newsletters/articles/analyzing-lasik-optical-data-using-zernike-functions.html
            %                                               Paul Fricker (2020). Zernike polynomials (https://www.mathworks.com/matlabcentral/fileexchange/7687-zernike-polynomials), MATLAB Central File Exchange. Retrieved April 6, 2020.
            %   Later modified for minor speedups by 
            %   Andrea Gatti 30/11/2019
            
            % Check and prepare the inputs:
            % -----------------------------
            if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
                error('zernfun:NMvectors','N and M must be vectors.')
            end
            
            if length(n)~=length(m)
                error('zernfun:NMlength','N and M must be the same length.')
            end
            
            n = n(:);
            m = m(:);
            if any(mod(n-m,2))
                error('zernfun:NMmultiplesof2', ...
                    'All N and M must differ by multiples of 2 (including 0).')
            end
            
            if any(m>n)
                error('zernfun:MlessthanN', ...
                    'Each M must be less than or equal to its corresponding N.')
            end
            
            if any( r>1 | r<0 )
                error('zernfun:Rlessthan1','All R must be between 0 and 1.')
            end
            
            if ( ~any(size(r)==1) ) || ( ~any(size(theta)==1) )
                error('zernfun:RTHvector','R and THETA must be vectors.')
            end
            
            r = r(:);
            theta = theta(:);
            length_r = length(r);
            if length_r~=length(theta)
                error('zernfun:RTHlength', ...
                    'The number of R- and THETA-values must be equal.')
            end
            
            % Check normalization:
            % --------------------
            if nargin==5 && ischar(nflag)
                isnorm = strcmpi(nflag,'norm');
                if ~isnorm
                    error('zernfun:normalization','Unrecognized normalization flag.')
                end
            else
                isnorm = false;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute the Zernike Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Determine the required powers of r:
            % -----------------------------------
            m_abs = abs(m);
            rpowers = [];
            for j = 1:length(n)
                rpowers = [rpowers m_abs(j):2:n(j)];
            end
            rpowers = unique(rpowers);
            
            % original implementation
            % Pre-compute the values of r raised to the required powers,
            % and compile them in a matrix:
            % -----------------------------
            if (rpowers(end)) ~= (numel(rpowers) - 1)
                if rpowers(1)==0
                    rpowern = arrayfun(@(p)r.^p,rpowers(2:end),'UniformOutput',false);
                    rpowern = cat(2,rpowern{:});
                    rpowern = [ones(length_r,1) rpowern];
                else
                    rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
                    rpowern = cat(2,rpowern{:});
                end
            else %faster approach
                % new lines
                rpowern = ones(size(r, 1), rpowers(end) + 1);
                if rpowers(end) > 1
                    rpowern(:, 2) = r;
                    for p = 3 : rpowers(end) + 1
                        rpowern(:, p) = rpowern(:, p - 1) .* r;
                    end
                end
                % end of new lines
            end
            
            
            % Compute the values of the polynomials:
            % --------------------------------------
            if length_r > 1e5                 % <= original implementation
                z = zeros(length_r,length(n));
                for j = 1:length(n)
                    s = 0:(n(j)-m_abs(j))/2;
                    pows = n(j):-2:m_abs(j);
                    for k = length(s):-1:1
                        p = (1-2*mod(s(k),2))* ...
                            prod(2:(n(j)-s(k)))/              ...
                            prod(2:s(k))/                     ...
                            prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                            prod(2:((n(j)+m_abs(j))/2-s(k)));
                        idx = (pows(k)==rpowers);
                        z(:,j) = z(:,j) + p * rpowern(:,idx);
                    end
                    
                    if isnorm
                        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1)/pi);
                    end
                end
            else
                % for small number of points this is faster
                z = zeros(length_r,length(n));
                for j = 1:length(n)
                    s = 0:(n(j)-m_abs(j))/2;
                    pows = n(j):-2:m_abs(j);
                    p_sum = zeros(size(rpowers))';
                    for k = length(s):-1:1
                        p = (1-2*mod(s(k),2))* ...
                            prod(2:(n(j)-s(k)))/              ...
                            prod(2:s(k))/                     ...
                            prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                            prod(2:((n(j)+m_abs(j))/2-s(k)));
                        idx = (pows(k)==rpowers);
                        p_sum(idx) = p_sum(idx) + p;
                    end
                    z(:,j) = sum(rpowern * p_sum, 2);
                    %z(:,j) = sum(fliplr(rpowern) * flipud(p_sum), 2);
                    
                    if isnorm
                        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1)/pi);
                    end
                end
            end
            
            % END: Compute the Zernike Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the Zernike functions:
            % ------------------------------
            idx_pos = m > 0;
            idx_neg = m < 0;
            
            m_abs = abs(m);
            m_list = 0 : max(m_abs);
            if any(idx_pos) || any(idx_neg)
                tt = [cos(theta * m_list) sin(theta * m_list)];
                m_abs(idx_neg) = m_abs(idx_neg) + max(m_abs) + 1;
                z = z .* tt(:, m_abs + 1);
            end
            
            % m_list = 1 : max(m_abs);
            % if any(idx_pos)
            %     ct = cos(theta * m_list);                                % <= new line
            %     z(:,idx_pos) = z(:,idx_pos) .* ct(:, m_abs(idx_pos));    % <= new line
            %     % z(:,idx_pos) = z(:,idx_pos).*cos(theta*m_abs(idx_pos)'); <= original implementation
            % end
            % if any(idx_neg)
            %     st = sin(theta * m_list);                                % <= new line
            %     z(:,idx_neg) = z(:,idx_neg) .* st(:, m_abs(idx_neg));    % <= new line
            %     % z(:,idx_neg) = z(:,idx_neg).*sin(theta*m_abs(idx_neg)'); <= original implementation
            % end
            
            % EOF zernfun
        end
        
        function z = getDiscr(l_list, m_list, r, theta, n_r_step)
            % Part of code extrapolated from:
            %  - https://scicomp.stackexchange.com/questions/30705/numerical-stability-of-higher-order-zernike-polynomials/30713#30713
            %
            % Implemented following this paper:
            %   https://www.researchgate.net/publication/255789472_Recursive_formula_to_compute_Zernike_radial_polynomials
            %   Honarvar, Barmak & Paramesran, Raveendran. (2013). Recursive formula to compute Zernike radial polynomials. Optics letters. 38. 2487-9. 10.1364/OL.38.002487.
            
            % Maximum degree to be computed
            l_max = max(l_list);
            
            % Discretization steps
            if nargin < 5 || isempty(n_r_step)
                n_r_step = 1000;
            end
            rho = [0 : n_r_step] / n_r_step;
            
            % Get Radial Zernike polynomial
            R = Zernike.getRadialPoly(l_max, rho);         
                        
            id_r = round(r * n_r_step) + 1;

            z = zeros(numel(r), numel(l_list));
            % Extract requested degrees: m, n
            for i = 1 : numel(l_list)
                z(:,i) = R{abs(m_list(i))+1, l_list(i)+1}(id_r);
            end
            clear R;
            
            % Apply azimuth dependence
            idx_pos = m_list > 0;
            idx_neg = m_list < 0;
            
            m_abs = abs(m_list);
            m_list = 0 : max(m_abs);
            if any(idx_pos) || any(idx_neg)
                tt = [cos(theta * m_list(:)') sin(theta * m_list(:)')];
                m_abs(idx_neg) = m_abs(idx_neg) + max(m_abs) + 1;
                for i = 1 : numel(m_abs)
                    z(:,i) = z(:,i) .* tt(:, m_abs(i) + 1);
                end
            end
        end
        
        function z = getDiscr2(l_list, m_list, r, az, n_r_step, n_t_step)
            % Part of code extrapolated from:
            %  - https://scicomp.stackexchange.com/questions/30705/numerical-stability-of-higher-order-zernike-polynomials/30713#30713
            %
            % Implemented following this paper:
            %   https://www.researchgate.net/publication/255789472_Recursive_formula_to_compute_Zernike_radial_polynomials
            %   Honarvar, Barmak & Paramesran, Raveendran. (2013). Recursive formula to compute Zernike radial polynomials. Optics letters. 38. 2487-9. 10.1364/OL.38.002487.
                        
            % Maximum degree to be computed
            z_cache = Zernike.getCache;
            
            l_max = max(l_list);
            
            % Discretization steps
            if nargin < 5 || isempty(n_r_step)
                n_r_step = 1000;
            end
            
            if nargin < 5 || isempty(n_t_step)
                n_t_step = 360/(0.1/2);
            end
            
            if isempty(z_cache) || (z_cache.l_max < l_max) || isempty(z_cache.R) || isempty(z_cache.T)
                rho = linspace(0, 1, n_r_step + 1);
                        
                % Get Radial Zernike polynomial
                R = Zernike.getRadialPoly(l_max, rho);
                                             
                theta = linspace(0, 360, n_t_step + 1)' / 180 * pi;
                m_full_list = 0 : l_max;
                T = [cos(theta * m_full_list(:)') sin(theta * m_full_list(:)')];
                z_cache = struct('l_max', l_max, 'n_r_step', n_r_step, 'n_t_step', n_t_step, 'T', T);
                z_cache.R = R; % R being a cell must be added in this way or MATLAB will generate a structure matrix
                
                Zernike.setCache(z_cache); % Set the cache in the persistent object
            end
            
            id_r = round(r * z_cache.n_r_step) + 1;
            
            % Extract requested degrees: m, n
            %tmp = cell2mat(z_cache.R(abs(m_list) + 1 + (size(z_cache.R, 1) .* l_list)));
            %z = tmp(:,id_r)';
            z = zeros(numel(r), numel(l_list));
            for i = 1 : numel(l_list)
                z(:,i) = z_cache.R{abs(m_list(i))+1, l_list(i)+1}(id_r);
            end
            % Apply azimuth dependence
            idx_pos = m_list > 0;
            idx_neg = m_list < 0;
            
            m_abs = abs(m_list);
            if any(idx_pos) || any(idx_neg)
                id_t = round(mod(az, 2*pi) / (2 * pi) * z_cache.n_t_step) + 1;
                m_abs(idx_neg) = m_abs(idx_neg) + z_cache.l_max + 1;
                tmp = z_cache.T(id_t, :);
                for i = 1 : numel(m_abs)
                    z(:,i) = z(:,i) .* tmp(:, m_abs(i) + 1);
                end
            end            
        end
        
        function z = getRecursive(l_list, m_list, rho, theta)
            
            % Maximum degree to be computed
            l_max = max(l_list);
            
            R = Zernike.getRadialPoly(l_max, rho);
                        
            z = zeros(numel(rho), numel(l_list));
            % Extract requested degrees: m, n
            for i = 1 : numel(l_list)
                z(:,i) = R{abs(m_list(i))+1, l_list(i)+1};
            end
            clear R;
            
            idx_pos = m_list > 0;
            idx_neg = m_list < 0;
                        
            m_list_bk = m_list;
            m_abs = abs(m_list);
            m_list = 0 : max(m_abs);
            if any(idx_pos) || any(idx_neg)
                tt = [cos(theta * m_list(:)') sin(theta * m_list(:)')];
                m_abs(idx_neg) = m_abs(idx_neg) + max(m_abs) + 1;
                for i = 1 : numel(m_abs)
                    z(:,i) = z(:,i) .* tt(:, m_abs(i) + 1);
                end
            end
        end
        
        function R = getRadialPoly(l_max, rho)
            % Part of code extrapolated from:
            %  - https://scicomp.stackexchange.com/questions/30705/numerical-stability-of-higher-order-zernike-polynomials/30713#30713
            %
            % Implemented following this paper:
            %   https://www.researchgate.net/publication/255789472_Recursive_formula_to_compute_Zernike_radial_polynomials
            %   Honarvar, Barmak & Paramesran, Raveendran. (2013). Recursive formula to compute Zernike radial polynomials. Optics letters. 38. 2487-9. 10.1364/OL.38.002487.
            
            rho = rho(:);
            rho_x_2 = 2 * rho;
            
            R = cell(l_max+1,l_max+1);
            R{0+1,0+1} = ones(numel(rho),1);                 % R^0_0  Unfortunately zero based cell indexing is not possible
            R{1+1,1+1} = R{0+1,0+1} .* rho;                   % R^1_1  ==>  R{...+1,...+1} etc.
            
            for n = 2:l_max
                if bitget(n,1) == 0                   % n is even
                    R{0+1,n+1} = -R{0+1, n-2+1} + rho_x_2 .* R{1+1, n-1+1};                  % R^0_n
                    m_lo = 2;
                    m_hi = n-2;
                else
                    m_lo = 1;
                    m_hi = n-1;
                end
                for m = m_lo:2:m_hi
                    R{m+1,n+1} = rho .* (R{m-1+1, n-1+1} + R{m+1+1, n-1+1}) - R{m+1,n-2+1};  % R^m_n
                end
                R{n+1,n+1} = rho .*R {n-1+1 , n-1+1};                                          % R^n_n
            end
        end        
    end
    
    methods (Static)
        function fh = showAllZernike(l, m, z_par, el_min, flag_hlist)
            % Show 3D plot of Zernike polynomials 
            %
            % INPUT 
            %   l, m        array of degrees to display
            %   z_par       value of the scaling factor to use
            %   el_min      cut_off
            %   flag_hlist  set to true to display expansion in horizontal
            %
            % SINTAX
            %   fh = showAllZernike(l, m, z_par)
            %
            % EXAMPLE
            %   [l, m] = Zernike.getAllDegrees(3, 3);
            %   Zernike.showAllZernike(l, m, ones(size(l)), 0); 

                       
            fh = figure(); Core_UI.beautifyFig(fh, 'light');
            
            %%% INTERNAL PARAMETER
            scale = 1;
            %%%
            
            if nargin < 5
                flag_hlist = false;
            end
            if numel(z_par) == 1
                z_par = ones(size(l)) * z_par;
            end
            
            x = -1 : 0.005 : 1;
            y = x;
            [X,Y] = meshgrid(x,y);
            [theta, r_prj] = cart2pol(X,Y); % This radius is the correct one for my polar projection             
            r_zern = r_prj;
            if nargin >= 4 && ~isempty(el_min)
                r_max = 1 - (2 * el_min / pi);
                idx = r_prj <= r_max;
            else
                idx = r_prj <= 1;
            end
            
            if flag_hlist
                dx = 2.2;
                for i = 1 : numel(l)
                    z = nan(size(X));
                    z(idx) = Zernike.get(l(i), m(i), r_zern(idx), theta(idx)) * z_par(i);
                    
                    h = imagesc(x + dx * i, y, z); hold on;
                    h.AlphaData = ~isnan(z);
                end
                xlim([-1 (dx*i+1)]);                
            else
                dx = 1.4; dy = 2.4; % All attached
                for i = 1 : numel(l)
                    z = nan(size(X));
                    z(idx) = zernfun(l(i), m(i), r_zern(idx), theta(idx)) * z_par(i);
                    
                    h = imagesc(x + dx * m(i), y + dy * -l(i),z); hold on;
                    h.AlphaData = ~isnan(z);
                    ylim(dy*[-1 0] * max(abs(m)) + [-1 1]);
                    xlim(dx*[-1 1] * max(abs(l)) + [-1 1]);                %dx = 1.2; dy = 1.6; % All attached
                    if i == 1
                        axis equal
                        ax = gca;
                        ax.YDir = 'normal';
                    end
                    lbl = text(dx * m(i) - 0.2, - dy * l(i) - 1.4, 0, sprintf('Z^{%d}_{%d}', m(i), l(i)), 'FontName', 'Serif', 'FontWeight', 'normal');
                end                
            end

            ax = gca;
            ax.YDir = 'normal';
            axis equal
            axis off
            set(gcf,'color','w');
            colormap(flipud(Cmap.get('RdBu', 1024)));
            cb = colorbar('FontName', 'Open Sans', 'FontWeight', 'normal', 'FontSize', 16);

            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh);
        end

        function fh = showMappingFunction(n, m)
            % Show the available elevation mapping functions
            % And the relative Zernike polynomial for az = 0
            %
            % SYNTAX
            %   fh = Zernike.showMappingFunction(n, m)
            if nargin == 0
                m = 0;
                n = 100;
            end
            
            el = 0:0.1:90;
            
            fh = figure;
            hold on;
            for f = -1:3
                fun = Zernike.getElFun(f);
                rho = fun(el/180*pi);
                subplot(1,2,1);
                plot(el, rho, 'LineWidth', 2); hold on;
                subplot(1,2,2);
                plot(el, f+Zernike.get(n, m, rho * 0, rho), 'LineWidth', 2); hold on;
            end
            %title(sprintf('Z_{%d}^{%d} \\fontsize{5} \n', m, n), 'FontName', 'Serif', 'FontSize', 14);
            subplot(1,2,1);
            title(sprintf('Mapping function \\fontsize{5} \n'), 'FontName', 'Serif', 'FontSize', 14);
            subplot(1,2,1);
            xlabel('Elevation [deg]');
            ylabel('Radius');
            legend({...
                'mp\_mode -1) (90 - el) / 90', ...
                'mp\_mode  0)  cos(el)', ...
                'mp\_mode  1)  cos(el)^2', ...
                'mp\_mode  2)  sin(90 * cos(el)^2)', ...
                'mp\_mode  3)  sin(90 * cos(el))'}, 'FontName', 'Serif', 'FontSize', 12, 'location', 'NorthEast');
            grid on;
            xlim([0 90]);
            subplot(1,2,2);
            title(sprintf('Zernike polynomial \\fontsize{5} \n'), 'FontName', 'Serif', 'FontSize', 14);
            xlabel('Elevation [deg]');
            
            Core_UI.beautifyFig(fh, 'light');
            
            ylabel(sprintf('Z_{%d}^{%d} \\fontsize{5} \n', m, n), 'FontName', 'Serif', 'FontSize', 16);
            grid on;
            xlim([0 90]);
            
            Core_UI.addExportMenu
        end
    end
    
    % Copy of methods from Core_Utils to make Zernike Class independent
    methods (Static, Access = private)        
        function [data_map, n_map, az_grid_out, el_grid_out] = hemiGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells, n_min)
            % Grid points on a regularly gridded semi sphere
            %
            % INPUT 
            %   az      azimuth
            %   el      elevation
            %
            % SYNTAX
            %   [data_map, n_map, az_grid, el_grid] = Zernike.hemiGridder(az, el, data, step_deg, step_deg_out, flag_congurent_cells)
            
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
        
        function [data, row, col] = hgrid2scatter(az, el, hmap)
            % Grid points on a regularly gridded semi sphere
            %
            % INPUT 
            %   az      azimuth list   [n x 1]
            %   el      elevation list [n x 1]
            %
            % SYNTAX
            %   [data, row, col] = Core_Utils.hgrid2scatter(az, el, hmap)
            
            % Define grid
            % az -180 : 180
            % el 0 : 90
            
            step_deg = fliplr([90 360] ./ size(hmap));
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
            
            % init data
            data = nan(size(az));
                        
            % fill maps
            for i = 1 : numel(data)
                data(i) = hmap(row(i), col(i));                
            end
        end        
    end
end

