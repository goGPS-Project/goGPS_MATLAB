%   CLASS Antenna
% =========================================================================
%
% DESCRIPTION
%   Collector of antenna models 
%   Container of the ATX file()
%
% EXAMPLE
%   atx = Antenna();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Antenna


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
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

classdef Antenna < handle

    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant)
    end
   
    
    %% PROPERTIES ANTENNA
    % ==================================================================================================================================================
    properties
        type     %    A20    | Antenna type: strict IGS rcvr_ant.tab antenna and radome code
        num      %    A20    | Serial number (blank: all representatives of the specified antenna type)
        
        cal_date %    A10    | Cal day
        dazi     %   F6.1    | aximuth increment
        az_grid  % Array of azimuth coordinates
        zen1     %   F6.1    | zenithal first limits for gri definition
        zen2     %   F6.1    | zenithal limits for gri definition
        dzen     %   F6.1    | zenithal increment
        n_freq   %     D1    | number of calibrated frequencies are present
        f_code   %   A3xN    | calibrated antenna frequencies, e.g. G01, G02, G05, E01, R01, ...
        
        start    % GPS_Time  | First epoch of calibration validity
        stop     % GPS_Time  | Last epoch of calibration validity
        
        sinex_code % A10     | antenna name for SINEX format
        
        pco      %   3FxN    | PCO Nort East Up for receivers XYZ for satellites
        pcv_noaz %           | no az dependent PCV
        pcv      %           | PCV
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = public)
        % Concrete implementation.  See Singleton superclass.
        function this = Antenna()
            % Core object creator
            this.init();
        end
    end    
    
    %% METHODS INIT & STATIC GETTERS & SETTERS
    % ==================================================================================================================================================
    methods (Static, Access = public)        
        function this = fromText(txt, lim, has_cr)
            % Import antenna parameters from a part of antex file 
            % that contains antenna information
            % From the first line START OF ANTENNA
            % To the laste line END OF ANTENNA
            %
            % INPUT
            %   txt     raw *char sequuence
            %   lim     optional contains start, stop, size, of each line
            %   has_cr  optional has carriage return (does it contain char 13?)
            %
            % SYNTAX
            %   Antenna.fromText
            this = Antenna();
            log = Core.getLogger();
            
            % Checking if some input to th function are missing
            if nargin < 3 || isempty(has_cr)
                if ~isempty(find(txt(1:min(1000, numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
            end
            
            if nargin < 2 || isempty(lim)
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  (numel(txt) - double(has_cr))
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                lim = [lim lim(:,2) - lim(:,1)];
                while lim(end,3) < 3
                    lim(end,:) = [];
                end
            else
                lim(:, 1:2) = lim(:, 1:2) + 1 - lim(1); % to be shure reset lim offset
            end
            % Start parsing ---------------------------------------------------------
            
            id = find(txt(lim(:,1) + 60) == 'T', 1, 'first'); % find line containg TYPE / SERIAL NO
            if ~isempty(id)
                this.type = txt(lim(id,1) + (0:19)); % antenna type
                this.num = txt(lim(id,1) + (20:39)); % antenna serial number
            end
            
            id = find(txt(lim(:,1) + 60) == 'M', 1, 'first'); % find line containg METH / BY / # / DATE
            if ~isempty(id)
                this.cal_date = GPS_Time(datenum(txt(lim(id,1) + (50:58)), 'dd-mmm-yy')); % calibration date            
            end
            
            id = find(txt(lim(:,1) + 60) == 'D', 1, 'first'); % find line containg DAZI
            if ~isempty(id)
                this.dazi = str2double(txt(lim(id,1) + (1:7)));
            else
                this.dazi = 0.0;
            end
            
            id = find(txt(lim(:,1) + 60) == 'Z', 1, 'first'); % find line containg ZEN1 / ZEN2 / DZEN
            if ~isempty(id)
                tmp = sscanf(txt(lim(id,1) + (2:19)), '%f');
                if numel(tmp) < 3
                    % Thats ot normal
                    log.addWarning(sprintf('Antenna %s have not valid ZEN1 / ZEN2 / DZEN definition', strtrim(this.type)))
                    this.zen1 = 0.0;
                    this.zen2 = 0.0;
                    this.dzen = 0.0;
                else
                    this.zen1 = tmp(1);
                    this.zen2 = tmp(2);
                    this.dzen = tmp(3);
                end
            else
                log.addWarning(sprintf('Antenna %s have missing ZEN1 / ZEN2 / DZEN definition', strtrim(this.type)))
                this.zen1 = 0.0;
                this.zen2 = 0.0;
                this.dzen = 0.0;
            end
            
            id = find(txt(lim(:,1) + 66) == 'F', 1, 'first'); % find line containg # VALID FROM
            if ~isempty(id)
                this.start = GPS_Time(txt(lim(id,1) + (0:42)), [], true);
                if this.start.isEmpty
                    % No validity start
                    this.start = GPS_Time(GPS_Time.GPS_ZERO);
                end                
            else
                % No validity start
                this.start = GPS_Time(GPS_Time.GPS_ZERO);
            end
            
            id = find(txt(lim(:,1) + 66) == 'U', 1, 'first'); % find line containg # VALID UNTIL
            if ~isempty(id)
                this.stop = GPS_Time(txt(lim(id,1) + (0:42)), [], true);
                if this.stop.isEmpty
                    % No validity stop
                    this.stop = GPS_Time('2100 01 01');
                end                
            else
                % No validity stop
                this.stop = GPS_Time('2100 01 01');
            end
            
            id = find(txt(lim(:,1) + 64) == 'X', 1, 'first'); % find line containg SINEX CODE
            if ~isempty(id)
                this.sinex_code = txt(lim(id,1) + (0:9));
            else
                this.sinex_code = '          ';
            end

            id = find(txt(lim(:,1) + 60) == '#', 1, 'first'); % find line containg # OF FREQUENCIES
            if ~isempty(id)
                this.n_freq = str2double(txt(lim(id,1) + (0:5)));
            else
                this.n_freq = 0;
            end

            id_fstart = find(txt(lim(:,1) + 61) == 'T'); % find lines containg START OF FREQUENCY
            id_pco = find(txt(lim(:,1) + 60) == 'N'); % find lines containg  NORTH / EAST / UP
            id_fstop = find(txt(lim(:,1) + 61) == 'N'); % find lines containg END OF FREQUENCY
            
            if this.n_freq ~= numel(id_fstart) || this.n_freq ~= numel(id_pco) || this.n_freq ~= numel(id_fstop)
                % Frequencies are currupted, no scanning is admissible for this antenna
                log.addWarning(sprintf('Antenna %s have missing currupted frequency definition', strtrim(this.type)))
                this.n_freq = 0;
            else
                % for each frequency                
                this.f_code = '   ';
                for f = 1 : this.n_freq
                    % get freq name
                    this.f_code(f,:) = txt(lim(id_fstart(f),1) + (3:5));
                    if ~strcmp(this.f_code(f,:), txt(lim(id_fstop(f),1) + (3:5)))
                        log.addWarning(sprintf('Antenna %s have non valid Frequency %s', this.f_code(f,:)));
                    end
                    
                    pco = nan2zero(sscanf(txt(lim(id_pco(f),1) + (0:29)), '%f')');
                    if numel(pco) ~= 3
                        pco = [0 0 0];
                    end
                    if this.dazi == 0
                        n_az = 0;
                    else
                        n_az = 360 / this.dazi + 1;
                    end
                    n_lin = id_fstop(f) - id_pco(f) - 1;
                    this.pco{f} = pco;
                    if n_lin ~= (n_az + 1)
                        % something bad happened
                        log.addWarning(sprintf('Antenna %s have non valid PCV/PCO', this.type))
                    else
                        
                        
                        n_zen = (this.zen2 - this.zen1) / this.dzen + 1; % Number of zenithal dependent values
                        % next line is NOAZI
                        pcv_noaz = sscanf(txt(lim(id_pco(f)+1,1) + (8:(8 * (n_zen + 1) - 1))), '%f')';
                        if numel(pcv_noaz) ~= n_zen
                            log.addWarning(sprintf('Antenna %s have non valid PCV noaz', this.type))
                            pcv_noaz = zeros(1, n_zen);
                        end
                        this.pcv_noaz{f} = pcv_noaz;
                        
                        if n_az > 0
                            n_char = 8 * (n_zen);
                            pcv = sscanf(serialize(txt(repmat(lim(id_pco(f)+(2:(n_az + 1)),1), 1, n_char + 1) + repmat((8:(8 * (n_zen + 1))), n_az, 1))')', '%f');
                            if isempty(this.az_grid)
                                this.az_grid = sscanf(serialize(txt(repmat(lim(id_pco(f)+(2:(n_az + 1)),1), 1, 8 + 1) + repmat((1:9), n_az, 1))')', '%f');
                            end
                            if numel(pcv) ~= (n_az * n_zen)
                                log.addWarning(sprintf('Antenna %s have non valid PCV', this.type))
                                pcv = zeros(n_az, n_zen);
                            else
                                pcv = reshape(pcv, n_zen, n_az)';
                            end
                            this.pcv{f} = pcv;
                        end
                    end
                end                    
            end
        end        
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods        
        function init(this)
            % Empty the antenna properties
            %
            % SYNTAX
            %   this.init();
            
            this.type = '';
            this.num = '';
            
            this.cal_date = GPS_Time(GPS_Time.GPS_ZERO);
            this.dazi = 0.0;
            this.zen1 = 0.0;
            this.zen2 = 0.0;
            this.dzen = 0.0;
            this.n_freq = 0;
            this.f_code = '';
            
            this.start = GPS_Time(GPS_Time.GPS_ZERO);
            this.stop = GPS_Time('2100 01 01'); % Any date in the future
            
            this.sinex_code = '';
            
            this.pco = {};
            this.pcv_noaz = {};
            this.pcv = {};            
        end
        
        function import(this, type, num, cal_date, dazi, az_grid, zen1, zen2, dzen, n_freq, f_code, start, stop, sinex_code, pco, pcv_noaz, pcv)
            % Function to import all the parameters to fill the antenna
            %
            % INPUT
            % type       char [20]          - type
            % num        char [20]          - num
            % cal_data   GPS_Time           - calibration date
            % dazi       double             - delta azimuth
            % az_grid    double array       - azimut grid values
            % zen1       double             - first zenith value
            % zen2       double             - last zenith value
            % dzen       double             - delta zenith
            % n_freq     double             - number of stored frequencies
            % f_code     char [nfreq x 3]   - frequency code (e.g. G01, G02, E01, ...)
            % start      GPS_Time           - first day of validity
            % stop       GPS_Time           - last day of validiti
            % sinex_code char [20]          - sinex code
            % pco        cell array [nfreq] - Phase Center Offset
            % pcv_noaz   cell array [nfreq] - Phase Center Variation non azimuth dependent
            % pcv        cell_array [nfreq] - Phase Center Variation azimuth dependent (table)
            %
            % SYNTAX
            %  this.import(type, num, cal_data, dazi, az_grid, zen1, zen2, dzen, n_freq, f_code, start, stop, sinex_code, pco, pcv_noaz, pcv)
            this.type = type;
            this.num = num;
            this.cal_date = cal_date;
            this.dazi = dazi;
            this.az_grid = az_grid;
            this.zen1 = zen1;
            this.zen2 = zen2;
            this.dzen = dzen;
            this.n_freq = n_freq;
            this.f_code = f_code;
            this.start = start;
            this.stop = stop;
            this.sinex_code = sinex_code;
            this.pco = pco;
            this.pcv_noaz = pcv_noaz;
            this.pcv = pcv;
        end
    end
    
    %% METHODS UTILITIES
    % ==================================================================================================================================================
    methods
        function is_empty = isEmpty(this)
            % Return true if the antenna does not contains values for PCV/PCO
            %
            % SYNTAX
            %   this.isEmpty
            is_empty = true(size(this(:), 1), 1);
            for a = 1 : size(this(:), 1)
                is_empty(a) = isempty(this(a).pco) && isempty(this(a).pcv) && isempty(this(a).pcv_noaz);
            end            
        end
        
        function [f_id] = getFreqId(this, f_code)
            % Get frequency id given the trequency code
            %
            % INPUT
            %   f_code  3ch frequency code e.g. G01,E12, ...
            %
            % SYNTAX
            %    f_id = this.getFreqId(f_code)
            
            f_id = find(Core_Utils.code3Char2Num(this.f_code) == Core_Utils.code3Char2Num(f_code));
            if isempty(f_id)
                % Search for the closest frequency
                log = Core.getLogger();
                cc = Core.getState.getConstellationCollector;
                % find closer GPS frequency
                gps_id = find(this.f_code(:, 1) == 'G');
                aval_frq = char(this.f_code(gps_id,3));
                fr_gps = zeros(size(aval_frq));
                for i = 1 : length(aval_frq)
                    ava_frq_idx = cc.getGPS.CODE_RIN3_2BAND == aval_frq(i);
                    fr_gps(i) = cc.getGPS.F_VEC(ava_frq_idx);
                end
                trg_fr_lid = cc.getSys(f_code(1)).CODE_RIN3_2BAND == f_code(3);
                fr_trg = cc.getSys(f_code(1)).F_VEC(trg_fr_lid);
                [~, id_near] = min(abs(fr_gps-fr_trg));
                f_id = gps_id(id_near);
                if isempty(f_id)
                    log.addWarning(sprintf('No corrections found for antenna type "%s" on frequency %s', strtrim(this.type), f_code));
                else
                    log.addWarning(sprintf('No corrections found for antenna type "%s" on frequency %s, using %s instead', strtrim(this.type), f_code, this.f_code(gps_id(id_near),:)));
                end
            end
        end
        
        function [pco, f_id] = getPCO(this, f_code)
            % Get PCO given the trequency code (NEU/XYZ)
            % 
            % INPUT
            %   f_code  3ch frequency code e.g. G01,E12, ...
            %
            % SYNTAX
            %   [pco, f_id] = this.getPCO(f_code)
            
            f_id = this.getFreqId(f_code);
            if isempty(f_id)
                pco = 0;
            else
                pco = this.pco{f_id}';
            end
        end
                
        function pcv_delay = getPCV(this, f_id, el, az, method)
            % get the pcv correction for a given satellite and a given
            % azimuth and elevations using linear or bilinear interpolation
            
            if ischar(f_id)
                f_id = this.getFreqId(f_id);
            end
            
            if isempty(f_id)
                pcv_delay = 0;
            else
                if (nargin < 5)
                    method = 'minterp';
                end
                pcv_delay = zeros(size(el));
                
                % transform el in zen
                zen = 90 - el;
                % get el idx
                zen_pcv = this.zen1 : this.dzen : this.zen2;
                
                if nargin < 4 || (numel(this.pcv) < f_id) || isempty(this.pcv{f_id}) % no azimuth change
                    pcv_val = this.pcv_noaz{f_id}; % extract the right frequency
                    %pcv_delay = d_f_r_el .* pcv_val(zen_idx)' + (1 - d_f_r_el) .* pcv_val(zen_idx + 1)';
                    
                    % Use polynomial interpolation to smooth PCV
                    pcv_delay = Core_Utils.interp1LS(zen_pcv', pcv_val', min(8,numel(pcv_val)), zen);
                else
                    % find azimuth indexes
                    az_pcv = this.az_grid';
                    min_az = az_pcv(1);
                    max_az = az_pcv(end);
                    d_az = (max_az - min_az)/(length(az_pcv)-1);
                    az_idx = min(max(floor((mod(az, 360) - min_az) / d_az) + 1, 1), length(az_pcv) - 1); % floor obs point id
                    d_f_r_az = min(max(mod(az, 360) - (az_idx-1)*d_az, 0)/d_az, 1); % normalized distance from floor obs point
                    
                    if strcmp(method, 'polyLS')
                        % Polynomial interpolation in elevation and griddedInterpolant in az
                        % filter PCV values in elevetion - slower
                        %%
                        pcv_val = this.pcv{f_id}'; % extract the right frequency
                        step = 0.25;
                        zen_interp = (((floor(min(zen)/step)*step) : step : (ceil(max(zen)/step)*step))' - this.zen1);
                        
                        % Interpolate with a 7th degree polinomial in elevation
                        pcv_val = Core_Utils.interp1LS(zen_pcv', pcv_val, min(8,size(pcv_val, 1)), zen_interp);
                        
                        n_az = size(pcv_val, 2);
                        b_dx = (n_az - floor(n_az / 3)) : (n_az - 1); % dx border
                        b_sx = 2 : (n_az - floor(n_az / 3));          % sx border
                        % replicate border for azimuth periodicity
                        pcv_val = [pcv_val(:, b_dx), pcv_val, pcv_val(:, b_sx)];
                        az_val = [az_pcv(b_dx) - 360, az_pcv, az_pcv(b_sx) + 360];
                        
                        % Using scattered interpolant for azimuth to have a smooth interpolation
                        [zen_m, az_m] = ndgrid(zen_interp, az_val);
                        fun = griddedInterpolant(zen_m, az_m, pcv_val);
                        pcv_delay = fun(zen, mod(az, 360));
                    elseif strcmp(method, 'minterp')
                        % Interpolation with griddedInterpolant
                        % linear interpolation - faster
                        pcv_val = this.pcv{f_id}'; % extract the right frequency
                        
                        n_az = size(pcv_val, 2);
                        b_dx = (n_az - floor(n_az / 3)) : (n_az - 1); % dx border
                        b_sx = 2 : (n_az - floor(n_az / 3));          % sx border
                        % replicate border for azimuth periodicity
                        pcv_val = [pcv_val(:, b_dx), pcv_val, pcv_val(:, b_sx)];
                        az_val = [az_pcv(b_dx) - 360, az_pcv, az_pcv(b_sx) + 360];
                        
                        % Using scattered interpolant for azimuth to have a smooth interpolation
                        [zen_m, az_m] = ndgrid(zen_pcv, az_val);
                        fun = griddedInterpolant(zen_m, az_m, pcv_val, 'linear');
                        pcv_delay = fun(zen, az);
                    elseif strcmp(method, 'lin')
                        
                        % linear interpolation between consecutive values
                        % bad performance
                        %%
                        zen_idx = min(max(floor((zen - this.zen1)/this.dzen) + 1 , 1),length(zen_pcv) - 1);
                        d_f_r_el = min(max(zen_idx*this.dzen - zen, 0)/ this.dzen, 1) ;
                        pcv_val = this.pcv{f_id}; % extract the right frequency
                        %interpolate along zenital angle
                        idx1 = sub2ind(size(pcv_val),az_idx,zen_idx);
                        idx2 = sub2ind(size(pcv_val),az_idx,zen_idx+1);
                        pcv_delay_lf =  d_f_r_el .* pcv_val(idx1) + (1 - d_f_r_el) .* pcv_val(idx2);
                        idx1 = sub2ind(size(pcv_val),az_idx+1,zen_idx);
                        idx2 = sub2ind(size(pcv_val),az_idx+1,zen_idx+1);
                        pcv_delay1_rg = d_f_r_el .* pcv_val(idx1) + (1 - d_f_r_el) .* pcv_val(idx2);
                        %interpolate alogn azimtuh
                        pcv_delay = (1 - d_f_r_az) .* pcv_delay_lf + d_f_r_az .* pcv_delay1_rg;
                    end
                end
            end
        end
        
        function showPCV(this)
            for f = 1 : this.n_freq
                %%
                fh = figure();
                fh.Name = sprintf('%03d: %s - %s', fh.Number, strtrim(this.type), strtrim(this.f_code(f, :))); fh.NumberTitle = 'off';
                if isempty(this.pcv)
                    pcv = repmat(this.pcv_noaz{f}', 1, 180);
                else
                    pcv = this.pcv{f}';
                end
                polarplot3d(pcv, 'PlotType', 'surf', 'RadialRange',[this.zen1 this.zen2] / 180 * pi, ...
                    'AxisLocation', 0, 'InterpMethod', 'cubic', ...
                    'PlotType', 'surfn', 'tickspacing', 15, ...
                    'GridColor', [0.7 0.7 0.7]);

                caxis(perc(abs(pcv(:)), 0.99) * [-1 1]); colormap(flipud(Cmap.get('PuOr', 256)));
                colorbar();
                axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38], ...
                          'Xlim', 1.5 * [-this.zen2 this.zen2] / 180 * pi, 'Ylim', 1.5 * [-this.zen2 this.zen2] / 180 * pi, ... 
                          % 'XTick', [], 'YTick', [], 'Color', 'none', 'XColor', 'none', 'YColor', 'none' 
                          };
                ax = gca;
                set(ax, axprop{:});
                %ax.Children(end).FaceAlpha = 0.95;
                title(sprintf('%s PCV\nFrequency %s', strtrim(this.type), strtrim(this.f_code(f, :))), 'interpreter', 'none')

            
            end
        end
    end

end
