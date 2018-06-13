%   CLASS Receiver_Output
% =========================================================================
%
% 
%   Class to store receiver outputs
%
% EXAMPLE
%   trg = Receiver_Output();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea, Giulio Tagliaferro ...
%  Contributors:
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
%--------------------------------------------------------------------------
classdef Receiver_Output < Receiver_Commons
    properties (SetAccess = public, GetAccess = public)
        parent         % habdle to parent object
    end
     % ==================================================================================================================================================
    %% PROPERTIES CELESTIAL INFORMATIONS
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        sat = struct( ...
            'outlier_idx_ph',   [], ...    % logical index of outliers
            'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
            'az',               [], ...    % double  [n_epoch x n_sat] azimuth
            'el',               [], ...    % double  [n_epoch x n_sat] elevation  
            'res',              [], ...    % residual per staellite
            'slant_td',         []  ...    % slant total delay (except ionosphere delay)
            )
    end
    % ==================================================================================================================================================
    %% PROPERTIES TROPO
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
     
        
        pressure      % pressure           double   [n_epoch x 1]
        temperature   % temperature           double   [n_epoch x 1]
        humidity      % humidity           double   [n_epoch x 1]
    end
      % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    
    methods
        function this = Receiver_Output(parent)
            this.parent = parent;
            this.init();
        end
        
        function init(this)
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.rf = Core_Reference_Frame.getInstance();
            this.w_bar = Go_Wait_Bar.getInstance();
            
            this.reset();
        end
        
        function reset(this)
            this.reset@Receiver_Commons();
                       
            this.sat = struct( ...
            'outlier_idx_ph',   [], ...    % logical index of outliers
            'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
            'az',               [], ...    % double  [n_epoch x n_sat] azimuth
            'el',               [], ...    % double  [n_epoch x n_sat] elevation  
            'res',              [], ...    % residual per staellite
            'slant_td',         []  ...    % slant total delay (except ionosphere delay)
            )
            
        end
        
    end
     % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        % standard utility
        function toString(this)
            % Display on screen information about the receiver
            % SYNTAX this.toString();
            for r = 1 : numel(this)
                if ~this(r).isempty
                    fprintf('----------------------------------------------------------------------------------\n')
                    this(r).log.addMarkedMessage(sprintf('Receiver (%d) %s', r, this(r).marker_name));
                    fprintf('----------------------------------------------------------------------------------\n')
                    this(r).log.addMessage(sprintf(' From     %s', this(r).time.first.toString()));
                    this(r).log.addMessage(sprintf(' to       %s', this(r).time.last.toString()));
                    this(r).log.newLine();
                    this(r).log.addMessage(sprintf(' Rate of the observations [s]:            %d', this(r).getRate()));
                    this(r).log.newLine();
              
                    fprintf(' ----------------------------------------------------------\n')
                    if ~isempty(this(r).xyz)
                        enu = zero2nan(this(r).xyz); [enu(:, 1), enu(:, 2), enu(:, 3)] = cart2plan(zero2nan(this(r).xyz(:,1)), zero2nan(this(r).xyz(:,2)), zero2nan(this(r).xyz(:,3)));
                        xyz_m = median(zero2nan(this(r).xyz), 1, 'omitnan');
                        enu_m = median(enu, 1, 'omitnan');
                        this(r).log.newLine();
                        this(r).log.addMessage(' Receiver median position:');
                        this(r).log.addMessage(sprintf('     X = %+16.4f m        E = %+16.4f m\n     Y = %+16.4f m        N = %+16.4f m\n     Z = %+16.4f m        U = %+16.4f m', ...
                            xyz_m(1), enu_m(1), xyz_m(2), enu_m(2), xyz_m(3), enu_m(3)));
                    end
                    fprintf(' ----------------------------------------------------------\n')
                end
            end
        end
        
        function time = getTime(this)
            % return the time stored in the object
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTime()
            
            time = this(1).time.getCopy();
        end
        
        function [pwv, time] = getPwv(this)
            % SYNTAX
            %  [pwv, time] = this.getPwv()
            
            pwv = {};
            time = {};
            for r = 1 : size(this, 2)
                time{r} = this(1, r).time.getEpoch(this(1, r).getIdSync); %#ok<AGROW>
                pwv{r} = this(1, r).pwv(this(1, r).getIdSync); %#ok<AGROW>
                
                for s = 2 : size(this, 1)
                    pwv_tmp = this(s, r).pwv(this(s, r).getIdSync);
                    time_tmp = this(s, r).time.getEpoch(this(s, r).getIdSync);
                    pwv{r} = [pwv{r}; pwv_tmp];
                    time{r} = time{r}.append(time_tmp);
                end
            end
            
            if numel(pwv) == 1
                pwv = pwv{1};
                time = time{1};
            end
        end
        
    end
    
     % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods
        function legacyImportResults(this, file_prefix, run_start, run_stop)
            % Import after reset a position and tropo file (if present)
            %
            % SYNTAX
            %   this.legacyImportResults(file_prefix, <run_start>, <run_stop>)
            %
            % INPUT
            %   file_name     it could include the key ${RUN} that will be substituted with a 3 digits number containing the run, from run_start to run_stop
            %   run_start     number of the first run to load
            %   run_stop      number of the last run to load
            %
            if (nargin == 1) || isempty(file_prefix)
                [file_prefix, file_path] = uigetfile('*.txt', 'Select a _position.txt or _tropo.txt file');
                file_prefix = [file_path file_prefix];
            end
            this.reset();
            if (length(file_prefix) > 13 && strcmp(file_prefix(end - 12 : end), '_position.txt'))
                file_prefix = file_prefix(1 : end - 13);
            end
            if (length(file_prefix) > 10 && strcmp(file_prefix(end - 9 : end), '_tropo.txt'))
                file_prefix = file_prefix(1 : end - 10);
            end
            
            if nargin < 4
                run_start = 0;
                run_stop = 0;
            end
            
            GPS_RUN = '${RUN}';
            r = 0;
            for run = run_start : run_stop
                r = r + 1;
                file_name = [strrep(file_prefix, GPS_RUN, sprintf('%03d', run)) '_position.txt'];
                marker_name = File_Name_Processor.getFileName(file_name);
                this.marker_name = marker_name(1:4);
                this.log.addMessage(this.log.indent(sprintf('Importing %s', File_Name_Processor.getFileName(file_name))));
                if exist(file_name, 'file')
                    this.legacyAppendPosition(file_name);
                    
                    file_name = [strrep(file_prefix, GPS_RUN, sprintf('%03d', run)) '_tropo.txt'];
                    if exist(file_name, 'file')
                        this.log.addMessage(this.log.indent(sprintf('Importing %s', File_Name_Processor.getFileName(file_name))));
                        this.legacyAppendTropo(file_name)
                    else
                        this.log.addMessage(sprintf('Error loading the tropo file, it does not exists'));
                    end
                else
                    this.log.addMessage(sprintf('Error loading the position file, it does not exists'));
                end
            end
            
            if isempty(this.id_sync)
                this.id_sync = 1 : this.time.length();
            end
        end
        
        function legacyAppendPosition (this, file)
            % import and append from a position file
            
            % Open position file as a string stream
            fid = fopen(file);
            txt = fread(fid,'*char')';
            txt(txt == 13) = [];
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            
            % corrupted lines
            ko_lines = find(lim(:, 3) ~= median(lim(lim(:,3)>400,3)));
            for l = numel(ko_lines) : -1 : 1
                txt(lim(ko_lines(l), 1) : lim(ko_lines(l), 2) + 1) = [];
            end
            
            % importing header informations
            eoh = 1; % end of header (the header of position files contains only 1 line)
            % File example:
            %    Date        GPS time           GPS week          GPS tow         Latitude        Longitude      h (ellips.)           ECEF X           ECEF Y           ECEF Z        UTM North         UTM East      h (orthom.)         UTM zone        Num. Sat.             HDOP            KHDOP      Local North       Local East          Local H    Ambiguity fix     Success rate              ZTD              ZWD              PWV
            %2017/04/03    00:00:00.000             1943        86400.000      45.80216141       9.09562643         291.5094     4398305.8406      704150.1081     4550153.9697     5072071.0952      507430.9212         244.4506             32 T               10            0.870            0.472           0.0000           0.0000           0.0000                0           0.0000          2.30187          0.04934          0.00000
            data = sscanf(txt(lim(2,1):end)','%4d/%2d/%2d    %2d:%2d:%6f             %4d %16f %16f %16f %16f %16f %16f %16f %16f %16f %16f %14d %c %16d %16f %16f %16f %16f %16f %16d %16f %16f %16f %16f\n');
            data = reshape(data, 30, numel(data)/30)';
            % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(data(:,1:6));
            else
                this.time.append6ColDate(data(:,1:6));
            end
            
            lat = data(:,9);
            lon = data(:,10);
            h_ellips = data(:,11);
            h_ortho = data(:,17);
            
            xyz = data(:,12:14);
            enu = data(:,[16 15 17]);
            n_sat = data(:,20);
            hdop =  data(:,21);
            khdop = data(:,22);
            a_fix = data(:,26);
            s_rate = data(:,27);
            
            % Append in obj
            this.xyz = [this.xyz; xyz];
            this.enu = [this.enu; enu];
            
            this.lat = [this.lat; lat];
            this.lon = [this.lon; lon];
            this.h_ellips = [this.h_ellips; h_ellips];
            this.h_ortho = [this.h_ortho; h_ortho];
            
            this.n_sat = [this.n_sat; n_sat];
            this.hdop = [this.hdop; hdop];
            this.khdop = [this.khdop; khdop];
            this.a_fix = [this.a_fix; a_fix];
            this.s_rate = [this.s_rate; s_rate];
            
            clear data;
        end
        
        function legacyAppendTropo (this, file)
            % import and append from a tropo file
            
            % Open tropo file as a string stream
            fid = fopen(file);
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            
            % importing header informations
            eoh = 1; % end of header (the header of tropo files contains only 1 line)
            
            % list and count satellites in view and not in view
            s = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d(?=-az)', 'match');
            num_sat = length(s);
            this.sat.id = reshape(cell2mat(s), 3, num_sat)';
            
            % extract all the epoch lines
            string_time = txt(repmat(lim(2:end,1),1,26) + repmat(0:25, size(lim,1)-1, 1))';
            % convert the times into a 6 col time
            date = cell2mat(textscan(string_time,'%4f/%2f/%2f    %2f:%2f:%6.3f'));
            
            % import it as a GPS_Time obj
            % it should be imported from the _position file
            time = GPS_Time(date, [], 1);
            [~, id_int, id_ext] = intersect(round(this.time.getMatlabTime() * 86400 * 1e3), round(time.getMatlabTime() * 86400 * 1e3));
            n_epo = length(id_int);
            
            % extract all the ZHD lines
            string_zhd = txt(repmat(lim(2:end,1),1,17) + repmat(62:78, size(lim,1)-1, 1))';
            tmp = sscanf(string_zhd,'%f')'; clear string_zhd
            this.apr_zhd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.apr_zhd(id_int, 1) = tmp(id_ext);
            
            % extract all the ZTD lines
            string_ztd = txt(repmat(lim(2:end,1),1,17) + repmat(78:94, size(lim,1)-1, 1))';
            tmp = sscanf(string_ztd,'%f'); clear string_ztd
            this.ztd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.ztd(id_int, 1) = tmp(id_ext);
            
            % extract all the TGN lines
            string_tgn = txt(repmat(lim(2:end,1),1,17) + repmat(95:111, size(lim,1)-1, 1))';
            tmp = sscanf(string_tgn,'%f'); clear string_tgn
            this.tgn(end + 1 : size(this.xyz, 1), 1) = nan;
            this.tgn(id_int, 1) = tmp(id_ext);
            
            % extract all the TGE lines
            string_tge = txt(repmat(lim(2:end,1),1,17) + repmat(112:128, size(lim,1)-1, 1))';
            tmp = sscanf(string_tge,'%f'); clear string_tge
            this.tge(end + 1 : size(this.xyz, 1), 1) = nan;
            this.tge(id_int, 1) = tmp(id_ext);
            
            % extract all the ZWD lines
            string_zwd = txt(repmat(lim(2:end,1),1,17) + repmat(129:145, size(lim,1)-1, 1))';
            tmp = sscanf(string_zwd,'%f'); clear string_zwd
            this.zwd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.zwd(id_int, 1) = tmp(id_ext);
            
            % extract all the PWV lines
            string_pwv = txt(repmat(lim(2:end,1),1,17) + repmat(146:162, size(lim,1)-1, 1))';
            tmp = sscanf(string_pwv,'%f'); clear string_pwv
            this.pwv(end + 1 : size(this.xyz, 1), 1) = nan;
            this.pwv(id_int, 1) = tmp(id_ext);
            
            %  extract all STD values if present
            slant_start = regexp(txt(lim(1,1) : lim(1,2)),'STD') - 6;
            num_sat = numel(slant_start);
            this.sat.slant_td(end + 1 : size(this.xyz, 1), 1 : num_sat) = nan;
            for s = 1 : numel(slant_start)
                tmp = sscanf(txt(bsxfun(@plus, repmat(slant_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1) - 1))', '%f');
                this.sat.slant_td(id_int, s) = tmp(id_ext);
            end
            
            % extract all azimuth and elevation lines in a matrix with 2 layers -
            % 1st is azimuth, 2nd is elevation
            az_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-az') - 6;
            el_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-el') - 6;
            num_sat = numel(az_start);
            this.sat.az = [this.sat.az; zeros(n_epo, num_sat)];
            this.sat.el = [this.sat.el; zeros(n_epo, num_sat)];
            for s = 1 : num_sat
                az = sscanf(txt(bsxfun(@plus, repmat(az_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))', '%f');;
                el = sscanf(txt(bsxfun(@plus, repmat(el_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))', '%f');;
                this.sat.az(id_int, s) = az(id_ext);
                this.sat.el(id_int, s) = el(id_ext);
            end
        end
    end
    
     %% METHODS PLOTTING FUNCTIONS
    % ==================================================================================================================================================
    
    % Various debug images
    % name variant:
    %   c cartesian
    %   s scatter
    %   p polar
    %   m mixed
    methods (Access = public)
        
        function showAll(this)
            this.toString
            this.showAll@Receiver_Commons();
        end
    end
end