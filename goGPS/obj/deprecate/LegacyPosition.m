classdef LegacyPosition < Exportable
        
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
    %  Written by:       Alice Bonfiglio
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
    
    properties (Access = private)
        log     % logger
    end

    properties
        time     % time as GPS_Time                    GPS_Time [1 x 1] stores n_epoch
        
        xyz      % point position
        enu      % point position
        
        lat      % ellipsoidal latitude
        lon      % ellipsoidal longitude
        h_ellips % ellipsoidal height
        h_ortho  % hortometric height
        
        n_sat
        hdop
        khdop
        a_fix
        s_rate
    end
        
    methods
        % Creator
        function this = LegacyPosition(file_name)
            % Core object creator
            this.log = Logger.getInstance();
            this.reset();
            if nargin == 1
                if iscell(file_name)
                    for f = 1 : length(file_name)
                        this.log.addMessage(sprintf('Importing %s', file_name{f}));
                        if exist(file_name{f}, 'file')
                            this.appendPosition(file_name{f});
                        else
                            this.log.addMessage(sprintf('Error loading the last file, the file does not exists'));
                        end
                    end
                else
                    this.log.addMessage(sprintf('Importing %s', file_name));
                    if exist(file_name, 'file')
                        this.appendPosition(file_name);
                    else
                        this.log.addMessage(sprintf('Error loading the file, it does not exists'));
                    end
                end
            end
        end
    end
    
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        function reset(this)
            this.time = GPS_Time();
            this.xyz = [];
            this.enu = [];
            this.lat = [];
            this.lon = [];
            
            this.h_ellips = [];
            this.h_ortho = [];
            
            this.n_sat = [];
            this.hdop =  [];
            this.khdop = [];
            this.a_fix = [];
            this.s_rate = [];
        end
        
        function importPosition(this, file)
            % Import after reset a position file
            % SYNTAX: this.importPosition(file)
            this.reset();
            this.appendPosition(file)
        end
        
        function appendPosition (this, file)
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
            ko_lines = find(lim(:, 3) ~= median(lim(:,3)));
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
    end
    
    
    methods (Static)
        function pos = loadBatch(file_name, run_start, run_stop)
            % Load all the position files of a session
            %
            % SYNTAX:  
            %   position = Position.loadBatch(file_name, run_start, run_stop)
            %
            % INPUT:
            %   file_name     it should include the key ${RUN} that will be substituted with a 3 digits number containing the run, from run_start to run_stop
            %   run_start     number of the first run to load
            %   run_stop      number of the last run to load
            %
            % OUTPUT:
            %   position      position object containing all the data loaded from the files
            %
            GPS_RUN = '${RUN}';
            file_name_list = {};
            r = 0;
            for run = run_start : run_stop
                r = r + 1;
                file_name_list{r} = strrep(file_name, GPS_RUN, sprintf('%03d', run));
            end
            if numel(file_name_list) > 0
                pos = Position(file_name_list);
            else
                pos = Position();
            end
        end
    end
end
