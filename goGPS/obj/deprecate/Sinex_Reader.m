classdef Sinex_Reader < handle
        
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
    
    properties
        time          % time as matlab time                     float   [n_epoch x 1]
        ztd           % zenith total delay [mm]                 float   [n_epoch x 1]
        zwd           % zenith wet delay [mm]                   float   [n_epoch x 1]
        error         % formal error [mm]                       float   [n_epoch x 1]        
    end
    
    properties (Access = private)
        log     % logger
    end
    
    methods
        % Creator
        function this = Sinex_Reader(file_name, site)
            % Core object creator
            this.log = Core.getLogger();
            this.reset();
            if nargin == 1
                if iscell(file_name)
                    for f = 1 : length(file_name)
                        this.log.addMessage(sprintf('Importing %s', file_name{f}));
                        if exist(file_name{f}, 'file')
                            this.appendData(file_name{f}, site);
                        else
                            this.log.addMessage(sprintf('Error loading the last file, the file does not exists'));
                        end
                    end
                else
                    this.log.addMessage(sprintf('Importing %s', file_name));
                    if exist(file_name, 'file')
                        this.appendData(file_name, site);
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
            this.ztd  = [];
            this.zwd  = [];
            this.error  = [];
        end
        
        function importData(this, file, site)
            % Import after reset a tropo file
            % SYNTAX: this.importData(file)
            this.reset();
            this.appendSinex_Reader(file, site)
        end
        
        function appendData(this, file, site)
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
            eoh = 0; % end of header (the ONS station file has no header)
            
            % corrupted lines
            ko_lines = find(lim(1:end, 3) ~= median(lim(1:end,3)));
            for l = numel(ko_lines) : -1 : 1
                txt(lim(ko_lines(l), 1) : lim(ko_lines(l), 2) + 1) = [];
            end
            
            % extract all the epoch lines
            data = sscanf(txt(lim(1,1):end)','%f %f %f %f\n');
            data = reshape(data, 4, numel(data)/4)';
            
            julian_time = data(:,1);
            ztd = data(:,2) * 1e-3;
            zwd = data(:,3) * 1e-3;
            error = data(:,4) * 1e-3;
            
            time = jd2MatTime(julian_time);
            
           % Append in obj
           
           % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(time);
            else
                this.time.appendMatTime(time);
            end
            
            this.ztd = [this.ztd; ztd];
            this.zwd = [this.zwd; zwd];
            this.error = [this.error; error];
                       
            clear data;
        end
    end
    
    methods (Static)
        function ons = loadBatch(file_name, run_start, run_stop)
            % Load all the tropo files of a session
            %
            % SYNTAX:  
            %   tropo = Tropo.loadBatch(file_name, run_start, run_stop)
            %
            % INPUT:
            %   file_name     it should include the key ${RUN} that will be substituted with a 3 digits number containing the run, from run_start to run_stop
            %   run_start     number of the first run to load
            %   run_stop      number of the last run to load
            %
            % OUTPUT:
            %   tropo         troposphere object containing all the data loaded from the files
            %
            GPS_RUN = '${RUN}';
            file_name_list = {};
            r = 0;
            for run = run_start : run_stop
                cur_file_name = strrep(file_name, GPS_RUN, sprintf('%02d', run));
                if exist(cur_file_name, 'file') == 2
                    r = r + 1;
                    file_name_list{r} = cur_file_name;
                end
            end
            
            ons = Sinex_Reader(file_name_list);
        end
    end
end
