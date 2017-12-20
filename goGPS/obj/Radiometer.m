classdef Radiometer < handle
        
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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
        time          % time as GPS_Time                                    GPS_Time [1 x 1] stores n_epoch
        
        zwd           % equivalent wet zenith delay [cm]                    double   [n_epoch x 1]
        sigma_zwd     % zwd sigma error [cm]                                double   [n_epoch x 1]
        zwd_21        % equivalent wet zenith delay - 21 GHz only [cm]      double   [n_epoch x 1]
        sigma_zwd_21  % zwd sigma error [cm]                                double   [n_epoch x 1]
        zwd_saast     % zenital wet delay from Saastamoinen model           double   [n_epoch x 1]
        air_mass       % airmass of the observation in the zenith            double   [n_epoch x 1]
        el            % elevation of the wvr observation [deg]              double   [n_epoch x 1]
        az            % azimuth of the wvr observation [deg]                double   [n_epoch x 1]
        
        bright_21     % observed brightness at 21.0 GHz [K]                 double   [n_epoch x 1]
        sigma_br_21   % sigma error of brightness at 21 GHz [K]             double   [n_epoch x 1]
        bright_314    % observed brightness at 31.4 GHz [K]                 double   [n_epoch x 1]
        sigma_br_314  % sigma error of brightness at 31.4 GHz [K]           double   [n_epoch x 1]
        zbright_21    % equivalent zenith brightness at 21.0 GHz [K]        double   [n_epoch x 1]
        sigma_zbr_21  % sigma error of eq zenith brightness at 21 GHz [K]   double   [n_epoch x 1]
        zbright_314   % equivalent zenith brightness at 31.4 GHz [K]        double   [n_epoch x 1]
        sigma_zbr_314 % sigma error of eq zenith brightness at 31.4 GHz [K] double   [n_epoch x 1]
        lwc           % liquid water content in zenith direction [mm]       double   [n_epoch x 1]
        sigma_lwc     % sigma error of liquid water content [mm]            double   [n_epoch x 1]
        
    end
    
    properties (Access = private)
        log     % logger
    end
    
    methods
        % Creator
        function this = Radiometer(file_name)
            % Core object creator
            this.log = Logger.getInstance();
            this.reset();
            if nargin == 1
                if iscell(file_name)
                    for f = 1 : length(file_name)
                        this.log.addMessage(sprintf('Importing %s', file_name{f}));
                        if exist(file_name{f}, 'file')
                            this.appendRadiometer(file_name{f});
                        else
                            this.log.addMessage(sprintf('Error loading the last file, the file does not exists'));
                        end
                    end
                else
                    this.log.addMessage(sprintf('Importing %s', file_name));
                    if exist(file_name, 'file')
                        this.appendRadiometer(file_name);
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
            this.zwd  = [];
            this.sigma_zwd  = [];
            this.zwd_21  = [];
            this.sigma_zwd_21  = [];
            this.zwd_saast  = [];
            this.air_mass  = [];
            this.el  = [];
            this.az  = [];
            
            this.bright_21  = [];    
            this.sigma_br_21  = [];
            this.bright_314  = [];
            this.sigma_br_314  = [];
            this.zbright_21  = [];
            this.sigma_zbr_21  = [];
            this.zbright_314  = [];
            this.sigma_zbr_314  = [];
            this.lwc  = [];
            this.sigma_lwc  = [];
        end
        
        function importRadiometer(this, file)
            % Import after reset a tropo file
            % SYNTAX: this.importTropo(file)
            this.reset();
            this.appendRadiometer(file)
        end
        
        function appendRadiometer (this, file)
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
            eoh = 3; % end of header (the header of radiometer files contains 3 lines)
            
            % corrupted lines
            ko_lines = find(lim(4:end, 3) ~= median(lim(4:end,3)));
            for l = numel(ko_lines) : -1 : 1
                txt(lim(ko_lines(l)+3, 1) : lim(ko_lines(l)+3, 2) + 1) = [];
            end
            
            % extract all the epoch lines
            data = sscanf(txt(lim(4,1):end)','%4d %2d %2d %2d %2d %2d %7f %6f %7f %6f %7f %6f %6f %6f %7f %6f %7f %6f %7f %6f %7f %6f %7f %6f\n');
            data = reshape(data, 24, numel(data)/24)';
            
            % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(data(:,1:6));
            else
                this.time.append6ColDate(data(:,1:6));
            end
                        
            zwd = data(:,7);
            sigma_zwd = data(:,8);
            zwd_21 = data(:,9);
            sigma_zwd_21 = data(:,10);
            zwd_saast = data(:,11);
            air_mass = data(:,12);
            el = data(:,13);
            az = data(:,14);
            az = az + 180 * ((90-el) < 0);
            for i = 1 : length(az)
                if az(i) > 360
                az(i) = mod(az(i),360);
                end
            end
            el = 90 - abs(90-el);
            
            bright_21  = data(:,15);    
            sigma_br_21 = data(:,16);
            bright_314 = data(:,17);
            sigma_br_314 = data(:,18);
            zbright_21 = data(:,19);
            sigma_zbr_21 = data(:,20);
            zbright_314 = data(:,21);
            sigma_zbr_314 = data(:,22);
            lwc = data(:,23);
            sigma_lwc = data(:,24);
            
           % Append in obj
            this.zwd = [this.zwd; zwd];
            this.sigma_zwd = [this.sigma_zwd; sigma_zwd];
            this.zwd_21 = [this.zwd_21; zwd_21];
            this.sigma_zwd_21 = [this.sigma_zwd_21; sigma_zwd_21];
            this.zwd_saast = [this.zwd_saast; zwd_saast];
            this.air_mass = [this.air_mass; air_mass];
            this.el = [this.el; el];
            this.az = [this.az; az];
            this.bright_21 = [this.bright_21; bright_21];
            this.sigma_br_21 = [this.sigma_br_21; sigma_br_21];
            this.bright_314 = [this.bright_314; bright_314];
            this.sigma_br_314 = [this.sigma_br_314; sigma_br_314];
            
            this.zbright_21 = [this.zbright_21; zbright_21];
            this.sigma_zbr_21 = [this.sigma_zbr_21; sigma_zbr_21];
            this.zbright_314 = [this.zbright_314; zbright_314];
            this.sigma_zbr_314 = [this.sigma_zbr_314; sigma_zbr_314];
            this.lwc = [this.lwc; lwc];
            this.sigma_lwc = [this.sigma_lwc; sigma_lwc];
            
            clear data;
        end
    end   
    
    methods (Static)
        function wvr = loadBatch(file_name, run_start, run_stop)
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
            
            if numel(file_name_list) > 0
                wvr = Radiometer(file_name_list);
            else
                wvr = Radiometer();
            end
        end
    end
end
