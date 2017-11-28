classdef Tropo < handle
        
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
        time     % time as GPS_Time                    GPS_Time [1 x 1] stores n_epoch
        
        zhd      % zenital hydrostatic delay           double   [n_epoch x 1]
        ztd      % zenital tropospheric delay          double   [n_epoch x 1]
        zwd      % zenital wet delay                   double   [n_epoch x 1]
        pwv      % precipitable water vapour           double   [n_epoch x 1]
        
        tgn      % tropospheric gradient north         double   [n_epoch x n_sat]
        tge      % tropospheric gradient east          double   [n_epoch x n_sat]

        sat = struct('id', [], ...    % id (e.g. G01)               char     [n_sat, 3]
                     'el', [], ...    % elevetion                   double   [n_epoch x n_sat]
                     'az', [], ...    % azimuth                     double   [n_epoch x n_sat]                     
                     'slant_td', []); % slant tropospheric delay    double   [n_epoch x n_sat]
    end
    
    properties (Access = private)
        log     % logger
    end
    
    methods
        % Creator
        function this = Tropo(file_name)
            % Core object creator
            this.log = Logger.getInstance();
            this.reset();
            if nargin == 1
                if iscell(file_name)
                    for f = 1 : length(file_name)
                        this.log.addMessage(sprintf('Importing %s', file_name{f}));
                        if exist(file_name{f}, 'file')
                            this.appendTropo(file_name{f});
                        else
                            this.log.addMessage(sprintf('Error loading the last file, the file does not exists'));
                        end
                    end
                else
                    this.log.addMessage(sprintf('Importing %s', file_name));
                    if exist(file_name, 'file')
                        this.appendTropo(file_name);
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
            this.zhd  = [];
            this.zwd  = [];
            this.ztd  = [];
            this.pwv  = [];
            
            this.tgn = [];
            this.tge = [];
            this.sat = struct('id', [], ...
                              'el', [], ...
                              'az', [], ...
                              'slant_td', []);
        end
        
        function importTropo(this, file)
            % Import after reset a tropo file
            % SYNTAX: this.importTropo(file)
            this.reset();
            this.appendTropo(file)
        end
        
        function appendTropo (this, file)
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
            n_epo = size(date, 1);
            
            % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(date, [], 1);
            else
                this.time.append6ColDate(date);
            end
            
            % extract all the ZHD lines
            string_zhd = txt(repmat(lim(2:end,1),1,17) + repmat(62:78, size(lim,1)-1, 1))';
            this.zhd = [this.zhd; sscanf(string_zhd,'%f')]; clear string_zhd
            
            % extract all the ZTD lines
            string_ztd = txt(repmat(lim(2:end,1),1,17) + repmat(78:94, size(lim,1)-1, 1))';
            this.ztd = [this.ztd; sscanf(string_ztd,'%f')]; clear string_ztd
            
            % extract all the TGN lines
            string_tgn = txt(repmat(lim(2:end,1),1,17) + repmat(95:111, size(lim,1)-1, 1))';
            this.tgn = [this.tgn; sscanf(string_tgn,'%f')]; clear string_tgn
            
            % extract all the TGE lines
            string_tge = txt(repmat(lim(2:end,1),1,17) + repmat(112:128, size(lim,1)-1, 1))';
            this.tge = [this.tge; sscanf(string_tge,'%f')]; clear string_tge
            
            % extract all the ZWD lines
            string_zwd = txt(repmat(lim(2:end,1),1,17) + repmat(129:145, size(lim,1)-1, 1))';
            this.zwd = [this.zwd; sscanf(string_zwd,'%f')]; clear string_zwd
            
            % extract all the PWV lines
            string_pwv = txt(repmat(lim(2:end,1),1,17) + repmat(146:162, size(lim,1)-1, 1))';
            this.pwv = [this.pwv; sscanf(string_pwv,'%f')]; clear string_pwv
            
            %  extract all STD values if present
            slant_start = regexp(txt(lim(1,1) : lim(1,2)),'STD') - 6;
            num_sat = numel(slant_start);
            this.sat.slant_td = [this.sat.slant_td; zeros(n_epo, num_sat)];
            for s = 1 : numel(slant_start)
                tmp = txt(bsxfun(@plus, repmat(slant_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1) - 1))';
                this.sat.slant_td((end - n_epo + 1) : end, s) = sscanf(tmp, '%f');
            end
            
            % extract all azimuth and elevation lines in a matrix with 2 layers -
            % 1st is azimuth, 2nd is elevation
            az_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-az') - 6;
            el_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-el') - 6;
            num_sat = numel(az_start);
            this.sat.az = [this.sat.az; zeros(n_epo, num_sat)];
            this.sat.el = [this.sat.el; zeros(n_epo, num_sat)];
            for s = 1 : num_sat
                az = txt(bsxfun(@plus, repmat(az_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))';
                el = txt(bsxfun(@plus, repmat(el_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))';
                this.sat.az((end - n_epo + 1) : end, s) = sscanf(az, '%f');
                this.sat.el((end - n_epo + 1) : end, s) = sscanf(el, '%f');
            end
        end
    end
    
    
    methods (Static)
        function tropo = loadBatch(file_name, run_start, run_stop)
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
                r = r + 1;
                file_name_list{r} = strrep(file_name, GPS_RUN, sprintf('%03d', run));
            end
            if numel(file_name_list) > 0
                tropo = Tropo(file_name_list);
            else
                tropo = Tropo();
            end
        end
    end
end
