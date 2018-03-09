%   CLASS File_Wizard
% =========================================================================
%
% DESCRIPTION
%   Class to check and prepare the files needed by the processing
%   (e.g. navigational files)
%
% EXAMPLE
%   settings = FTP_Server();
%
% FOR A LIST OF CONSTANTs and METHODS use doc File_Wizard
%
% REQUIRES:
%   goGPS settings;
%
% COMMENTS
%   Server structure:

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, Giulio Tagliaferro, ...
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

classdef File_Wizard < handle
    
    properties (SetAccess = protected, GetAccess = public)
        % custom ftp downloader
        ftpd_custom;
        % download sources
        source = struct('cddis', struct('ftpd', FTP_Downloader('ftp://cddis.gsfc.nasa.gov/','21'), ...
            'par',  struct('gps', struct('path', '/gnss/', ...
            'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
            'eph', struct('final', 'products/${WWWW}/igs${WWWWD}.sp3', ...
            'rapid', 'products/${WWWW}/igr${WWWWD}.sp3', ...
            'ultra', 'products/${WWWW}/igu${WWWWD}_${6H}.sp3', ...
            'broadcast',  'data/daily/${YYYY}/brdc/brdc${DOY}0.${YY}n'), ...
            'clk', struct('final', 'products/${WWWW}/igs${WWWWD}.clk', ...
            'rapid', 'products/${WWWW}/igr${WWWWD}.clk'), ...
            'clk_30s', struct('final', 'products/${WWWW}/igs${WWWWD}.clk_30s'), ...
            'erp', struct('final', 'products/${WWWW}/igs${WWWW}7.erp', ...
            'rapid', 'products/${WWWW}/igr${WWWWD}.erp', ...
            'ultra', 'products/${WWWW}/igu${WWWWD}_${6H}.erp')),...
            'emx', struct('name', 'Natural Resources Canada - NRCAN GNSS solution', ...
            'eph', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.clk'), ...
            'erp', struct('final', 'products/${WWWW}/igs${WWWW}7.erp')), ...
            'cod', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
            'eph', struct('final', 'products/${WWWW}/cod${WWWWD}.eph'), ...
            'erp', struct('final', 'products/${WWWW}/igs${WWWW}7.erp'), ...
            'clk_05s', struct('final', 'products/${WWWW}/cod${WWWWD}.clk_05s')), ...
            'jpl', struct('name', 'Jet Propulsion Laboratory (JPL) - NASA (USA)', ...
            'eph', struct('final', 'products/${WWWW}/jpl${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'products/${WWWW}/jpl${WWWWD}.clk'), ...
            'erp', struct('final', 'products/${WWWW}/jpl${WWWW}7.erp')), ...
            'esa', struct('name', 'European Space Operations Centre (ESA) - ESA (Europe)', ...
            'eph', struct('final', 'products/${WWWW}/esa${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'products/${WWWW}/esa${WWWWD}.clk'), ...
            'erp', struct('final', 'products/${WWWW}/esa${WWWW}7.erp')), ...
            'gbm', struct('name', 'GeoForschungsZentrum Potsdam (GBM) - GFZ (Germany)', ...
            'eph', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
            'broadcast',  'data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
            'clk_30s', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
            'erp', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.erp'), ...
            'bds', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.bias')) ...
            )), ...
            'glo', struct('path', '/', ...
            'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
            'eph', struct('final', 'glonass/products/${WWWW}/igl${WWWWD}.sp3', ...
            'rapid', 'glonass/products/${WWWW}/igv${WWWWD}.sp3', ...
            'broadcast',  'pub/gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
            'erp', struct('final', 'glonass/products/${WWWW}/igs${WWWW}7.erp', ...
            'rapid', 'glonass/products/${WWWW}/igr${WWWWD}.erp', ...
            'ultra', 'glonass/products/${WWWW}/igu${WWWWD}_${6H}.erp')), ...
            'emx', struct('name', 'Natural Resources Canada - NRCAN GNSS solution', ...
            'eph', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.clk'), ...
            'erp', struct('final', 'products/${WWWW}/igs${WWWW}7.erp')), ...
            'cod', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
            'erp', struct('final', 'products/${WWWW}/cod${WWWWD}.eph'), ...
            'eph', struct('final', 'products/${WWWW}/cod${WWWWD}.eph'), ...
            'clk_05s', struct('final', 'products/${WWWW}/cod${WWWWD}.clk_05s')), ...
            'esa', struct('name', 'European Space Operations Centre (ESA) - ESA (Europe)', ...
            'eph', struct('final', 'products/${WWWW}/esa${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'products/${WWWW}/esa${WWWWD}.clk'), ...
            'erp', struct('final', 'products/${WWWW}/esa${WWWW}7.erp')), ...
            'gbm', struct('name', 'GeoForschungsZentrum Potsdam (GFZ)', ...
            'eph', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
            'broadcast',  '/pub/gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
            'clk_30s', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
            'erp', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWW}7.erp'), ...
            'bds', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.bias')))), ...
            'mxd', struct('path', '/', ...
            'center', struct('gbm', struct('name', 'GeoForschungsZentrum Potsdam (GFZ)', ...
            'eph', struct('final', 'gnss/products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
            'broadcast', 'gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
            'clk_30s', struct('final', 'gnss/products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
            'erp', struct('final', 'gnss/products/mgex/${WWWW}/gbm${WWWWD}.erp'), ...
            'bds', struct('final', 'gnss/products/mgex/${WWWW}/gbm${WWWWD}.bias')),...
            'emx', struct('name', 'Natural Resources Canada - NRCAN GNSS solution', ...
            'eph', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.sp3'), ...
            'clk_30s', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.clk'), ...
            'erp', struct('final', 'gnss/products/${WWWW}/igs${WWWW}7.erp')), ...
            'cod', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
            'erp', struct('final', 'gnss/products/${WWWW}/cod${WWWWD}.eph'), ...
            'eph', struct('final', 'gnss/products/${WWWW}/cod${WWWWD}.eph'), ...
            'clk_05s', struct('final', 'gnss/products/${WWWW}/cod${WWWWD}.clk_05s')), ...
            'codexp', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
            'erp', struct('final', 'gps/products/mgex/${WWWW}/COD0MGXFIN_${YYYY}${DOY}0000_03D_12H_ERP.ERP'), ...
            'eph', struct('final', 'gps/products/mgex/${WWWW}/COD0MGXFIN_${YYYY}${DOY}0000_01D_05M_ORB.SP3'), ...
            'bds', struct('final', 'gps/products/mgex/${WWWW}/COD0MGXFIN_${YYYY}${DOY}0000_01D_01D_OSB.BIA'), ...
            'clk_30s', struct('final', 'gps/products/mgex/${WWWW}/COD0MGXFIN_${YYYY}${DOY}0000_01D_30S_CLK.CLK')), ...
            'tum', struct('name', 'Technische Universitaet Muenchen', ...
            'eph', struct('final', 'gnss/products/mgex/${WWWW}/tum${WWWWD}.sp3', ...
            'broadcast',  'gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
            'clk', struct('final', 'gnss/products/mgex/${WWWW}/tum${WWWWD}.sp3'), ...
            'erp', struct('final', 'gnss/products/mgex/${WWWW}/tum${WWWWD}.erp'), ...
            'bds', struct('final', 'gnss/products/mgex/${WWWW}/tum${WWWWD}.bias'))...
            ))))... %);
            ,'ign',struct('ftpd', FTP_Downloader('ftp://igs.ign.fr/','21'), ...
            'par',  struct('mxd', struct('path', '/', ...
            'center', struct('cas', struct('name', 'Chinese Academy of Sciences (CAS) ', ...
            'dcb', struct('final', 'pub/igs/products/mgex/dcb/${YYYY}/CAS0MGXRAP_${YYYY}${DOY}0000_01D_01D_DCB.BSX.gz')))))) );
        
        date_start; % first epoch of common observations (among all the obs files)
        date_stop;  % last epoch of common observations (among all the obs files)
        center_name;% center code
        rm;         % resource manager
        sys_c;      % system collector
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        state = Global_Configuration.getCurrentSettings();     %  Global state, to import custom server and service preferences
    end
    
    properties (SetAccess = private, GetAccess = private)
        log = Logger.getInstance(); % Handler to the log object
        fnp = File_Name_Processor();
        ftp_downloaders;
    end
    
    methods
        function this = File_Wizard(state)
            % Constructor
            %  SYNTAX File_Wizard(<state>)
            % Uses state for getting settings info
            % Modify state to update eph_name and clk_name
            
            if (nargin >= 1)
                this.state = handle(state);
            else
                this.state = Global_Configuration.getCurrentSettings();
            end
            this.rm = Remote_Resource_Manager(this.state.getRemoteSourceFile);
            this.sys_c = this.state.cc.SYS_C(this.state.cc.active_list);
        end
        
        function [status] = conjureResource(this, resource_name, date_start, date_stop, center_name)
            if nargin < 3
                date_start = this.date_start;
                date_stop = this.date_stop;
            else
                this.date_start = date_start;
                this.date_stop = date_stop;
            end
            if nargin < 5
                center_name  = this.center_name;
            end
            [file_tree, latency] = this.rm.getFileStr(center_name, resource_name);
            if isempty(file_tree)
                status = false;
                return;
            end
            if isempty(latency)
                latency = [-Inf +Inf];
            end
            n_h_passed = (GPS_Time.now() - date_stop)/3600;
            % check local
            this.log.addMessage(this.log.indent('Checking local folders ...\n'))
            [status, file_tree] = this.navigateTree(file_tree, 'local_check');
            if status
                this.log.addMessage(this.log.indent('All files has been found locally\n'))
            else
                this.log.addMessage(this.log.indent('Some files not found locally\n'))
            end
            % check remote
            if  this.state.flag_check_remote && not(status)
                if n_h_passed  < latency(1)
                    status = false;
                else
                    this.log.addMessage(this.log.indent('Checking remote folders ...\n'))
                    [status, file_tree] = this.navigateTree(file_tree, 'remote_check');
                    if status
                        this.log.addMessage(this.log.indent('All files has been found remotely\n'));
                    else
                        this.log.addMessage(this.log.indent('Some files not found remotely\n'))
                    end
                end
                if status
                    this.log.addMessage(this.log.indent('Downloading Resources ...\n'));
                    [status, ~] = this.navigateTree(file_tree, 'download');
                    if not(status)
                        this.log.addWarning('Not all file have been found or uncopress\n');
                    end
                end
            end
            this.sys_c = this.state.cc.SYS_C(this.state.cc.active_list); % set sys_c again from constellation collector
        end
        
        function idx = getServerIdx(this, address , port)
            % get idx of server if not present open the connection
            if nargin < 3
                port = 21;
            end
            idx = 0;
            for i = 1 : length(this.ftp_downloaders)
                if strcmp(this.ftp_downloaders{i}.getAddress , address)
                    idx = i;
                    return
                end
            end
            if idx == 0
                this.ftp_downloaders{end+1} = FTP_Downloader(address, port);
                idx = length(this.ftp_downloaders);
            end
        end
        
        function [status, file_tree] = navigateTree(this, file_tree, mode)
            % navigate into file tree and perform operations
            % INPUT:
            %     file_tree: structure containing the file tree and the
            %     logical operators
            %     mode: 'local_check' , 'remote_check' , 'download'
            status = false;
            if iscell(file_tree) % if is a leaf
                if strcmp(file_tree{1}, 'null') || file_tree{2}
                    status = true;
                    if strcmp(mode, 'download') && ~strcmp(file_tree{1},'null') && file_tree{3} ~=0
                        loc_n = file_tree{3};
                        f_struct = this.rm.getFileLoc(file_tree{1});
                        f_name = f_struct.filename;
                        f_path = [f_struct.(['loc' sprintf('%03d',loc_n)]) f_name];
                        step_s = this.fnp.getStepSec(f_path);
                        dsa = this.date_start.getCopy();
                        dso = this.date_stop.getCopy();
                        dsa.addIntSeconds(-step_s);
                        dso.addIntSeconds(+step_s);
                        file_name_lst = this.fnp.dateKeyRepBatch(f_path, dsa, dso);
                        status = true;
                        f_status_lst = file_tree{4};
                        for i = 1 : length(file_name_lst)
                            if ~f_status_lst(i)
                                file_name = file_name_lst{i};
                                [server] = regexp(file_name,'(?<=\?{)\w*(?=})','match'); % saerch for ?{server_name} in paths
                                server = server{1};
                                file_name = strrep(file_name,['?{' server '}'],'');
                                [s_ip, port] = this.rm.getServerIp(server);
                                idx = this.getServerIdx(s_ip, port);
                                out_dir = this.state.getFileDir(file_name);
                                out_dir =  this.fnp.dateKeyRepBatch(out_dir, dsa, dso);
                                out_dir = out_dir{1};
                                status = status && this.ftp_downloaders{idx}.downloadUncompress(file_name, out_dir);
                            end
                        end
                    end
                elseif ~strcmp(mode, 'download')
                    f_struct = this.rm.getFileLoc(file_tree{1});
                    if isempty(f_struct.filename)
                        this.log.addError(sprintf('File resource "%s" not found: remote_resource.ini seems to be corrupted', file_tree{1}));
                    else
                        f_name = f_struct.filename;
                        this.state.setFile(f_name);
                        if strcmp(mode, 'local_check')
                            f_path = this.fnp.checkPath([this.state.getFileDir(f_name) filesep f_name]);
                            step_s = this.fnp.getStepSec(f_path);
                            dsa = this.date_start.getCopy();
                            dso = this.date_stop.getCopy();
                            dsa.addIntSeconds(-step_s);
                            dso.addIntSeconds(+step_s);
                            file_name_lst = this.fnp.dateKeyRepBatch(f_path, dsa, dso);
                            status = true;
                            f_status_lst = false(length(file_name_lst),1); %file list to be saved in tree with flag of downloaded or not
                            for i = 1 : length(file_name_lst)
                                f_status = exist(file_name_lst{i}, 'file') == 2;
                                f_status_lst(i) = f_status;
                                status = status && f_status;
                                if status
                                    this.log.addStatusOk(sprintf('%s have been found locally',this.fnp.getFileName(file_name_lst{i})));
                                end
                            end
                            if status
                                file_tree{3} = 0;
                            end
                            file_tree{4} = f_status_lst;
                        elseif strcmp(mode, 'remote_check')
                            
                            for i = 1 : f_struct.loc_number
                                f_path = [f_struct.(['loc' sprintf('%03d',i)]) f_name];
                                step_s = this.fnp.getStepSec(f_path);
                                dsa = this.date_start.getCopy();
                                dso = this.date_stop.getCopy();
                                dsa.addIntSeconds(-step_s);
                                dso.addIntSeconds(+step_s);
                                file_name_lst = this.fnp.dateKeyRepBatch(f_path, dsa, dso);
                                status = true;
                                f_status_lst = file_tree{4};
                                for j = 1 : length(file_name_lst)
                                    if ~f_status_lst(j)
                                        file_name = file_name_lst{j};
                                        [server] = regexp(file_name,'(?<=\?{)\w*(?=})','match'); % saerch for ?{server_name} in paths
                                        server = server{1};
                                        file_name = strrep(file_name,['?{' server '}'],'');
                                        [s_ip, port] = this.rm.getServerIp(server);
                                        idx = this.getServerIdx(s_ip, port);
                                        status = status && this.ftp_downloaders{idx}.check(file_name);
                                        if status
                                            this.log.addStatusOk(sprintf('%s have been found remotely', this.fnp.getFileName(file_name)));
                                        else
                                            this.log.addStatusOk(sprintf('%s have not been found remotely', this.fnp.getFileName(file_name)));
                                            break
                                        end
                                    end
                                end
                                if status
                                    file_tree{3} = i;
                                    break
                                end
                            end
                        end
                    end
                end
                file_tree{2} = status;
            else % is if a brach go deeper into the branches / leafs
                b_name = fieldnames(file_tree);
                b_name = b_name{1};
                or_flag = strcmp(b_name, 'or');
                and_flag = strcmp(b_name, 'and');
                if or_flag
                    status = false;
                elseif and_flag
                    status = true;
                end
                branch = fieldnames(file_tree.(b_name));
                for i = 1 :length(branch)
                    [status_b, file_tree_b] = this.navigateTree(file_tree.(b_name).(branch{i}), mode);
                    if or_flag
                        status  = status || status_b;
                        if status
                            file_tree.(b_name).(branch{i}) = file_tree_b;
                            break
                        end
                    elseif and_flag
                        status  = status && status_b;
                    end
                    file_tree.(b_name).(branch{i}) = file_tree_b;
                end
                
            end
        end
        
        function conjureFiles(this, date_start, date_stop, center_name)
            % Prepare all the files needed for processing
            if (nargin == 1)
                [date_start, date_stop] = this.conjureObsFile();
            end
            dsa = date_start.getCopy();
            dso = date_stop.getCopy();
            
            if nargin < 3
               center_name = this.state.getRemoteCenter();
            end
            % check if selected center os compatible with selected
            % constellation
            centers = this.rm.getData('CENTER','available');
            is_ok = false;
            for i = 1 : length(centers)
                if ~is_ok
                    split = strsplit(centers{i},'@');
                    %centername = split{2};
                    %centername = strsplit(centerconst,'_');
                    %centername = centername{1};
                    sys_c = split{1};
                    if length(intersect(this.sys_c,sys_c)) == length(this.sys_c)
                        if strcmp(split{2}, center_name{1})
                            is_ok = true;
                            this.center_name = center_name{1};
                        end
                    end
                end
            end
            
            if ~is_ok
                this.log.addError(['Selected center: ' this.state.preferred_center{1} ' not compatible with selected constellations: ' this.sys_c]);
            end
                
            % Prepare all the files needed for processing
            
            this.state.setProcessingTime(dsa, dso, false);
            this.state.updateNavFileName();
            this.state.updateErpFileName();
            this.conjureNavFiles(dsa, dso);
            
            this.conjureDCBFiles(dsa, dso);
            this.conjureCRXFiles(dsa, dso);
            this.conjureIonoFiles(dsa, dso);
        end
        
        function [first_epoch, last_epoch] = conjureObsFile(this)
            % Prepare the extended file name of the files to be used in goGPS
            % In a future here I'll download the required navigational files of a station in a network
            
            first_rec_files = this.state.getRecPath(1);
            fh = File_Rinex(first_rec_files, 100);
            this.log.newLine();
            first_epoch = fh.first_epoch.first;
            last_epoch = fh.last_epoch.last;
        end
        
        %{
        function conjureDCBFiles(this, date_start, date_stop)
            status = this.conjureResource('dcb',date_start, date_stop, 'cas');
            if status
                this.log.addMarkedMessage('All DCB files present')
            else
                this.log.addMarkedMessage('Not all DCB files found program might misbehave')
            end
        end
        %}
        
        function conjureDCBFiles(this, date_start, date_stop)
            
            % SYNTAX:
            %   this.conjureDCBFiles(gps_week, gps_time);
            %
            % INPUT:
            %   date_start = starting GPS_Time
            %   date_stop = ending GPS_Time
            %
            % OUTPUT:
            %
            % DESCRIPTION:
            %   Download of CAS .DCB files from the IGN server.
            if date_start.getCalEpoch >= 2013 % use CAS DCB
                dcb_ok = true;
                % check if file are present
                fnp = File_Name_Processor();
                ss = 'mxd';
                archive = 'ign';
                provider = 'cas';
                dcb_type = 'final';
                dcb_name = this.source.(archive).par.(ss).center.(provider).dcb.(dcb_type);
                tmp_date_start = date_start.getCopy;
                tmp_date_stop = date_stop.getCopy;
                file_list = fnp.dateKeyRepBatch(dcb_name, tmp_date_start, tmp_date_stop);
                names = {};
                for i = 1 : length(file_list)
                    [~, name, ext] = fileparts(file_list{i});
                    names{end+1} = name;
                    if exist(checkPath([this.state.getDcbDir '/' name]), 'file') ~= 2
                        dcb_ok = false;
                    end
                end
                this.state.setDcbFile(names);
                if (~dcb_ok)
                    this.source.(archive).ftpd.download(this.source.(archive).par.(ss).path, file_list, this.state.getDcbDir());
                    for i = 1 : length(file_list)
                        [~, name, ext] = fileparts(file_list{i});
                        if (isunix())
                            system(['gzip -fd ' this.state.getDcbDir() '/' name ext]);
                        else
                            try
                                [status, result] = system(['".\utility\thirdParty\7z1602-extra\7za.exe" -y x ' '"' this.state.getDcbDir() '/' name ext '"' ' -o' '"' this.state.getDcbDir() '"']); %#ok<ASGLU>
                                delete([this.state.getDcbDir() '/' name ext]);
                            catch
                                this.log.addWarning(sprintf(['Please decompress the ' name ext ' file before trying to use it in goGPS.\n']));
                                compressed = 1;
                            end
                        end
                    end
                    
                    
                else
                    this.log.addStatusOk('Dcb files are present ^_^');
                    this.log.newLine();
                end
                
                
                
                
            else % use DCB from CODE
                gps_week = double([date_start.getGpsWeek; date_stop.getGpsWeek ]);
                gps_time = [date_start.getGpsTime; date_stop.getGpsTime ];
                %[file_dcb, compressed] = download_dcb(gps_weeks, gps_times);
                
                % Pointer to the global settings:
                state = Global_Configuration.getCurrentSettings();
                
                file_dcb = {};
                compressed = 0;
                
                %AIUB FTP server IP address
                % aiub_ip = '130.92.9.78'; % ftp.aiub.unibe.ch
                aiub_ip = 'ftp.aiub.unibe.ch';
                
                %download directory
                down_dir = state.dcb_dir;
                
                %convert GPS time to time-of-week
                gps_tow = weektime2tow(gps_week, gps_time);
                
                % starting time
                date_f = gps2date(gps_week(1), gps_tow(1));
                
                % ending time
                date_l = gps2date(gps_week(end), gps_tow(end));
                
                % Check / create output folder
                if not(exist(down_dir, 'dir'))
                    mkdir(down_dir);
                end
                
                log.addMessage(log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                
                year_orig  = date_f(1) : 1 : date_l(1);
                if (length(year_orig) < 1)
                    %fprintf('ERROR: Data range not valid.\n')
                    return
                elseif (length(year_orig) == 1)
                    month = date_f(2) : 1 : date_l(2);
                    year = year_orig;
                else
                    month = date_f(2) : 1 : 12;
                    year  = date_f(1).*ones(size(month));
                    for y = 2 : length(year_orig)-1
                        month = [month 1 : 1 : 12];
                        year = [year (year+y-1).*ones(1,12)];
                    end
                    month = [month 1 : 1 : date_l(2)];
                    year  = [year date_l(1).*ones(1,date_l(2))];
                end
                
                %connect to the DCB server
                try
                    ftp_server = ftp(aiub_ip);
                catch
                    fprintf(' connection failed.\n');
                    this.state.setDcbFile({''});
                    return
                end
                
                fprintf('\n');
                
                m = 0;
                
                for y = 1 : length(year_orig)
                    
                    %target directory
                    s = ['/CODE/', num2str(year_orig(y))];
                    
                    cd(ftp_server, '/');
                    cd(ftp_server, s);
                    
                    while(m <= length(month)-1)
                        
                        m = m + 1;
                        
                        ff = {'P1C1','P1P2'};
                        
                        for p = 1 : length(ff)
                            %target file
                            s2 = [ff{p} num2str(two_digit_year(year(m)),'%02d') num2str(month(m),'%02d') '.DCB.Z'];
                            if not(exist([down_dir '/' s2(1:end-2)]) == 2)
                                try
                                    mget(ftp_server,s2,down_dir);
                                    if (isunix())
                                        system(['uncompress -f ' down_dir '/' s2]);
                                    else
                                        try
                                            [status, result] = system(['".\utility\thirdParty\7z1602-extra\7za.exe" -y x ' '"' down_dir '/' s2 '"' ' -o' '"' down_dir '"']); %#ok<ASGLU>
                                            delete([down_dir '/' s2]);
                                            s2 = s2(1:end-2);
                                        catch
                                            fprintf(['Please decompress the ' s2 ' file before trying to use it in goGPS.\n']);
                                            compressed = 1;
                                        end
                                    end
                                    fprintf(['Downloaded DCB file: ' s2 '\n']);
                                catch
                                    cd(ftp_server, '..');
                                    s1 = [ff{p} '.DCB'];
                                    mget(ftp_server,s1,down_dir);
                                    cd(ftp_server, num2str(year_orig(y)));
                                    s2 = [s2(1:end-2) '_TMP'];
                                    movefile([down_dir '/' s1], [down_dir '/' s2]);
                                    fprintf(['Downloaded DCB file: ' s1 ' --> renamed to: ' s2 '\n']);
                                end
                            else
                                fprintf([s2(1:end-2) ' already present\n']);
                            end
                            %cell array with the paths to the downloaded files
                            entry = {[down_dir, '/', s2]};
                            file_dcb = [file_dcb; entry]; %#ok<AGROW>
                            this.state.setDcbFile(file_dcb);
                        end
                        
                        if (month(m) == 12)
                            break
                        end
                    end
                end
                
                close(ftp_server);
                
                fprintf('Download complete.\n')
            end
        end
        
        function conjureNavFiles(this, date_start, date_stop)
            list_preferred = this.state.getPreferredEph();
            for i = 1 : length(list_preferred)
                status = this.conjureResource(list_preferred{i}, date_start, date_stop);
                if status
                    break
                end
            end
            if status
                this.log.addMarkedMessage('All ephemerids files present')
            else
                this.log.addMarkedMessage('Not all ephemerids files found program might misbehave')
            end
        end
        
        function conjureErpFiles(this, date_start, date_stop)
            list_preferred = this.state.preferred_erp;
            for i = 1 : length(list_preferred)
                status = this.conjureResource(list_preferred{i},date_start, date_stop);
                if status
                    break
                end
            end
            if status
                this.log.addMarkedMessage('All erp files present')
            else
                this.log.addMarkedMessage('Not all erp files found program might misbehave')
            end
        end
        
        function conjureIonoFiles(this, date_start, date_stop)
            list_preferred = this.state.preferred_iono;
            for i = 1 : length(list_preferred)
                status = this.conjureResource(['iono_' list_preferred{i}], date_start, date_stop);
                if status
                    break
                end
            end
            if status
                this.log.addMarkedMessage('All iono files present')
            else
                this.log.addMarkedMessage('Not all iono files found program might misbehave')
            end
            
        end
        
        function conjureCRXFiles(this, date_start, date_stop)
            % SYNTAX:
            %   this.conjureCRXFiles(gps_week, gps_time);
            %
            % INPUT:
            %   date_start = starting GPS_Time
            %   date_stop = ending GPS_Time
            %
            % OUTPUT:
            %
            % DESCRIPTION:
            %   Download of .CRX files from the AIUB FTP server.
            date_start = date_start.getCopy();
            date_stop = date_stop.getCopy();
            gps_week = double([date_start.getGpsWeek; date_stop.getGpsWeek ]);
            gps_time = [date_start.getGpsTime; date_stop.getGpsTime ];
            %[file_crx] = download_crx(gps_weeks, gps_times);
            
            % Pointer to the global settings:
            state = Global_Configuration.getCurrentSettings();
            
            %AIUB FTP server IP address
            % aiub_ip = '130.92.9.78'; % ftp.aiub.unibe.ch
            aiub_ip = 'ftp.aiub.unibe.ch';
            
            %download directory
            down_dir = state.crx_dir;
            
            %convert GPS time to time-of-week
            [~, start_sow] = date_start.getGpsWeek;
            [~, stop_sow] = date_stop.getGpsWeek;
            gps_tow = [start_sow; stop_sow];
            
            % starting time
            date_f = gps2date(gps_week(1), gps_tow(1));
            
            % ending time
            date_l = gps2date(gps_week(end), gps_tow(end));
            
            % Check / create output folder
            if not(exist(down_dir, 'dir'))
                mkdir(down_dir);
            end
                        
            year  = date_f(1) : 1 : date_l(1);
            file_crx = cell(length(year),1);
            
            %connect to the CRX server
            try
                ftp_server = ftp(aiub_ip);
            catch
                this.log.addMessage(this.log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                this.log.addError('Connection failed.\n');
                return
            end
            
            fprintf('\n');
            
            for y = 1 : length(year)
                
                %target directory
                s = '/BSWUSER52/GEN';
                
                cd(ftp_server, '/');
                cd(ftp_server, s);
                
                %target file
                s2 = ['SAT_' num2str(year(y),'%04d') '.CRX'];
                
                % read the last modification of the CRX
                d = dir([down_dir, '/', s2]);
                t = GPS_Time(d.datenum);
                
                % If there's no CRX or the CRX is older than the day of the processing and it has not been downloaded in the last day
                % do not do
                if isempty(d) || ((t < date_stop.addSeconds(10*86400)) && (GPS_Time.now - t > 43200))
                    this.log.addMessage(this.log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                    %if not(exist([down_dir, '/', s2]) == 2)
                    mget(ftp_server,s2,down_dir);
                    this.log.addMessage(this.log.indent(sprintf(['Downloaded CRX file: ' s2 '\n'])));
                end
                
                % cell array with the paths to the downloaded files
                file_crx{y} = [down_dir, '/', s2];                                
            end
            
            try
                close(ftp_server);
            catch
            end
            
            this.log.addStatusOk(this.log.indent(('CRX ready!\n')))
            
        end
        
    end
    
end
