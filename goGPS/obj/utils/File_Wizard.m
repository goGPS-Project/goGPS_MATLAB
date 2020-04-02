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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, Giulio Tagliaferro, ...
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
        current_resource; % current resource being processed 
        
        nrt = false;%near real time flag % continue even if all the files have not been found
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        state = Core.getCurrentSettings();     %  Global state, to import custom server and service preferences
    end
    
    properties (SetAccess = private, GetAccess = private)
        log = Core.getLogger(); % Handler to the log object
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
                this.state = Core.getState();
            end
            this.rm = Remote_Resource_Manager(this.state.getRemoteSourceFile());
            this.sys_c = this.state.cc.SYS_C(this.state.cc.active_list);
            this.center_name = this.state.selected_center{1};
        end
        
        function [status] = conjureResource(this, resource_name, date_start, date_stop, center_name)
            % Conjure the deisdered resource giveng the desidered center
            % and the times bounds
            %
            % SYNTAX: 
            %   [status] = this.conjureResource(resource_name, date_start, date_stop, center_name)
            %
            % INPUT:
            %      resource_name = name of the resource (e.g. final_erp)
            %      date_start = GPS_Time start
            %      date_stop = GPS_Time stop
            %      center_name = name of the center (e.g. code)
            %
            % OUPUT:
            %     status = 1 everything has been found 0 no
            
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
            this.current_resource = resource_name;
            [file_tree, latency] = this.rm.getFileStr(center_name, resource_name);
            if isempty(file_tree)
                [file_tree, latency] = this.rm.getFileStr('default', resource_name);
                if isempty(file_tree)
                    status = false;
                    return;
                end
            end
            if isempty(latency)
                latency = [-Inf +Inf];
            end
            n_h_passed = (GPS_Time.now() - date_stop)/3600;
            % check local
            this.log.addMessage(this.log.indent('Checking local folders ...'));
            [status, file_tree] = this.navigateTree(file_tree, 'local_check');
            if status
                this.log.addStatusOk('All files have been found locally', 10);
            else
                this.log.addWarning('Some files not found locally');
            end
            % check remote
            if  this.state.isAutomaticDownload && ~status
                if latency(1)~=0 && n_h_passed  < latency(1)
                    this.log.addWarning(this.log.indent(sprintf('Not enough latency for finding all the %s orbits...', resource_name)));
                    status = false;
                else
                    this.log.addMessage(this.log.indent('Checking remote folders ...'));
                    [status, file_tree, ext] = this.navigateTree(file_tree, 'remote_check');
                    if status
                        this.log.addStatusOk('All files have been found remotely', 10);
                    else
                        this.log.addError('Some files not found remotely, starting to download what it''s available');
                        status = true;
                    end
                
                    if status || this.nrt
                        this.log.addMessage(this.log.indent('Downloading Resources ...'));
                        [status, ~] = this.navigateTree(file_tree, 'download');
                        if not(status)
                            this.log.addWarning('Not all file have been found or unconpressed');
                        end
                    end
                end
            end
            this.sys_c = this.state.cc.SYS_C(this.state.cc.active_list); % set sys_c again from constellation collector
        end
        
        function idx = getServerIdx(this, address , port)
            % Get idx of server int the FTP_downlader if not present open
            % the connection and append a FTP_downloder object
            %
            % SYNTAX:
            %   idx = this.getServerIdx(address, port)
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
        
        function [status, file_tree, ext] = navigateTree(this, file_tree, mode)
            % Navigate into th logical file tree (see Remote_resource_manager.getFileStr) and perform operations
            %
            % SYNTAX:
            %      [status, file_tree] = this.navigateTree(file_tree, mode)
            % INPUT:
            %     file_tree: structure containing the file tree and the
            %     logical operators
            %     mode: 'local_check' , 'remote_check' , 'download'
            status = false;
            ext = '';
            if iscell(file_tree) % if is a leaf
                if strcmp(file_tree{1}, 'null') || file_tree{2}
                    status = true;
                    if strcmp(mode, 'download') && ~strcmp(file_tree{1},'null') && file_tree{3} ~=0
                        loc_n = file_tree{3};
                        f_struct = this.rm.getFileLoc(file_tree{1});
                        f_name = f_struct.filename;
                        f_path = [f_struct.(['loc' sprintf('%03d',loc_n)]) f_name];
                        step_s = min(3*3600, this.fnp.getStepSec(f_path)); %supposing a polynomial of degree 12 and SP3 orbit data every 15 min (at worst)
                        dsa = this.date_start.getCopy();
                        dso = this.date_stop.getCopy();
                        dsa.addIntSeconds(-step_s);
                        dso.addIntSeconds(+step_s);
                        [file_name_lst, date_list] = this.fnp.dateKeyRepBatch(f_path, dsa, dso);
                        file_name_lst = flipud(file_name_lst);
                        status = true;
                        f_status_lst = file_tree{4};
                        f_ext_lst = file_tree{5};
                        
                        f_status_lst = Core_Utils.aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, date_list);
                        
                        for i = 1 : length(file_name_lst)
                            if isempty(f_status_lst) || ~f_status_lst(i)
                                file_name = file_name_lst{i};
                                server = regexp(file_name,'(?<=\?{)\w*(?=})','match','once'); % search for ?{server_name} in paths
                                file_name = strrep(file_name,['?{' server '}'],'');
                                [s_ip, port] = this.rm.getServerIp(server);
                                out_dir = this.state.getFileDir(file_name);
                                out_dir =  this.fnp.dateKeyRep(out_dir, date_list.getEpoch(date_list.length - i + 1));
                                if ~(exist(out_dir, 'file') == 7)
                                    mkdir(out_dir);
                                end
                                %out_dir = out_dir{1};
                                if instr(port,'21')
                                    idx = this.getServerIdx(s_ip, port);
                                    if ~this.nrt
                                        status = status && this.ftp_downloaders{idx}.downloadUncompress(file_name, out_dir);
                                    else
                                        status = this.ftp_downloaders{idx}.downloadUncompress(file_name, out_dir) && status;
                                    end
                                else
                                    if ~this.nrt
                                        status = status && Core_Utils.downloadHttpTxtResUncompress([s_ip file_name], out_dir);
                                    else
                                        status = Core_Utils.downloadHttpTxtResUncompress([s_ip file_name], out_dir) && status;
                                    end
                                end
                            end
                        end
                    end
                elseif ~strcmp(mode, 'download')
                    f_struct = this.rm.getFileLoc(file_tree{1});
                    if isempty(f_struct.filename)
                        this.log.addError(sprintf('File resource "%s" not found: remote_resource.ini seems to be corrupted', file_tree{1}));
                    else
                        f_name = f_struct.filename;
                        this.state.setFile(f_name, this.current_resource);
                        if strcmp(mode, 'local_check')
                            f_path = this.fnp.checkPath([this.state.getFileDir(f_name) filesep f_name]);
                            step_s = min(3*3600, this.fnp.getStepSec(f_path)); %supposing a polynomial of degree 12 and SP3 orbit data every 15 min (at worst)
                            dsa = this.date_start.getCopy();
                            dso = this.date_stop.getCopy();
                            dsa.addIntSeconds(-step_s);
                            dso.addIntSeconds(+step_s);
                            file_name_lst = flipud(this.fnp.dateKeyRepBatch(f_path, dsa, dso));
                            status = true;
                            f_status_lst = false(length(file_name_lst),1); % file list to be saved in tree with flag of downloaded or not
                            for i = 1 : length(file_name_lst)
                                file_info = dir(file_name_lst{i});
                                if isempty(file_info)
                                    f_status = false;
                                else
                                    % if the file is empty delete it
                                    if (file_info.bytes == 0)
                                        f_status = false;
                                        this.log.addError(sprintf('"%s" file is empty => deleting it...', [file_name_lst{i}]));
                                        delete(file_name_lst{i});
                                    else
                                        f_status = true;
                                    end
                                end
                                        
                                f_status_lst(i) = f_status;
                                status = status && f_status;
                                if f_status
                                    this.log.addStatusOk(sprintf('%s ready',this.fnp.getFileName(file_name_lst{i})), 20); % logging on 10 (default is 9, if ok do not show this)
                                else
                                    this.log.addWarning(sprintf('%s have not been found locally', this.fnp.getFileName(file_name_lst{i})));
                                end
                            end
                            if status
                                file_tree{3} = 0;
                            end
                            file_tree{4} = f_status_lst;
                        elseif strcmp(mode, 'remote_check')
                            
                            old_server = struct('name', '', 's_ip', '', 'port', '');
                            for i = 1 : f_struct.loc_number
                                f_path = [f_struct.(['loc' sprintf('%03d',i)]) f_name];
                                step_s = min(3*3600, this.fnp.getStepSec(f_path)); %supposing a polynomial of degree 12 and SP3 orbit data every 15 min (at worst)
                                dsa = this.date_start.getCopy();
                                dso = this.date_stop.getCopy();
                                dsa.addIntSeconds(-step_s);
                                dso.addIntSeconds(+step_s);
                                file_name_lst = flipud(this.fnp.dateKeyRepBatch(f_path, dsa, dso));
                                status = true;
                                f_status_lst = file_tree{4};
                                f_ext_lst = cell(numel(f_status_lst),1);
                                for j = 1 : length(file_name_lst)
                                    if ~f_status_lst(j)
                                        file_name = file_name_lst{j};
                                        [server] = regexp(file_name,'(?<=\?{)\w*(?=})','match'); % search for ?{server_name} in paths
                                        if isempty(server)
                                            this.log.addWarning(sprintf('No server is configured to download "%s"\nCheck remote_resources.ini', file_name));
                                            status = false;
                                        else
                                            server = server{1};
                                            file_name = strrep(file_name,['?{' server '}'],'');
                                            
                                            if strcmp(server, old_server.name)
                                                s_ip = old_server.s_ip;
                                                port = old_server.port;
                                            else
                                                [s_ip, port] = this.rm.getServerIp(server);
                                                old_server = struct('name', server, 's_ip', s_ip, 'port', port);
                                            end
                                            
                                            if instr(port,'21')
                                                idx = this.getServerIdx(s_ip, port);
                                                [stat, ext] = this.ftp_downloaders{idx}.check(file_name);
                                                if ~this.nrt
                                                    status = status && stat;
                                                else
                                                    status = stat;
                                                end
                                            else
                                                [stat, ext] = Core_Utils.checkHttpTxtRes([s_ip file_name]);
                                                if ~this.nrt
                                                    status = status && stat;
                                                else
                                                    status = stat;
                                                end
                                            end
                                            f_ext_lst{j} = ext;
                                            if status
                                                this.log.addStatusOk(sprintf('%s found (on remote server %s)', this.fnp.getFileName(file_name), server), 20);
                                            else
                                                if instr(port,'21')
                                                    this.log.addWarning(sprintf('"ftp://%s:%s%s" have not been found remotely', s_ip, port, file_name));
                                                else
                                                    this.log.addWarning(sprintf('"http://%s:%s%s" have not been found remotely', s_ip, port, file_name));
                                                end
                                                if ~this.nrt
                                                    break
                                                end
                                            end
                                        end
                                        
                                    end
                                end
                                file_tree{5} = f_ext_lst;                                
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
                    [status_b, file_tree_b, ext] = this.navigateTree(file_tree.(b_name).(branch{i}), mode);
                    if or_flag
                        status  = status || status_b;
                        if status
                            file_tree.(b_name).(branch{i}) = file_tree_b;
                            if ~this.nrt
                                break
                            end
                        end
                    elseif and_flag
                        status  = status && status_b;
                    end
                    file_tree.(b_name).(branch{i}) = file_tree_b;
                end
                
            end
        end
        
        function conjureFiles(this, date_start, date_stop, center_name)
            % Get all the files needed for processing
            %
            % SYNTAX:
            %     this.conjureFiles(date_start, date_stop, center_name)
            dsa = date_start.getCopy();
            dso = date_stop.getCopy();
            if (GPS_Time.now() - dso) < (24*3600)
                this.nrt = true;
            end
            if nargin < 4
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
                        if strcmp(split{2}, center_name)
                            is_ok = true;
                            this.center_name = center_name;
                        end
                    end
                end
            end
            
            if ~strcmp(center_name, 'none')
                this.state.setNoResources(false);
            end

            if ~is_ok
                if strcmp(center_name, 'none')
                    this.state.setNoResources(true);
                    this.log.addWarning('Resource center have not been selected, orbits will not be computed!!!');                    
                else
                    this.log.addError(['Selected center: ' center_name ' not compatible with selected constellations: ' this.sys_c]);
                    error('Ending execution: missing valid orbits')
                end
            end
                
            % Prepare all the files needed for processing
            
            if ~this.state.isNoResources()
                %this.state.updateNavFileName();
                %this.state.updateErpFileName();
                this.conjureNavFiles(dsa, dso);
                if this.state.isAutomaticDownload()
                    this.conjureDCBFiles(dsa, dso);
                    this.conjureCRXFiles(dsa, dso);
                end
                if this.state.needIonoMap() || this.state.isIonoKlobuchar()
                    this.conjureIonoFiles(dsa, dso, this.state.isIonoKlobuchar());
                end
                if this.state.isAtmLoading()
                    this.conjureAtmLoadFiles(dsa, dso);
                end
                if this.state.isVMF()
                    this.conjureVmfFiles(dsa, dso);
                end
            end
        end
        
        function conjureAtmLoadFiles(this, date_start, date_stop)
            this.log.addMarkedMessage('Checking Athmospheric loading files');
            status = this.conjureResource('atm_load',date_start, date_stop);

            if status
                this.log.addStatusOk('Atmospheric loading files are present ^_^');
            else
                this.log.addWarning('Not all atmospheric files founds, Atmospheric loading will not be applyied');
            end
            
        end
        
        function conjureVmfFiles(this, date_start, date_stop)
            this.log.addMarkedMessage('Checking VMF files');
            date_stop = date_stop.getCopy();
            %date_stop.addSeconds(6*3600);
            status = this.conjureResource('vmf',date_start, date_stop);

            if status
                this.log.addStatusOk('Vienna Mapping Function files are present ^_^');
            else
                this.log.addWarning('Not all vmf files founds');
            end
            
        end
        
        function conjureDCBFiles(this, date_start, date_stop)
            % Download of CAS .DCB files from the IGN server.
            %
            % SYNTAX:
            %   this.conjureDCBFiles(gps_week, gps_time);
            %
            % INPUT:
            %   date_start = starting GPS_Time
            %   date_stop = ending GPS_Time
            %
            % OUTPUT:
            
            legacy = false;
            if ~legacy
                status = this.conjureResource('bias', date_start, date_stop);
            else
                this.log.addMarkedMessage('Checking DCB files');
                if date_start.getCalEpoch >= 2013 % use CAS DCB
                    dcb_ok = true;
                    % check if file are present
                    fnp = File_Name_Processor();
                    ss = 'mxd';
                    archive = 'ign';
                    provider = 'cas';
                    dcb_type = 'final';
                    dcb_name = this.source.(archive).par.(ss).center.(provider).dcb.(dcb_type);
                    [~, dcb_file_name, ext] = fileparts(dcb_name);
                    this.state.setDcbFile([this.state.getDcbDir filesep dcb_file_name]);
                    tmp_date_start = date_start.getCopy;
                    tmp_date_stop = date_stop.getCopy;
                    file_list = fnp.dateKeyRepBatch(dcb_name, tmp_date_start, tmp_date_stop);
                    names = {};
                    dcb_ok = true(length(file_list),1);
                    for i = 1 : length(file_list)
                        [~, name, ext] = fileparts(file_list{i});
                        names{end+1} = name;
                        if exist(fnp.checkPath([this.state.getDcbDir filesep name]), 'file') ~= 2
                            dcb_ok(i) = false;
                        end
                    end
                    
                    if any(~dcb_ok)
                        aria_try = true;
                        try % ARIA2C download
                            aria_file_list = file_list;
                            f_ext_lst = cell(size(file_list));
                            for i = 1 : numel(file_list)
                                aria_file_list{i} = [this.source.(archive).ftpd.getFullAddress this.source.(archive).par.(ss).path file_list{i}];
                                [~, name, ext] = fileparts(file_list{i});
                                if strcmpi(ext, '.gz') || strcmpi(ext, '.Z')
                                    f_ext_lst{i} = ext;
                                    aria_file_list{i} = aria_file_list{i}(1 : end - length(ext));
                                end
                            end
                            dcb_ok = Core_Utils.aria2cDownloadUncompress(aria_file_list, f_ext_lst, dcb_ok, [], this.state.getDcbDir());
                            if any(~dcb_ok)
                                id_ko = find(dcb_ok == 0);
                                for i = id_ko'
                                    this.log.addWarning(sprintf('download of %s from %s failed, file not found or not accessible', [this.source.(archive).ftpd.getFullAddress this.source.(archive).par.(ss).path], names{i}));
                                end
                            end
                        catch ex
                            this.source.(archive).ftpd.download(this.source.(archive).par.(ss).path, file_list, this.state.getDcbDir());
                            for i = 1 : length(file_list)
                                if not(dcb_ok(i))
                                    [~, name, ext] = fileparts(file_list{i});
                                    if (isunix())
                                        system(['gzip -fd ' this.state.getDcbDir() filesep name ext]);
                                    else
                                        try
                                            [status, result] = system(['".\utility\thirdParty\7z1602-extra\7za.exe" -y x ' '"' this.state.getDcbDir() filesep name ext '"' ' -o' '"' this.state.getDcbDir() '"']); %#ok<ASGLU>
                                            delete([this.state.getDcbDir() filesep name ext]);
                                        catch
                                            this.log.addWarning(sprintf(['Please decompress the ' name ext ' file before trying to use it in goGPS.']));
                                            compressed = 1;
                                        end
                                    end
                                end
                            end
                        end
                    else
                        this.log.addStatusOk('DCB files are present ^_^');
                    end
                else % use DCB from CODE
                    gps_week = double([date_start.getGpsWeek; date_stop.getGpsWeek ]);
                    gps_time = [date_start.getGpsTime; date_stop.getGpsTime ];
                    
                    % Pointer to the global settings:
                    state = Core.getCurrentSettings();
                    
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
                    
                    this.log.addMessage(this.log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                    
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
                            year = [year (year_orig(y)).*ones(1,12)];
                        end
                        month = [month 1 : 1 : date_l(2)];
                        year  = [year date_l(1).*ones(1,date_l(2))];
                    end
                    
                    %connect to the DCB server
                    try
                        ftp_server = ftp(aiub_ip);
                    catch
                        this.log.addWarning('connection failed.\n');
                        this.state.setDcbFile({''});
                        return
                    end
                    
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
                                s2 = [ff{p} num2str(two_digit_year(year(y)),'%02d') num2str(month(m),'%02d') '.DCB.Z'];
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
                                                this.log.addWarning(sprintf(['Please decompress the ' s2 ' file before trying to use it in goGPS.']));
                                                compressed = 1;
                                            end
                                        end
                                        this.log.addMessage(this.log.indent(['DCB file downloaded: ' s2 ]));
                                    catch
                                        cd(ftp_server, '..');
                                        s1 = [ff{p} '.DCB'];
                                        mget(ftp_server,s1,down_dir);
                                        cd(ftp_server, num2str(year_orig(y)));
                                        s2 = [s2(1:end-2) '_TMP'];
                                        [move_success, message] = movefile([down_dir '/' s1], [down_dir '/' s2], 'f');
                                        if ~move_success
                                            log.addError(message);
                                        end
                                        this.log.addWarning(['Downloaded DCB file: ' s1 ' --> renamed to: ' s2]);
                                    end
                                else
                                    this.log.addMessage(this.log.indent([s2(1:end-2) ' already present\n']));
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
                    
                    this.log.addStatusOk('Dcb files have been downloded');
                end
            end
        end
        
        function conjureNavFiles(this, date_start, date_stop)
            % Wrapper of conjureResources for navigational files
            %
            % SYNTAX:
            %   this.conjureNavFiles(date_start, date_stop)
            %
            this.log.addMarkedMessage('Checking ephemerides / clocks / ERPs');
            list_preferred = this.state.getPreferredEph();
            for i = 1 : length(list_preferred)
                status = this.conjureResource(list_preferred{i}, date_start, date_stop);
                if status || (this.nrt && (strcmp(list_preferred{i},'ultra') || strcmp(list_preferred{i},'broadcast')))
                    break
                end
            end
            if status
                this.log.addStatusOk('Ephemerides files are present ^_^')
            else
                if ~this.nrt
                    this.log.addError('Not all ephemerides files have been found program might misbehave');
                else
                    this.log.addError('Not all ephemerides files have been found');
                end
            end
        end               
        
        function conjureIonoFiles(this, date_start, date_stop, flag_brdc)
            % Wrapper of conjureResources for iono files
            %
            % SYNTAX:
            %   this.conjureIonoFiles(date_start, date_stop)
            %
            this.log.addMarkedMessage('Checking ionospheric resources files');
            list_preferred = this.state.preferred_iono;
            for i = 1 : length(list_preferred)
                if strcmp(list_preferred{i}, 'broadcast')
                    status = this.conjureResource(['iono_' list_preferred{i}], date_start, date_stop);
                else
                    status = this.conjureResource(['iono_' list_preferred{i}], date_start, date_stop);
                    if flag_brdc
                        status = this.conjureResource('iono_broadcast', date_start, date_stop);
                    end
                end
                if status
                    break
                end
            end
            if (this.state.iono_model == 2) && (this.state.iono_management == 3 || this.state.flag_apr_iono)
                status = this.conjureResource(['iono_broadcast'], date_start, date_stop);
            end
            if status
                this.log.addStatusOk('Ionosphere resource files are present ^_^')
            else
                this.log.addWarning('Not all iono files found program might misbehave')
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
            
            log = Core.getLogger;
            this.log.addMarkedMessage('Checking CRX (Satellite problems file)');
            date_start = date_start.getCopy();
            date_stop = date_stop.getCopy();
            gps_week = double([date_start.getGpsWeek; date_stop.getGpsWeek ]);
            %gps_time = [date_start.getGpsTime; date_stop.getGpsTime ];
            
            % Pointer to the global settings:
            state = Core.getCurrentSettings();
            
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
                log.addMessage(log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                log.addError('Connection failed.');
                return
            end
            
            
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
                    log.addMessage(log.indent(sprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])));
                    %if not(exist([down_dir, '/', s2]) == 2)
                    try
                        move_success = false;
                        if exist(fullfile(down_dir, s2), 'file') == 2
                            [move_success, message] = movefile(fullfile(down_dir, s2), fullfile(down_dir, [s2 '.old']), 'f');
                            if ~move_success
                                log.addError(message);
                            end
                        end
                        try % ARIA2C download
                            clear file_name_lst f_ext_lst;
                            file_lst{1} = ['ftp://' aiub_ip s '/' s2];
                            f_ext_lst{1} = '';
                            f_status_lst = false;
                            f_status_lst = Core_Utils.aria2cDownloadUncompress(file_lst, f_ext_lst, f_status_lst, [], down_dir);
                        catch ex
                            % fprintf(ex.message)
                            mget(ftp_server,s2,down_dir);
                            f_status_lst = true;
                        end
                        if ~f_status_lst
                            throw(MException('Verify CRX download', 'download error'));
                        end
                    catch ex
                        log.addWarning(sprintf('CRX file have not been updated due to connection problems: %s', ex.message))
                        if exist(fullfile(down_dir, [s2, '.old']), 'file') == 2
                            if move_success
                                [move_success, message] = movefile(fullfile(down_dir, [s2 '.old']), fullfile(down_dir, s2));
                                if ~move_success
                                    log.addError(message);
                                end
                            end
                        end
                    end
                    log.addMessage(log.indent(sprintf(['Downloaded CRX file: ' s2 '\n'])));
                end
                
                % cell array with the paths to the downloaded files
                file_crx{y} = [down_dir, '/', s2];                                
            end
            
            try
                close(ftp_server);
            catch
            end
            log.addStatusOk('CRX files are present ^_^')
        end
        
    end
    
end
