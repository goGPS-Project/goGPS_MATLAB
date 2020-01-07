%   CLASS FTP_Downloader
% =========================================================================
%
% DESCRIPTION
%   Class to manage downloads from FTP of the various required files
%
% EXAMPLE
%   settings = FTP_Downloader();
%
% FOR A LIST OF CONSTANTs and METHODS use doc FTP_Downloader
%
% COMMENTS
%   Nothing to report

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
%  Contributors:     Andrea Gatti, Giulio Tagliaferro
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

classdef FTP_Downloader < handle
    properties (Constant, GetAccess = public)
        OK           =  0;   % All the files have been downloaded properly
        ERR_NEP      = -1;   % Not Enough parameters!
        ERR_NIC      = -2;   % No Internet Connection!
        ERR_FTP_FAIL = -3;   % Download failed for FTP problems
        ERR_FNF      = -4;   % at least one file not found
        ERR_FNV      = -5;   % at least one file not valid (it can not be uncompressed)
        W_FP         =  1    % Warning File Present - at least one file to be dowloaded is already present

    end

    properties (SetAccess = protected, GetAccess = protected)
        port = '21';    % PORT to access the FTP server, stored as number
        remote_dir;     % Base dir path on the remote server to locate the file to download
        file_name;      % name of the file
        local_dir;      % Download directory: location on the local machine for the storage of the file to be downloaded
        ftp_server;     % object containing the connector
        addr;           % IP address of the FTP server, stored as string
        f_name_pool = {};    % list with checked folder ant its names
    end


    properties (SetAccess = private, GetAccess = public)
        log = Core.getLogger();
        fnp = File_Name_Processor();
    end

    methods
        function addr = getAddress(this)
            addr = this.addr;
        end
        
        function addr = getFullAddress(this)
            addr = ['ftp://' this.addr ':' this.port '/' this.remote_dir];
        end
        
        function this = FTP_Downloader(ftp_addr, ftp_port, remote_dir, file_name,  local_dir)
            % Constructor
            % EXAMPLE: FTP_Downloader('')

            this.addr = ftp_addr;
            % Cleaning address
            if strcmp(this.addr(end),'/')
                this.addr(end) = [];
            end
            if strcmp(this.addr(end),'\')
                this.addr(end) = [];
            end
            this.addr = strrep(this.addr,'ftp://','');

            if (nargin > 1) && (~isempty(ftp_port))
                if ~ischar(ftp_port)
                    ftp_port = num2str(ftp_port);
                end
                this.port = ftp_port;
            end
            if (nargin > 2)
                this.remote_dir = remote_dir;
            end
            if (nargin > 3)
                this.file_name = file_name;
            end
            if (nargin > 4)
                this.local_dir = local_dir;
            end

            % Open the connection with the server
            if (this.checkNet)
                try
                    this.ftp_server = ftp(strcat(this.addr, ':', this.port));
                    cd(this.ftp_server);
                    warning('off')
                    sf = struct(this.ftp_server);
                    warning('on')
                    sf.jobject.enterLocalPassiveMode();
                catch
                    this.ftp_server = [];
                    this.log.addWarning(['Could not connect to: ' this.addr]);
                end
            end
        end

        function delete(this)
            % destructor close the connection with the server
            close(this.ftp_server);
        end
        
        function [status, ext]  = check(this, filepath)
            folder = this.fnp.getPath(filepath);
            f_name = this.fnp.getFileName(filepath);
            not_in_cache = true;
            for i = 1 : length(this.f_name_pool)
                folder_c = this.f_name_pool{i};
                if strcmp(folder, folder_c{1})
                    idx = i;
                    not_in_cache = false;
                end
            end
            ext = '';
            if not_in_cache
                try
                    files = dir(this.ftp_server, folder);
                    this.f_name_pool{end+1} = {folder , files};
                    idx = length(this.f_name_pool);
                catch
                    status = false;
                    this.log.addWarning(['Could not connect to: ' this.addr]);
                    return
                end
            end
            folder_s = this.f_name_pool{idx};
            files = folder_s{2};
            for i = 1 : length(files)
                if ~isempty(strfind(files(i).name, f_name))
                    [~, ~, ext_s] = fileparts(files(i).name); 
                    [~, ~, ext] = fileparts(f_name);
                    % If the file on the server is compressed return the extension of the compression
                    ext = iif(strcmp(ext, ext_s), '', ext_s);
                    status = true;
                    return
                end
            end
            status = false;
            
        end
        
        function [status]  = downloadUncompress(this, filepath, out_dir)
            %try
                path = File_Name_Processor.getPath(filepath);
                fname = File_Name_Processor.getFileName(filepath);
                cd(this.ftp_server, path);
                this.log.addMessage(this.log.indent(sprintf('downloading %s ...',fname)));
                if ~(7 ==exist(out_dir,'dir'))
                    mkdir(File_Name_Processor.getFullDirPath(out_dir));
                end
                retry = 0;
                status = false;
                while(retry < 2)
                    try
                        try
                            fpath = mget(this.ftp_server, [fname '*'], out_dir); % some ftp do not work properly with "*" wildcard -> try standard names too
                        catch
                            fpath = [];
                        end
                        if isempty(fpath)
                            try
                                fpath = mget(this.ftp_server, [fname '.Z'], out_dir);
                            catch
                                fpath = [];
                            end
                        end
                        if isempty(fpath)
                            try
                                fpath = mget(this.ftp_server, [fname '.gz'], out_dir);
                            catch
                                fpath = [];
                            end
                        end
                        if isempty(fpath)
                            try
                                fpath = mget(this.ftp_server, fname, out_dir);
                            catch
                                fpath = [];
                            end
                        end
                        retry = 10;
                        
                        if isempty(fpath)
                            status = false;
                            if this.log.isScreenOut; fprintf('\b'); end
                            this.log.addMessage(' Failed');
                            return
                        else
                            status = true;
                            if this.log.isScreenOut; fprintf('\b'); end
                            this.log.addMessage(' Done');
                            [~, ~, fext] = fileparts(fpath{1});
                            if strcmp(fext,'.Z') || strcmp(fext,'.gz')
                                if (isunix())
                                    system(['gzip -d -f ' fpath{1} '&> /dev/null &']);
                                else
                                    try
                                        [status, result] = system(['.\utility\thirdParty\7z1602-extra\7za.exe -y x '  fpath{1} ' -o'  out_dir ]); %#ok<ASGLU>
                                        if (status == 0)
                                            status = true;
                                            delete([fpath{1}]);
                                        end
                                    catch
                                        this.log.addError(sprintf('Please decompress the %s file before trying to use it in goGPS!!!', fname));
                                        status = false;
                                    end
                                end
                            end
                            
                        end
                    catch ex
                        pause(0.5 * (retry - 1));
                        this.log.addError(sprintf('%s - %s', fname, ex.message));
                        retry = retry + 1;
                    end
                end
                %             catch
%             end
        end

        function [status, compressed] = download(this, remote_dir, file_name, local_dir, force_overwrite)
            % function to download a file (or a list of files) from a ftp a server
            % SYNTAX:
            %   status = this.download(<remote_dir>, file_name, local_dir, <force_overwrite = true>)
            %   status = this.download(remote_dir, file_name, local_dir)
            %   status = this.download(file_name, local_dir)
            %   status = this.download(file_name, local_dir, force_overwrite)
            %
            % EXAMPLE:
            %   ftpd = FTP_Downloader('127.0.0.1','25','/');
            %   ftpd.download('file.txt', './');

            % managing function overloading
            if (nargin < 3)
                this.log.addError('FTP_Downloader.download, not enough input parameters');
                status = FTP_Downloader.ERR_NEP;
                return;
            elseif (nargin == 3)
                force_overwrite = false;
                local_dir = file_name;
                this.file_name = remote_dir;
            elseif (nargin == 4)
                if ~ischar(local_dir)
                    force_overwrite = local_dir;
                    local_dir = file_name;
                    this.file_name = remote_dir;
                    this.remote_dir = FTP_Downloader.remote_dir;
                else
                    force_overwrite = false;
                    this.remote_dir = remote_dir;
                    this.file_name = file_name;
                end
            else
                this.remote_dir = remote_dir;
                this.file_name = file_name;
            end

            % function start
            if (this.checkNet)
                % connect to the server
                try
                    this.log.addMarkedMessage(sprintf('Initializing download process from %s', strcat(this.addr, ':', this.port)));
                    this.log.newLine();

                    % Try to get the current directory to check the ftp connection
                    try
                        cd(this.ftp_server);
                    catch
                        this.log.addWarning('connected with remote FTP has been closed, trying to re-open it');
                        this.ftp_server = ftp(strcat(this.addr, ':', this.port));
                    end
                    
                    warning('off')
                    sf = struct(this.ftp_server);
                    warning('on')
                    sf.jobject.enterLocalPassiveMode();
                    
                    % convert file_name in a cell array
                    if ~iscell(this.file_name)
                        this.file_name = {this.file_name};
                    end

                    n_files = numel(this.file_name);
                    status = FTP_Downloader.OK;
                    i = 0;
                    while (i < n_files) && (status >= FTP_Downloader.OK)
                        i = i + 1;

                        % get the exact remote path / file name
                        full_path = strcat(this.remote_dir, this.file_name{i});
                        if iscell(full_path)
                            full_path = full_path{1};
                        end
                        if (full_path(1) ~= '/')
                            full_path = strcat('/', full_path);
                        end
                        [remote_dir, file_name, file_ext] = fileparts(full_path);
                        compressed = strcmpi(file_ext,'.Z');

                        % check for a local version
                        local_path = fullfile(local_dir, iif(compressed, file_name, strcat(file_name, file_ext)));
                        local_file = exist(local_path, 'file');

                        if ~local_file || force_overwrite
                            try
                                % check file existence (without extension)
                                file_exist = ~isempty(dir(this.ftp_server, full_path));
                                if (~compressed && ~file_exist)
                                    full_path = strcat(full_path,'.Z');
                                    file_ext = strcat(file_ext,'.Z');
                                    file_exist = ~isempty(dir(this.ftp_server, full_path));
                                    compressed = true;
                                end
                                file_name = strcat(file_name, file_ext);

                                if (file_exist)
                                    this.log.addStatusOk(sprintf('%s found in %s, downloading...',file_name, remote_dir));

                                    try
                                        if ~(exist(local_dir,'dir'))
                                            mkdir(local_dir);
                                        end
                                        % move to the remote dir of the file
                                        cd(this.ftp_server, remote_dir);
                                        mget(this.ftp_server, file_name, local_dir);
                                        close(this.ftp_server);
                                        this.ftp_server = ftp(strcat(this.addr, ':', this.port));
                                        if compressed
                                            try
                                                if (isunix())
                                                    system(['uncompress -f ' local_dir filesep file_name]);
                                                    compressed = false;
                                                else
                                                    try
                                                        [status, result] = system(['".\utility\thirdParty\7z1602-extra\7za.exe" -y x ' '"' local_dir filesep file_name '"' ' -o' '"' local_dir '"']); %#ok<ASGLU>
                                                        delete([local_dir filesep file_name]);
                                                    catch
                                                        this.log.addError(sprintf('Please decompress the %s file before trying to use it in goGPS!!!', file_name));
                                                        compressed = 1;
                                                    end
                                                end
                                                this.file_name{i} = file_name(1:end-2);
                                                this.log.addMessage(sprintf('\b\b file ready!'));
                                            catch
                                                this.log.addWarning(sprintf('decompression of %s from %s failed', file_name, strcat(this.addr, ':', this.port)));
                                                status = FTP_Downloader.ERR_FNV;
                                            end
                                        else
                                            this.log.addMessage(sprintf('\b\b file ready!'));
                                        end
                                    catch
                                        this.log.addWarning(sprintf('download of %s from %s failed', file_name, strcat(this.addr, ':', this.port, remote_dir)));
                                        this.log.addWarning('file not found or not accessible');
                                        status = FTP_Downloader.ERR_FNF;
                                    end
                                else
                                    this.log.addWarning(sprintf('download of %s from %s failed', file_name, strcat(this.addr, ':', this.port, remote_dir)));
                                    this.log.addWarning('file not found or not accessible');
                                    status = FTP_Downloader.ERR_FNF;
                                end
                            catch ex
                                this.log.addWarning(sprintf('connection to %s failed', strcat(this.addr, ':', this.port, remote_dir)));
                                this.log.addWarning(sprintf('%s', ex.message));
                                status = FTP_Downloader.ERR_FNF;
                            end
                        else
                            this.log.addStatusOk(sprintf('%s has been found locally', strcat(file_name, file_ext)));
                            status = FTP_Downloader.W_FP;
                        end
                    end
                catch ex
                    this.log.addWarning(sprintf('connection to %s failed - %s', strcat(this.addr, ':', this.port), ex.message));
                    status = FTP_Downloader.ERR_FTP_FAIL;
                end
            else
                status = FTP_Downloader.ERR_NIC;
            end
            this.log.newLine();
        end
    end

    methods (Static)
        function flag = checkNet()
            % Check whether internet connection is available
            if ispc
                flag = ~system('ping -n 1 www.fast.com > nul');
            elseif isunix
                flag = ~system('ping -c 1 www.fast.com > /dev/null 2>&1');
            end
        end

        function [file, status, compressed] = urlRead(ftp_addr, ftp_port, remote_dir, file_name, local_dir)
            % download and read an ftp file
            % SYNTAX: [file, status, compressed] = urlRead(ftp_addr, ftp_port, remote_dir, file_name, local_dir)
            status = -1;
            compressed = false;
            try
                ftpd = FTP_Downloader(ftp_addr, ftp_port, remote_dir, file_name,  local_dir);
                [status, compressed] = ftpd.download(remote_dir, file_name, local_dir);
                fid = fopen(fullfile(local_dir, file_name), 'r');
                file = fread(fid);
                fclose(fid);
            catch ex
                file = '';
                log = Core.getLogger();
                log.addError(ex.message);
            end
        end
    end

end
