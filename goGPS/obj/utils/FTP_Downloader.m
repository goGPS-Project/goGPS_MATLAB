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
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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
        addr;           % IP address of the FTP server, stored as string
        port = '21';    % PORT to access the FTP server, stored as number
        remote_dir;     % Base dir path on the remote server to locate the file to download
        file_name;      % name of the file
        local_dir;      % Download directory: location on the local machine for the storage of the file to be downloaded
        ftp_server;     % object containing the connector
    end


    properties (SetAccess = private, GetAccess = public)
        log = Logger.getInstance(); % Handler to the log object
    end

    methods
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
                catch
                    this.ftp_server = [];
                end
            end
        end

        function delete(this)
            % destructor close the connection with the server
            close(this.ftp_server);
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
            url =java.net.URL('http://igs.org');

            % read the URL
            try
                link = openStream(url);
                flag = ~isempty(link);
            catch
                flag = false;
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
                log = Logger.getInstance();
                log.addError(ex.message);
            end
        end
    end

end
