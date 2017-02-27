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
%    |___/                    v 0.5.0
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
    
    properties (SetAccess = protected, GetAccess = protected)
        ftp_ip;       % IP address of the FTP server, stored as string
        remote_dir;   % Base dir path on the remote server to locate the file to download
        file_name;    % name of the file
        local_dir;    % Download directory: location on the local machine for the storage of the file to be downloaded
        
    end
                          
                              
    properties (SetAccess = private, GetAccess = public)
        logger = Logger.getInstance(); % Handler to the logger object
    end
    
   
    methods        
        function this = FTP_Downloader(ftp_ip, remote_dir, file_ame,  local_dir)
            % Constructor

        end        
    end
    
    methods (Static)
        function flag = checkNet()
            % Check whether internet connection is available            
            url =java.net.URL('http://www.google.org');
            
            % read the URL
            try
                link = openStream(url);
                flag = ~isempty(link);
            catch
                flag = false;
            end
        end
    end
    
end
