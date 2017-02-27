%   CLASS File_Rinex
% =========================================================================
%
% DESCRIPTION
%   Class to store file_paths for RINEX files
%
% EXAMPLE
%   settings = File_Rinex();
%
% FOR A LIST OF CONSTANTS and METHODS use doc File_Rinex

%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011 
%--------------------------------------------------------------------------

classdef File_Rinex < handle
    
    properties (SetAccess = protected, GetAccess = protected)
        logger = Logger.getInstance(); % Handler to the logger object
    end
    
    properties (SetAccess = protected, GetAccess = public)
       is_valid = false;                            % flag, if true it means that the object contains at least one valid rinex file
       base_dir = '../data/data_RINEX';             % directory containing all the files
       file_name_list = {'yamatogawa_rover', 'yamatogawa_master'};       % file names (they can be multiple files for different days)
       ext = '.obs';
       is_valid_list = false;                       % for each element of file_name_list check the validity of the file
       
       is_composed = false;                         % when this flag is set, it means that the file_name depends on variables such as DOY DOW YYYY SSSS MM ecc...              
       
       first_epoch = GPS_Time();                    % first epoch stored in the RINEX (updated after checkValidity)
       last_epoch = GPS_Time();                     % last epoch stored in the RINEX (updated after checkValidity)
    end
    
    properties (SetAccess = private, GetAccess = public)
    end
        
    methods
        function obj = File_Rinex(base_dir, file_name, ext)            
            % Creator of File_Rinex (base_dir, file_name, ext)
            %            File_Rinex (base_dir, file_name)
            %            File_Rinex (file_name)
            
            % fill the path with the imported file names
            switch (nargin)
                case 0 % only instantiate the object
                    return
                case 1 % populate from (file_name)
                    if iscellstr(base_dir)
                        for f = 1 : numel(base_dir)
                            [obj.base_dir, obj.file_name_list{f}, obj.ext] = fileparts(checkPath(base_dir{f}));
                        end
                    else
                        [obj.base_dir, file_name, obj.ext] = fileparts(checkPath(fullfile(base_dir)));
                        obj.file_name_list = {file_name};
                    end
                case 2 % populate from (base_dir, file_name)
                    if iscellstr(file_name)
                        for f = 1 : numel(file_name)
                            [obj.base_dir, obj.file_name_list{f}, obj.ext] = fileparts(checkPath(fullfile(base_dir, file_name{f})));
                        end
                    else
                        [obj.base_dir, file_name, obj.ext] = fileparts(checkPath(fullfile(base_dir, file_name)));
                        obj.file_name_list = {file_name};
                    end
                case 3 % populate from (base_dir, file_name, ext)
                    if (ext(1) ~= '.')
                        ext = ['~' ext];
                    end                        
                    if iscellstr(file_name)
                        for f = 1 : numel(file_name)
                            [obj.base_dir, obj.file_name_list{f}, obj.ext] = fileparts(checkPath(fullfile(base_dir, [file_name{f} ext])));
                        end
                    else
                        [obj.base_dir, file_name, obj.ext] = fileparts(checkPath(fullfile(base_dir, [file_name ext])));
                        obj.file_name_list = {file_name};
                    end
            end
            obj.checkValidity();
        end        
    end
    
    methods
        function checkValidity(obj)
            % Update the status of validity of the files here pointed
            
            % for each file present in the list
            for f = 1 : numel(obj.file_name_list)
                full_path = fullfile(obj.base_dir, [obj.file_name_list{f} obj.ext]);

                % check the existence
                obj.is_valid_list(f) = exist(full_path, 'file');
                if obj.is_valid_list(f)
                    % try to find the first and the last epoch stored in the file
                    try
                        fid = fopen(fullfile(obj.base_dir, [obj.file_name_list{f} obj.ext]));
                        line = fgetl(fid);
                        while isempty(strfind(line,'END OF HEADER')) && ischar(line)
                            line = fgetl(fid);
                        end
                        epoch_line = fgetl(fid);
                        
                        % this data conversion lines must be moved into a class GPS_Time
                        date   = cell2mat(textscan(epoch_line(2:28),'%f%f%f%f%f%f'));
                        if (date(1) < 80), date(1) = date(1) + 2000; end
                        obj.first_epoch.addEpoch(datenum(date), [], true);
                        obj.logger.addStatusOk(['"' obj.file_name_list{f} obj.ext '" appears to be a valid RINEX']);
                        obj.logger.addMessage(sprintf('        first epoch found at: %s', obj.first_epoch.last.toString()));
                        
                        % go to the end of the file to search for the last epoch
                        % to be sure to find at least one line containing a valid epoch, go to the end of the file minus 5000 characters
                        fseek(fid,-5000,'eof'); 
                        fgetl(fid); % Probably i'm not at the beginning of a line -> disregard the first reading
                        % Start searching for a valid epoch
                        line = fgetl(fid);
                        while ischar(line)
                            if (line(2) ~= ' ')
                                epoch_line = line;
                            end
                            line = fgetl(fid);
                        end
                        fclose(fid);
                        
                        % this data conversion lines must be moved into a class GPS_Time
                        date   = cell2mat(textscan(epoch_line(2:28),'%f%f%f%f%f%f'));
                        if (date(1) < 80), date(1) = date(1) + 2000; end
                        obj.first_epoch.addEpoch(datenum(date), [], true);
                        obj.logger.addMessage(sprintf('        last  epoch found at: %s', obj.first_epoch.last.toString()));

                    catch
                        obj.logger.addWarning(['"' full_path '" appears to be a corrupted RINEX file']);
                    end
                end
            end
            obj.is_valid = any(obj.is_valid_list);
            if (~obj.is_valid)
                obj.logger.addWarning('No valid RINEX found!!!');
            end
        end        
    end    
end
