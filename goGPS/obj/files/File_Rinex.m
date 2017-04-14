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
        base_dir = '../data/project/default_DD/RINEX/';                  % directory containing all the files
        file_name_list = {'yamatogawa_rover', 'yamatogawa_master'};      % file names (they can be multiple files for different days)
        ext = {'.obs', '.obs'};                                          % file names extension (they can be multiple files for different days)
        is_valid_list = false(1, 2);                 % for each element of file_name_list check the validity of the file
        
        is_composed = false;                         % when this flag is set, it means that the file_name depends on variables such as DOY DOW YYYY SSSS MM ecc...
        
        first_epoch = GPS_Time();                    % first epoch stored in the RINEX (updated after checkValidity)
        last_epoch = GPS_Time();                     % last epoch stored in the RINEX (updated after checkValidity)
        
        eoh = 0;                                     % end of header line - store the last line of the header
    end
    
    properties (SetAccess = private, GetAccess = private)
        id_date = 2:28;                              % Last character containing the 6 fields of the date (in case of 4digits year), it the pends on the type of rinex files (28 -> OBS RINEX 3)
    end
    
    properties (SetAccess = private, GetAccess = public)
    end
    
    methods
        function this = File_Rinex(base_dir, file_name, ext)
            % Creator of File_Rinex (base_dir, file_name, ext)
            %            File_Rinex (base_dir, file_name)
            %            File_Rinex (file_name)
            
            % fill the path with the imported file names
            switch (nargin)
                case 0 % only instantiate the object
                    return
                case 1 % populate from (file_name)
                    if iscellstr(base_dir)
                        this.file_name_list = {};
                        this.ext = {};
                        for f = 1 : numel(base_dir)
                            [this.base_dir, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(base_dir{f}));
                        end
                    else
                        [this.base_dir, file_name, ext] = fileparts(checkPath(fullfile(base_dir)));
                        this.file_name_list = {file_name};
                        this.ext = {ext};
                    end
                case 2 % populate from (base_dir, file_name)
                    if iscellstr(file_name)
                        this.file_name_list = {};
                        this.ext = {};
                        for f = 1 : numel(file_name)
                            [this.base_dir, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(fullfile(base_dir, file_name{f})));
                        end
                    else
                        [this.base_dir, file_name, ext] = fileparts(checkPath(fullfile(base_dir, file_name)));
                        this.file_name_list = {file_name};
                        this.ext = {ext};
                    end
                case 3 % populate from (base_dir, file_name, ext)
                    if (ext(1) ~= '.')
                        ext = ['~' ext];
                    end
                    if iscellstr(file_name)
                        this.file_name_list = {};
                        this.ext = {};
                        for f = 1 : numel(file_name)
                            [this.base_dir, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(fullfile(base_dir, [file_name{f} ext])));
                        end
                    else
                        [this.base_dir, file_name, ext] = fileparts(checkPath(fullfile(base_dir, [file_name ext])));
                        this.file_name_list = {file_name};
                        this.ext = {ext};
                    end
            end
            this.checkValidity();
        end
    end
    
    methods
        function checkValidity(this)
            % Update the status of validity of the files here pointed
            % SYNTAX: this.checkValidity();
            
            % pre allocate
            this.is_valid_list = false(1, numel(this.file_name_list));
            this.eoh = zeros(1, numel(this.file_name_list));
            this.first_epoch = GPS_Time();
            this.last_epoch = GPS_Time();
            % for each file present in the list
            for f = 1 : numel(this.file_name_list)
                full_path = fullfile(this.base_dir, [this.file_name_list{f} this.ext{f}]);
                
                % check the existence
                this.is_valid_list(f) = exist(full_path, 'file');
                if this.is_valid_list(f)
                    % try to find the first and the last epoch stored in the file
                    try
                        fid = fopen(fullfile(this.base_dir, [this.file_name_list{f} this.ext{f}]));
                        l = 1;
                        line = fgetl(fid);
                        while isempty(strfind(line,'END OF HEADER')) && ischar(line)
                            l = l + 1;
                            line = fgetl(fid);
                        end
                        this.eoh(f) = l;
                        epoch_line = fgetl(fid);
                        
                        % try to guess the time format
                        [id_start, id_stop] = regexp(epoch_line, '[.0-9]*');
                        this.id_date = id_start(1) : id_stop(6); % save first and last char limits of the date in the line -> suppose it composed by 6 fields
                        
                        this.first_epoch.addEpoch(epoch_line(this.id_date), [], true);
                        this.logger.addStatusOk(['"' this.file_name_list{f} this.ext{f} '" appears to be a valid RINEX']);
                        this.logger.addMessage(sprintf('        first epoch found at: %s', this.first_epoch.last.toString()));
                        
                        % go to the end of the file to search for the last epoch
                        % to be sure to find at least one line containing a valid epoch, go to the end of the file minus 5000 characters
                        fseek(fid,-10000,'eof');
                        fgetl(fid); % Probably i'm not at the beginning of a line -> disregard the first reading
                        % Start searching for a valid epoch
                        line = fgetl(fid);
                        while ischar(line)
                            % An epoch line has the second character containing the year of the observation
                            if (numel(line) > 2) && ~isempty(regexp(line(1:3),'(> [0-9].)|( [0-9].)', 'once'))
                                epoch_line = line;
                            end
                            line = fgetl(fid);
                        end
                        fclose(fid);
                        this.last_epoch.addEpoch(epoch_line(this.id_date), [], true);
                        this.logger.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()));
                        this.is_valid_list(f) = true;
                    catch ex
                        if this.first_epoch.lenght < f
                            this.first_epoch.addEpoch(0);
                        end
                        if this.last_epoch.lenght < f
                            this.last_epoch.addEpoch(0);
                        end
                        this.logger.addWarning(sprintf('"%s" appears to be a corrupted RINEX file - %s', full_path, ex.message()));
                        this.is_valid_list(f) = false;
                    end
                end
            end
            this.is_valid = all(this.is_valid_list);
            if (~this.is_valid)
                this.logger.addWarning('No valid RINEX found!!!');
            end
        end
        
        function file_name = getFileName(this, file_number)
            % Get the full path of a file (if the object contains a list of files, the id can be specified)
            % SYNTAX: file_name = this.getFileName(<file_number = 1>)
            if nargin == 1
                file_number = 1;
            end
            file_name = fullfile(this.base_dir, [this.file_name_list{file_number} this.ext{file_number}]);
        end        

        function line_num = getEOH(this, file_number)
            % Get the end of header line of a RINEX file (if the object contains a list of files, the id can be specified)
            % SYNTAX: line_num = this.getEOH(<file_number = 1>)
            if nargin == 1
                file_number = 1;
            end
            line_num = this.eoh(file_number);
        end
        
        function validity = isValid(this, file_number)
            % Get the validity of a RINEX file or the object (if the object contains a list of files, the id can be specified)
            % SYNTAX: validity = isValid(<file_number>)
            if (nargin == 1)
                validity = this.is_valid;
            else
                validity = this.is_valid_list(file_number);
            end
        end
                        
    end
end
