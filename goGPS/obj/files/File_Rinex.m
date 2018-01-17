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
        log = Logger.getInstance(); % Handler to the log object
    end

    properties (SetAccess = protected, GetAccess = public)
        is_valid = false;                            % flag, if true it means that the object contains at least one valid rinex file
        base_dir = {'../data/project/default_DD/RINEX/'};                  % directory containing all the files
        file_name_list = {'yamatogawa_rover', 'yamatogawa_master'};      % file names (they can be multiple files for different days)
        ext = {'.obs', '.obs'};                                          % file names extension (they can be multiple files for different days)
        is_valid_list = false(1, 2);                 % for each element of file_name_list check the validity of the file

        is_composed = false;                         % when this flag is set, it means that the file_name depends on variables such as DOY DOW YYYY SSSS MM ecc...

        first_epoch = GPS_Time();                    % first epoch stored in the RINEX (updated after checkValidity)
        last_epoch = GPS_Time();                     % last epoch stored in the RINEX (updated after checkValidity)

        verbosity_lev = 0;                          % Level of verbosity  0 = all messages, see Logger for more information

        eoh = 0;                                     % end of header line - store the last line of the header
    end

    properties (SetAccess = private, GetAccess = private)
        id_date = 2:28;                              % Last character containing the 6 fields of the date (in case of 4digits year), it the pends on the type of rinex files (28 -> OBS RINEX 3)
    end

    methods
        function this = File_Rinex(file_name, verbosity_lev)
            % Creator of File_Rinex (file_name)

            % fill the path with the imported file names
            if ~iscellstr(file_name)
                file_name = {file_name};
            end
            this.file_name_list = {};
            this.ext = {};
            for f = 1 : numel(file_name)
                [this.base_dir{f}, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(file_name{f}));
            end

            if nargin == 2
                this.verbosity_lev = verbosity_lev;
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
                full_path = fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]);

                % check the existence
                this.is_valid_list(f) = exist(full_path, 'file');
                if this.is_valid_list(f)
                    % try to find the first and the last epoch stored in the file
                    try
                        fid = fopen(fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]));
                        l = 1;
                        line = fgetl(fid);
                        while isempty(strfind(line,'END OF HEADER')) && ischar(line) %#ok<*STREMP>
                            l = l + 1;
                            line = fgetl(fid);
                        end
                        this.eoh(f) = l;
                        epoch_line = fgetl(fid);

                        % try to guess the time format
                        [id_start, id_stop] = regexp(epoch_line, '[.0-9]*');
                        this.id_date = id_start(1) : id_stop(6); % save first and last char limits of the date in the line -> suppose it composed by 6 fields

                        this.first_epoch.addEpoch(epoch_line(this.id_date), [], true);
                        this.log.addStatusOk(['"' this.file_name_list{f} this.ext{f} '" appears to be a valid RINEX'], this.verbosity_lev);
                        this.log.addMessage(sprintf('        first epoch found at: %s', this.first_epoch.last.toString()), this.verbosity_lev);

                        % go to the end of the file to search for the last epoch
                        % to be sure to find at least one line containing a valid epoch, go to the end of the file minus 5000 characters
                        fseek(fid,-10000,'eof');
                        fgetl(fid); % Probably i'm not at the beginning of a line -> disregard the first reading
                        % Start searching for a valid epoch
                        line = fgetl(fid);
                        while ischar(line)
                            % An epoch line has the second character containing the year of the observation
                            % e.g. RINEX 3:     " 15  8 23  0  0  0.0000000  0  8G05G07G28G02G06G09G30G13"
                            %      RINEX 3 NAV: "G01 2006 10 01 00 00 00 0.798045657575E-04 0.227373675443E-11 0.000000000000E+00"
                            %      RINEX 2:     "> 2016  7 18  0  1  0.0000000  0 36"
                            %      RINEX 2 NAV: " 3 98  2 15  0 15  0.0 0.163525342941D-03 0.363797880709D-11 0.108000000000D+05"
                            %      RINEX 3 MET: " 1996  4  1  0  0 15  987.1   10.6   89.5"
                            % this check could not work when comment are present after the header
                            if (numel(line) > 20) && ~isempty(regexp(line(1:10),'( [0-9]{2,4} )', 'once'))
                                epoch_line = line;
                            end
                            line = fgetl(fid);
                        end
                        fclose(fid);
                        this.last_epoch.addEpoch(epoch_line(this.id_date), [], true);
                        this.log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                        this.is_valid_list(f) = true;
                    catch ex
                        if this.first_epoch.length < f
                            this.first_epoch.addEpoch(0);
                        end
                        if this.last_epoch.length < f
                            this.last_epoch.addEpoch(0);
                        end
                        this.log.addWarning(['"' this.file_name_list{f} this.ext{f} '" appears to be a corrupted RINEX file']);
                        this.is_valid_list(f) = false;
                    end
                else
                    this.log.addError(['"' this.file_name_list{f} this.ext{f} '" appears to be missing']);
                    this.is_valid_list(f) = false;
                end
            end
            this.is_valid = all(this.is_valid_list);
            if (~this.is_valid)
                this.log.addWarning('Some or all the RINEX files are corrupted or missing!!!');
            end
        end

        function file_name = getFileName(this, file_number)
            % Get the full path of a file (if the object contains a list of files, the id can be specified)
            % SYNTAX: file_name = this.getFileName(<file_number = 1>)
            if nargin == 1
                file_number = 1;
            end
            file_name = fullfile(this.base_dir{file_number}, [this.file_name_list{file_number} this.ext{file_number}]);
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
