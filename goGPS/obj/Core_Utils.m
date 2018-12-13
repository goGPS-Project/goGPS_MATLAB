
%   CLASS Core_Utils
% =========================================================================
%
% DESCRIPTION
%   Class to manages utilities
%
% EXAMPLE
%   % set of static utilities
%   Core_Utils.diffAndPred
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, Giulio Tagliaferro, ...
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

classdef Core_Utils < handle
    properties
    end
    
    methods (Static)
        function exportCurFig(out_path)
            % EXAMPLE
            %   Core_Utilis.exportCurFig(fullfile('/Users/Andrea/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/GIMS/title.png'))
            fh = gcf; fh.WindowStyle = 'normal'; export_fig(fh, out_path, '-transparent', '-r150'); fh.WindowStyle = 'docked';
        end
        
        function diff_data = diffAndPred(data, n_order, t_ref)
            % compute diff predicting epoch 0 of each arc
            % using interp 1 pchip method
            %
            % SYNTAX
            %   Core_Utils.diffAndPred(data, t_ref)
            
            if nargin < 3
                t_ref = 1 : size(data,1);
            end
            if nargin < 2
                n_order = 1;
            end
            data = [repmat(data(1,:), n_order, 1); data];
            for s = 1 : size(data, 2)
                tmp = data(1 + n_order : end, s);
                id_ok = ~isnan(tmp);
                if sum(id_ok) > 2
                    data(1 : n_order, s) = interp1(t_ref(id_ok), tmp(id_ok), 1 - n_order : 0, 'pchip', 'extrap');
                    % interpolate the "left" of the first element of an arc
                    % because diff "eat" the first value
                    lim = getOutliers(id_ok);
                    lim_short = lim(lim(:,2) - lim(:,1) < 2 & lim(:,1) > 1, :);
                    for l = 1 : size(lim_short, 1)
                        data(lim_short(l, 1), s) = data(lim_short(l, 1)+1, s);
                    end
                    lim = lim(lim(:,2) - lim(:,1) > 2 & lim(:,1) > 1, :);
                    for l = 1 : size(lim, 1)
                        id_data = lim(l, 1) : lim(l, 2);
                        data(lim(l, 1), s) = interp1(t_ref(id_data), tmp(id_data), lim(l, 1) - 1, 'pchip', 'extrap');
                    end
                end
            end
            diff_data = diff(data, n_order);
        end
        
        function idx = findMO(find_list, to_find_el)
            % find the postion of the elements of to_find_el into find_list
            % find list should have unique elements
            idx = zeros(size(to_find_el));
            for i = 1: length(to_find_el)
                idx(i) = find(find_list == to_find_el(i),1);
            end
        end
        
        function [amb_idx, n_amb] = remEmptyAmbIdx(amb_idx, n_amb)
            % remove emtpy amb_idx
            %
            % SYNTAX:
            % amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx, <n_amb>)
            if nargin < 2
                n_amb = max(max(amb_idx));
            end
            i = 1;
            while (i <= n_amb)
                n_ep_amb = sum(sum(amb_idx == i));
                if n_ep_amb == 0
                    n_amb = n_amb - 1;
                    amb_idx(amb_idx > i) = amb_idx(amb_idx > i) - 1;
                else
                    i = i + 1;
                end
            end
        end
        
        function num = round_even(num)
           num = round((num-2)/2)*2+2;
        end
        
        function num = round_odd(num)
           num = round((num-1)/2)*2+1;
        end

        function num = code2Char2Num(str2)
            % Convert a 2 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str2(:,1:2) * [2^8 1]';
        end
        
        function str2 = num2Code2Char(num)
            % Convert a numeric value (float) of a 2 char string
            % SYNTAX
            %   str3 = Core_Utils.num2Code3ch(num)
            str2 = char(zeros(numel(num), 2));
            str2(:,1) = char(floor(num / 2^8));
            num = num - str2(:,2) * 2^8;
            str2(:,2) = char(num);
        end

        function num = code3Char2Num(str3)
            % Convert a 3 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code3ch2Num(str3);
            
            num = str3(:,1:3) * [2^16 2^8 1]';
        end
        
        function str3 = num2Code3Char(num)
            % Convert a numeric value (float) of a 3 char string
            % SYNTAX
            %   str3 = Core_Utils.num2Code3ch(num)
            str3 = char(zeros(numel(num), 3));
            str3(:,1) = char(floor(num / 2^16));
            num = num - str3(:,1) * 2^16;
            str3(:,2) = char(floor(num / 2^8));
            num = num - str3(:,2) * 2^8;
            str3(:,3) = char(num);
        end
        
        function num = code4Char2Num(str4)
            % Convert a 4 char string into a numeric value (float)
            % SYNTAX
            %   num = Core_Utils.code4ch2Num(str4);
            
            num = str4(:,1:4) * [2^24 2^16 2^8 1]';
        end
        
        function str4 = num2Code4Char(num)
            % Convert a numeric value (float) of a 4 char string
            % SYNTAX
            %   str4 = Core_Utils.num2Code4Char(num)
            str4 = char(zeros(numel(num), 4));
            str4(:,1) = char(floor(num / 2^24));
            num = num - str4(:,1) * 2^24;
            str4(:,2) = char(floor(num / 2^16));
            num = num - str4(:,2) * 2^16;
            str4(:,3) = char(floor(num / 2^8));
            num = num - str4(:,3) * 2^8;
            str4(:,4) = char(num);
        end
        
        function str4 = unique4ch(str4)
            % Perform unique on an array of 4 char codes
            %
            % SYNTAX
            %   str4 = Core_Utilis.unique4ch(str4)
            str4 = Core_Utils.num2Code4Char(unique(Core_Utils.code4Char2Num(str4)));
        end
        
        function str3 = unique3ch(str3)
            % Perform unique on an array of 3 char codes
            %
            % SYNTAX
            %   str3 = Core_Utilis.unique3ch(str3)
            str3 = Core_Utils.num2Code3Char(unique(Core_Utils.code3Char2Num(str3)));
        end
        
        function str2 = unique2ch(str2)
            % Perform unique on an array of 2 char codes
            %
            % SYNTAX
            %   str2 = Core_Utilis.unique2ch(str3)
            str2 = Core_Utils.num2Code2Char(unique(Core_Utils.code2Char2Num(str2)));
        end
        
        function [antenna_PCV] = readAntennaPCV(filename, antmod, date_limits)
            % SYNTAX:
            %   [antPCV] = this.readAntennaPCV(filename, antmod, date_start, date_stop);
            %
            % INPUT:
            %   filename    = antenna phase center offset/variation file
            %   antmod      = cell-array containing antenna model strings
            %   date_limits = in GPS_Time to filter useful data (matlab format preferred)
            %
            % OUTPUT:
            %   antenna_PCV (see description below)
            %
            % DESCRIPTION:
            %   Extracts antenna phase center offset/variation values from a PCO/PCV file in ATX format.
            %
            % RESOURCES:
            %   ftp://igs.org/pub/station/general/antex14.txt
            %
            % antenna_PCV struct definition
            % antenna_PCV.name           : antenna name (with radome code)
            % antenna_PCV.n_frequency    : number of available frequencies
            % antenna_PCV.frequency_name : array with name of available frequencies ({'G01';'G02';'R01',...})
            % antenna_PCV.frequency      : array with list of frequencies (carrier number) corresponding to the frequencies name ({'1';'2';'1',...})
            % antenna_PCV.sys            : array with code id of the system constellation of each frequency (1: GPS, 2: GLONASS, ...)
            % antenna_PCV.sysfreq        : array with codes of the system constellation and carrier of each frequency (11: GPS L1, 12: GPS L2, 21: GLONASS L1, ...)
            % antenna_PCV.offset         : ENU (receiver) or NEU (satellite) offset (one array for each frequency)
            % antenna_PCV.dazi           : increment of the azimuth (0.0 for non-azimuth-dependent phase center variations)
            % antenna_PCV.zen1           : Definition of the grid in zenith angle: minimum zenith angle
            % antenna_PCV.zen2           : Definition of the grid in zenith angle: maximum zenith angle
            % antenna_PCV.dzen           : Definition of the grid in zenith angle: increment of the zenith angle
            % antenna_PCV.tableNOAZI     : PCV values for NOAZI, in a cell array with a vector for each frequency [m]
            % antenna_PCV.tablePCV       : PCV values elev/azim depentend, in a cell array with a matrix for each frequency [m]
            % antenna_PCV.tablePCV_zen   : zenith angles corresponding to each column of antenna_PCV.tablePCV
            % antenna_PCV.tablePCV_azi   : azimutal angles corresponding to each row of antenna_PCV.tablePCV
            
            log = Core.getLogger();
            
            for m = numel(antmod) : -1 : 1
                antenna_PCV(m) = struct('name', antmod{m}, ...
                    'sat_type',[] ,...
                    'n_frequency', 0, ...
                    'available', 0, ...
                    'type', '', ...
                    'dazi', 0, ...
                    'zen1', 0, ...
                    'zen2', 0, ...
                    'dzen', 0, ...
                    'offset', [], ...
                    'frequency_name', [], ...
                    'frequency', [], ...
                    'sys', [], ...
                    'sysfreq', [], ...
                    'tableNOAZI', [], ...
                    'tablePCV_zen', [], ...
                    'tablePCV_azi', [], ...
                    'tablePCV', []);
            end
            antenna_found = zeros(length(antmod),1);
            
            % for each PCV file
            for file_pcv = 1 : size(filename, 1)
                if sum(antenna_found) < length(antmod)
                    if (~isempty(filename))
                        fid = fopen(char(filename(file_pcv, :)),'r');
                        if (fid ~= -1)
                            atx_file = textscan(fid,'%s','Delimiter', '\n', 'whitespace', '');
                            atx_file = atx_file{1};
                            fclose(fid);
                            
                            found = 0;
                            format = 0;
                            % get format (1: ATX, 2: Bernese 5.0, 3: Bernese 5.2)
                            l = 1;
                            line = atx_file{l};
                            if ~isempty(strfind(line, 'ANTEX VERSION / SYST'))
                                format = 1;
                            end
                            if ~isempty(strfind(line, 'MODEL NAME:'))
                                format = 2;
                            end
                            if ~isempty(strfind(line, 'ANTENNA PHASE CENTER VARIATIONS DERIVED FROM ANTEX FILE'))
                                format = 3;
                            end
                            
                            switch format
                                % ATX
                                case 1
                                    flag_stop = 0;
                                    ant_char = strcat(antmod{:});
                                    while (l < numel(atx_file) && found < length(antmod) && ~flag_stop)
                                        % go to the next antenna
                                        line = atx_file{l};
                                        while (l < numel(atx_file)-1) && ((length(line) < 76) || isempty(strfind(line(61:76),'START OF ANTENNA')))
                                            l = l + 1; line = atx_file{l};
                                        end
                                        l = l + 1; line = atx_file{l};
                                        
                                        if ~isempty(strfind(line,'TYPE / SERIAL NO')) %#ok<*STREMP> % antenna serial number
                                            if (nargin == 2) % receiver
                                                id_ant = strfind(strrep(ant_char,' ',''),strrep(line(1:20),' ',''));
                                                sat_type=[];
                                                 
                                            else
                                                id_ant = strfind(ant_char, line(21:23));
                                                sat_type=strtrim(line(1:20));
                                            end
                                            if ~isempty(id_ant)
                                                if (nargin == 2) % receiver
                                                    m = (id_ant - 1) / 20 + 1; % I'm reading the antenna
                                                else
                                                    m = (id_ant - 1) / 3 + 1; % I'm reading the antenna
                                                end
                                                
                                                if ~(antenna_PCV(m(1)).available)
                                                    
                                                    for a = 1:length(m)
                                                        log.addMessage(sprintf('Reading antenna %d => %s', m(a), antmod{m(a)}),100);
                                                    end
                                                    
                                                    invalid_date = 0;
                                                    
                                                    validity_start = [];
                                                    validity_end   = [];
                                                    
                                                    l_start = l; % line at the beginng of the antenna section
                                                    % look for "VALID FROM" and "VALID UNTIL" lines (if satellite antenna)
                                                    if (nargin > 2)
                                                        while (isempty(strfind(line,'VALID FROM')))
                                                            l = l + 1; line = atx_file{l};
                                                        end
                                                        validity_start = [str2num(line(3:6)) str2num(line(11:12)) str2num(line(17:18)) str2num(line(23:24)) str2num(line(29:30)) str2num(line(34:43))]; %#ok<*ST2NM>
                                                        l = l + 1; line = atx_file{l};
                                                        if (strfind(line, 'VALID UNTIL')) %#ok<*STRIFCND>
                                                            validity_end = [str2num(line(3:6)) str2num(line(11:12)) str2num(line(17:18)) str2num(line(23:24)) str2num(line(29:30)) str2num(line(34:43))];
                                                        else
                                                            validity_end = Inf;
                                                        end
                                                    end
                                                    
                                                    if (~isempty(validity_start)) % satellite antenna
                                                        if ~(date_limits.first.getMatlabTime() > datenum(validity_start) && (date_limits.last.getMatlabTime() < datenum(validity_end)))
                                                            invalid_date = 1;
                                                            antenna_PCV(m(1)).n_frequency = 0;
                                                            if isinf(validity_end)
                                                                log.addMessage(sprintf(' - out of range -> (%s : %s) not after %s', date_limits.first.toString(), date_limits.last.toString(), datestr(validity_start)), 100)
                                                            else
                                                                log.addMessage(sprintf(' - out of range -> (%s : %s) not intersecting (%s : %s)', date_limits.first.toString(), date_limits.last.toString(), datestr(validity_start), datestr(validity_end)), 100)
                                                            end
                                                        end
                                                    else  %receiver antenna
                                                    end
                                                    
                                                    if ~(invalid_date) % continue parsing
                                                        for a = 1:length(m)
                                                            log.addMessage(sprintf('Found a valid antenna %s', antmod{m(a)}), 50);
                                                        end
                                                        l = l_start;
                                                        
                                                        % get TYPE
                                                        antenna_PCV(m(1)).type = line(1:20);
                                                        
                                                        % PUT SATELLITE
                                                        antenna_PCV(m(1)).sat_type = sat_type;
                                                        
                                                        % get DAZI
                                                        while (isempty(strfind(line,'DAZI')))
                                                            l = l + 1; line = atx_file{l};
                                                        end
                                                        antenna_PCV(m(1)).dazi=sscanf(line(1:8),'%f');
                                                        
                                                        % get ZEN1 / ZEN2 / DZEN
                                                        while (isempty(strfind(line,'ZEN1 / ZEN2 / DZEN')))
                                                            l = l + 1; line = atx_file{l};
                                                        end
                                                        antenna_PCV(m(1)).zen1 = sscanf(line(1:8),'%f');
                                                        antenna_PCV(m(1)).zen2 = sscanf(line(9:14),'%f');
                                                        antenna_PCV(m(1)).dzen = sscanf(line(15:20),'%f');
                                                        
                                                        % get FREQUENCIES
                                                        while (isempty(strfind(line,'# OF FREQUENCIES')))
                                                            l = l + 1; line = atx_file{l};
                                                        end
                                                        antenna_PCV(m(1)).n_frequency=sscanf(line(1:8),'%d');
                                                        antenna_PCV(m(1)).offset = zeros(1,3,antenna_PCV(m(1)).n_frequency);
                                                        
                                                        %get information of each frequency
                                                        frequencies_found = 0;
                                                        
                                                        while frequencies_found < antenna_PCV(m(1)).n_frequency
                                                            while (isempty(strfind(line,'START OF FREQUENCY')))
                                                                l = l + 1; line = atx_file{l};
                                                            end
                                                            frequencies_found=frequencies_found+1;
                                                            antenna_PCV(m(1)).frequency_name(frequencies_found,:) = sscanf(line(4:6),'%s');
                                                            antenna_PCV(m(1)).frequency(frequencies_found) = sscanf(line(6),'%d');
                                                            
                                                            switch sscanf(line(4),'%c')
                                                                case 'G'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 1;
                                                                case 'R'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 2;
                                                                case 'E'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 3;
                                                                case 'J'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 4;
                                                                case 'C'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 5;
                                                                case 'I'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 6;
                                                                case 'S'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 7;
                                                            end
                                                            antenna_PCV(m(1)).sysfreq(frequencies_found) = antenna_PCV(m(1)).sys(frequencies_found) * 10 + antenna_PCV(m(1)).frequency(frequencies_found);
                                                            
                                                            while (isempty(strfind(line,'NORTH / EAST / UP')))
                                                                l = l + 1; line = atx_file{l};
                                                            end
                                                            if (~isempty(validity_start)) %satellite antenna
                                                                antenna_PCV(m(1)).offset(1,1:3,frequencies_found) = [sscanf(line(1:10),'%f'),sscanf(line(11:20),'%f'),sscanf(line(21:30),'%f')].*1e-3; % N,E,U
                                                                if (frequencies_found == antenna_PCV(m(1)).n_frequency)
                                                                    antenna_PCV(m(1)).available = 1;
                                                                end
                                                            else
                                                                antenna_PCV(m(1)).offset(1,1:3,frequencies_found) = [sscanf(line(11:20),'%f'),sscanf(line(1:10),'%f'),sscanf(line(21:30),'%f')].*1e-3; %E,N,U
                                                                antenna_PCV(m(1)).available = 1;
                                                            end
                                                            
                                                            number_of_zenith=(antenna_PCV(m(1)).zen2-antenna_PCV(m(1)).zen1)/antenna_PCV(m(1)).dzen+1;
                                                            if antenna_PCV(m(1)).dazi~=0
                                                                number_of_azimuth=(360-0)/antenna_PCV(m(1)).dazi+1;
                                                            else
                                                                number_of_azimuth=0;
                                                            end
                                                            
                                                            % NOAZI LINE
                                                            l = l + 1; line = atx_file{l};
                                                            antenna_PCV(m(1)).tableNOAZI(1,:,frequencies_found)=sscanf(line(9:end),'%f')'.*1e-3;
                                                            antenna_PCV(m(1)).tablePCV_zen(1,1:number_of_zenith,1)=antenna_PCV(m(1)).zen1:antenna_PCV(m(1)).dzen:antenna_PCV(m(1)).zen2;
                                                            
                                                            % TABLE AZI/ZEN DEPENDENT
                                                            if number_of_azimuth ~= 0
                                                                antenna_PCV(m(1)).tablePCV_azi(1,1:number_of_azimuth,1)=NaN(number_of_azimuth,1);
                                                                antenna_PCV(m(1)).tablePCV(:,:,frequencies_found)=NaN(number_of_azimuth,number_of_zenith);
                                                            else
                                                                antenna_PCV(m(1)).tablePCV_azi(1,1:number_of_azimuth,1)=NaN(1,1);
                                                                antenna_PCV(m(1)).tablePCV(:,:,frequencies_found)=NaN(1,number_of_zenith);
                                                            end
                                                            
                                                            l = l + 1; line = atx_file{l};
                                                            if (isempty(strfind(line,'END OF FREQUENCY')))
                                                                tablePCV=zeros(number_of_azimuth,number_of_zenith);
                                                                for i=1:number_of_azimuth
                                                                    tablePCV(i,:)=sscanf(line(9:end),'%f')'.*1e-3;
                                                                    l = l + 1; line = atx_file{l};
                                                                end
                                                                antenna_PCV(m(1)).tablePCV(:,:,frequencies_found)=tablePCV;
                                                                antenna_PCV(m(1)).tablePCV_azi(:,1:number_of_azimuth,1)=0:antenna_PCV(m(1)).dazi:360;
                                                            end
                                                            if number_of_azimuth == 0
                                                                antenna_PCV(m(1)).tablePCV(:,:,frequencies_found)=NaN(1,number_of_zenith);
                                                            end
                                                        end
                                                        found = found + length(m);
                                                        antenna_found(m) = 1;
                                                        for a = 2 : length(m)
                                                            antenna_PCV(m(a)) = antenna_PCV(m(1));
                                                        end
                                                    else % invalid_date
                                                        while (isempty(strfind(line,'END OF ANTENNA')))
                                                            l = l + 1; line = atx_file{l};
                                                        end
                                                    end
                                                elseif (nargin > 2) && strcmp(line(41:44),'    ')
                                                    flag_stop = true;
                                                    log.addMessage('There are no more antenna!!!',100);
                                                end
                                            end
                                        end
                                    end
                                case 2
                                    
                                    
                                case 3
                                    
                                    
                                case 0
                                    
                            end
                        else
                            log.addWarning('PCO/PCV file not loaded.\n');
                        end
                    else
                        log.addWarning('PCO/PCV file not loaded.\n');
                    end
                end
            end
            
            idx_not_found = find(~antenna_found);            
            if ~isempty(idx_not_found)
                w_msg = '';
                if numel(idx_not_found) == 1 && isempty(cell2mat(antmod(idx_not_found)))
                    w_msg = sprintf('No antenna model found for this station');
                else                   
                    no_ant = true;
                    for a = 1 : length(idx_not_found)
                        if ~isempty(antmod) && (length(antmod{idx_not_found(a)}) >=4 && ~strcmp(antmod{idx_not_found(a)}(1:4), 'NONE'))
                            no_ant = false;
                        end
                    end
                    if ~no_ant
                        w_msg = sprintf('The PCO/PCV model for the following antennas has not been found,\nsome models are missing or not defined at the time of processing');
                        for a = 1 : length(idx_not_found)
                            if ~isempty(antmod{idx_not_found(a)}) && ~strcmp(antmod{idx_not_found(a)}(1:4), 'NONE')
                                w_msg = sprintf('%s\n -  antenna model for "%s" is missing', w_msg, cell2mat(antmod(idx_not_found(a))));
                            end
                        end
                    end
                end
                if ~isempty(w_msg)
                    log.addWarning(w_msg);
                end
            end
        end
        
        function [status] = downloadHttpTxtResUncompress(filename, out_dir)
            log = Core.getLogger();
            fnp = File_Name_Processor();
            try
                options = weboptions;
                options.ContentType = 'text';
                options.Timeout = 15;
                [remote_location, filename, ext] = fileparts(filename);
                filename = [filename ext];
                log.addMessage(log.indent(sprintf('downloading %s ...',filename)));
                compressed_name = '';
                status = true;
                if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                    mkdir(out_dir);
                end
                try
                    txt = websave(fullfile(out_dir, filename), ['http://' remote_location '/' filename]);
                catch ex
                    if instr(ex.message, '404')
                        try
                            compressed_name = [filename, '.gz'];
                            txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name]);
                        catch ex
                            if instr(ex.message, '404')
                                try
                                    compressed_name = [filename, '.Z'];
                                    txt = websave(fullfile(out_dir, compressed_name), ['http://' remote_location '/' compressed_name]);
                                catch
                                    status = false;
                                end
                            end
                        end
                    end
                end
                if status
                    status = false; %#ok<NASGU>
                    if ~isempty(compressed_name)
                        compressed_name = fnp.checkPath(fullfile(out_dir, compressed_name));
                        if (isunix())
                            system(['gzip -d -f ' compressed_name '&> /dev/null &']);
                        else
                            try
                                [status, result] = system(['.\utility\thirdParty\7z1602-extra\7za.exe -y x '  compressed_name ' -o'  out_dir ]); %#ok<ASGLU>
                                if (status == 0)
                                    status = true;
                                end
                                delete(compressed_name);
                            catch
                                this.log.addError(sprintf('Please decompress the %s file before trying to use it in goGPS!!!', compressed_name));
                                status = false;
                            end
                        end
                    end
                    status = true;
                    log.addMessage(' Done');
                end
            catch
                status = false;
            end
        end
        
        function [status] = checkHttpTxtRes(filename)
            if isunix() || ismac()
                [resp, txt] = system(['curl --head ' filename]);
                if strfind(txt,'HTTP/1.1 200 OK')
                    status = true;
                else
                    [resp, txt] = system(['curl --head ' filename '.gz']);
                    if strfind(txt,'HTTP/1.1 200 OK')
                        status = true;
                    else
                        [resp, txt] = system(['curl --head ' filename '.Z']);
                        if strfind(txt,'HTTP/1.1 200 OK')
                            status = true;
                        else
                            status = false;
                        end
                    end
                end
            else
                log = Logger.getInstance.addWarning('HTTP check is implemeted only for Unix systems')
                status = true; % !!! to be implemented                
            end
        end
        
        function station_list = getStationList(dir_path, file_ext)
            % Get the list of stations present in a folder (with keys substituted)
            %
            % SYNTAX
            %   station_list = Core_Utilis.getStationList(dir_path)
            
            try
                % Calling dos is faster than dir with large directories
                if isunix
                    [~, d] = dos(['ls ' dir_path]); dir_list = strsplit(d);
                else
                    [~, d] = dos(['dir ' dir_path]); dir_list = strsplit(d);
                end
            catch
                dir_list = dir(dir_path);
                dir_list = {dir_list.name};
            end                                   
            
            % search for station files STAT${DOY}${S}${QQ}.${YY}
            if nargin == 1
                file_ext = '.';
            else
                file_ext = ['[' file_ext ']'];
            end
            file_list = [];
            for d = 1 : numel(dir_list)
                file_name_len = numel(dir_list{d});
                rin3_start = regexp(dir_list{d}, '\_R\_[0-9]{4}[0-9]{3}[0-9]{4}\_', 'once');
                if (file_name_len == 14) && ~isempty(regexp(dir_list{d}, ['.{4}[0-9]{3}.{1}[0-9]{2}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once'))
                    file_list = [file_list; [dir_list{d}(1:4) '${DOY}${S}${QQ}.${YY}' dir_list{d}(end)]]; %#ok<AGROW>
                    %file_list = [file_list; dir_list{d}(1:4)];
                elseif (file_name_len == 38) && ~isempty(rin3_start)
                    %file_list = [file_list; [dir_list{d}(1:rin3_start+2) '${YYYY}${DOY}${HH}${QQ}' dir_list{d}(rin3_start + 14 : end)]]; %#ok<AGROW>
                    file_list = [file_list; [dir_list{d}(1:rin3_start+2) '${YYYY}${DOY}' dir_list{d}(rin3_start + 10 : end)]]; %#ok<AGROW>
                end
            end
            station_list = {};
            if size(file_list, 2) > 1
                station_num = Core_Utils.code4Char2Num(file_list(:,1:4));
                station_name = unique(station_num);
                for s = 1 : numel(station_name)
                    station_list = [station_list; {file_list(find(station_num == station_name(s), 1, 'first'),:)}]; %#ok<AGROW>
                end
            end
            
            % search for station files STAT${DOY}${S}.${YY}
            file_list = [];
            for d = 1 : numel(dir_list)
                file_name_len = numel(dir_list{d});
                if (file_name_len == 12) && ~isempty(regexp(dir_list{d}, ['.{4}[0-9]{3}.{1}[\.]{1}[0-9]{2}' file_ext '{1}'], 'once'))
                    file_list = [file_list; [dir_list{d}(1:4) '${DOY}${S}.${YY}' dir_list{d}(end)]]; %#ok<AGROW>
                    %file_list = [file_list; dir_list{d}(1:4)];
                elseif (file_name_len == 38) && ~isempty(rin3_start)
                    file_list = [file_list; [dir_list{d}(1:rin3_start+2) '${YYYY}${DOY}' dir_list{d}(rin3_start + 10 : end)]]; %#ok<AGROW>
                end
            end
            if size(file_list, 2) > 1
                station_num = Core_Utils.code4Char2Num(file_list(:,1:4));
                station_name = unique(station_num);
                for s = 1 : numel(station_name)
                    station_list = [station_list {file_list(find(station_num == station_name(s), 1, 'first'),:)}]; %#ok<AGROW>
                end
            end
        end
    
        function data = injectData(data1, data2, idx1, idx2, data_size)
            % isert data2 into data1 at the position definied by idx1 and idx2
            % idx1 - 1 is the last element of data1 to be putted before data2 (0  if none)
            % idx2 + 1 is the first element of data1 to be put after data2 (data1 length if none)
            %
            % SYNTAX
            %   data = Core_Utils.injectData(data1, data2, idx1, idx2)
            if nargin == 5
                % check data size
                [m, n] = size(data2);
                if not(m == data_size(1) && n == data_size(2))
                    data2 = nan(data_size(1), data_size(2));
                end
            end
            if ~isempty(data2)
                data = [data1(1 : idx1 - 1, :); data2; data1(idx2 + 1 : end, :)];
            else
                data = data1;
            end
        end
        
        function data = injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start,interpolate)
            % inject smoothed data
            % 
            % INPUT:
            %   data_lft : all data left
            %   data_right : all data right
            %   idx_smt1 : which data of data_lft are to be smoothed
            %   idx_smt2 : which data of data_rgt are to be smoothed
            %   time_1: time of the left data to be smoothed
            %   time_2: time of the right data to be smoothed
            %   id_stop: first epoch of time_1 that should not be kept
            %   id_start: first epoch of time_2 to keep
            %
            % SYNTAX:
            %   data = Core_Utils.injectSmtData(data_lft, data_rgt, idx_smt1, idx_smt2, time_1, time_2, id_start)
            if nargin <9
            interpolate = true;
            end
            data_tosmt_lft = data_lft(idx_smt1);
            data_tosmt_rgt = data_rgt(idx_smt2);
            % we use mat time, is easier and we do not need extreme precision
            time_1 = time_1.getMatlabTime();
            time_2 = time_2.getMatlabTime();
            [idx1, idx2, time_tot] = Core_Utils.intersectOrderedDouble(time_1, time_2, median([diff(time_1); diff(time_2)])/4); % 1/4 the rate tolerance
            
            mix_len = min(0.007, abs((time_2(1) - time_1(end)))/20); % <= empirically found
            w2 = 1 ./ (1 + exp(-((time_tot - mean(time_tot)) / mix_len))); % todo: scale to ensure [0 1]
            w1 = 1 - w2;
            n_out = size(time_tot);
            data1 = nan(n_out);
            data2 = nan(n_out);
            data1(idx1) = data_tosmt_lft;
            data2(idx2) = data_tosmt_rgt;
            %id_start = idx1(id_start);            
            %id_ko = ((isnan(data1) & (1 : n_out)' < id_start) |(isnan(data2) & (1 : n_out)' >= id_start)) & ~(isnan(data1) & isnan(data2)); %?? should be used -> yes beacuse time is injected deleting overlapping times
            id_keep = unique([idx1(1:(id_stop-1)); idx2(id_start:end)]);
            % Interpolate missing data
            if interpolate
                is_nan = find(isnan(data1));
                extr_lft = is_nan(time_tot(is_nan) <= min(time_1));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time_1));
                data1(extr_lft) = data_tosmt_lft(1);
                data1(extr_rgh) = data_tosmt_lft(end);
                if any(~isnan(data_tosmt_lft)) &&  any(isnan(data1))
                    data1 = simpleFill1D(data1, isnan(data1), 'linear');
                end
                is_nan = find(isnan(data2));
                extr_lft = is_nan(time_tot(is_nan) <= min(time_2));
                extr_rgh = is_nan(time_tot(is_nan) >= max(time_2));
                data2(extr_lft) = data_tosmt_rgt(1);
                data2(extr_rgh) = data_tosmt_rgt(end);
                if any(~isnan(data_tosmt_rgt)) && any(isnan(data2))
                    data2 = simpleFill1D(data2, isnan(data2), 'linear');
                end
            else
                data1(isnan(data1)) = data2(isnan(data1));
                data2(isnan(data2)) = data1(isnan(data2));
            end
            % Merge
            data = w1.*data1 + w2.*data2;
            %data(id_ko) = [];
            data = data(id_keep);
            data = [data_lft(~idx_smt1); data; data_rgt(~idx_smt2)];            
        end
        
        function [idx1, idx2, double_tot] = intersectOrderedDouble(double_1, double_2, threshold)
            % given two ordered double give the index of the two vector in the joint vector considering the threshold
            % 
            % SYNTAX
            % [idx1, idx2] = Core_Utils.intersectOrderedDouble(double_1, double_2, threshold)
            l1 = length(double_1);
            l2 = length(double_2);
            idx1 = zeros(l1,1);
            idx2 = zeros(l2,1);
            i = 1;
            j = 1;
            tot = 1;
            while  i <=  l1 && j <= l2
                if abs(double_1(i) - double_2(j)) < threshold
                    idx1(i) = tot;
                    idx2(j) = tot;
                    i = i + 1;
                    j = j + 1;
                    tot = tot + 1;
                elseif double_1(i) < double_2(j)
                    idx1(i) = tot;
                    i = i + 1;
                    tot = tot + 1;
                else
                    idx2(j) = tot;
                    j = j + 1;
                    tot = tot +1;
                end
            end
            if j > l2 && i <= l1
                idx_end = (i : l1) -l1 + tot;
                idx1(i : l1) = idx_end;
            elseif i > l1 && j <= l2
                idx_end = (j : l2) -l2 + tot;
                idx2(j : l2) = idx_end;
            end
            double_tot = zeros(max(max(idx1), max(idx2)), 1);
            double_tot(idx1) = double_1;
            double_tot(idx2) = double_2;
        end
        
        function [wl_cyle_out, frac_bias] = getFracBias(wl_cycle, weigth)
            % get the common frac bias between cycles
            % NOTE: very coarse/nobrain/empirical solution - > a simpler one should be found
            %
            % SYNTAX
            %   [wl_cyle, frac_bias] = Core_Utils.getFracBias(wl_cycle)
            if nargin < 2
                weigth = ones(size(wl_cycle));
            end
            frac_bias = zeros(1,3);
            wl_cycle_frac = zeros(size(wl_cycle,1),3);
            % get receiver wsb
            wl_cycle_frac(:,1) = zero2nan(wl_cycle) - floor(zero2nan(wl_cycle));
            wl_cycle_frac(:,2) = zero2nan(wl_cycle) - round(zero2nan(wl_cycle));
            wl_cycle_frac(:,3) = zero2nan(wl_cycle) - ceil(zero2nan(wl_cycle));
            frac_bias(1) = median(wl_cycle_frac(:,1),'omitnan');
            frac_bias(2) = median(wl_cycle_frac(:,2),'omitnan');
            frac_bias(3) = median(wl_cycle_frac(:,3),'omitnan');
            wl_cyle_var = mean(abs(wl_cycle_frac-repmat(frac_bias,size(wl_cycle_frac,1),1)),'omitnan');
            [~,idx] = min(wl_cyle_var);
            frac_bias = frac_bias(idx);
            wl_cycle_frac = wl_cycle_frac(:,idx);
            a = 0;
            idx_rw = ones(size(wl_cycle_frac));
            while sum(zero2nan(wl_cycle_frac - frac_bias) < -0.5) > 0  || a < 4
                wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) = wl_cycle_frac((wl_cycle_frac-frac_bias) < -0.5) + 1;
                idx_rw(:) = 1;
                e = abs((wl_cycle_frac-frac_bias)) > 0.2;
                idx_rw(e) = 1./((wl_cycle_frac(e)-frac_bias)/0.2).^2; %idx_reweight
                idx_rw = idx_rw  .* weigth;
                idx_rw = idx_rw / sum(idx_rw);
                frac_bias = sum((wl_cycle_frac).*idx_rw);
                a = a+1;
            end
            wl_cyle_out = wl_cycle - frac_bias;
        end
        
        function [response] = timeIntersect(time1_st, time1_end, time2_st, time2_end)
            % check whether times of time bound 1 intersect with times of time bound 2
            %
            % SYNTAX
            % [response] = Core_Utils.timeIntersect(time1_st,time1_en, time2_st, time2_en)
            response = time1_st <= time2_end & time1_end >= time2_st;
            %(time2_st <= time1_st & time2_end >= time1_st) | (time2_st <= time1_en & time2_end >= time1_en);
        end
        
        function y_out = interp1LS(x_in, y_in, degree, x_out)
            % Least squares interpolant of a 1D dataset
            %
            % SYNTAX
            %   y_out = interp1LS(x_in, y_in, degree, x_out)
            
            if nargin < 4
                x_out = x_in;
            end
            
            for c = 1 : iif(min(size(y_in)) == 1, 1, size(y_in,2))
                if size(y_in, 1) == 1
                    y_tmp = y_in';
                else
                    y_tmp = y_in(:, c);
                end
                id_ok = ~isnan(y_tmp(:)) & ~isnan(x_in(:));
                x_tmp = x_in(id_ok);
                y_tmp = y_tmp(id_ok);
                
                n_obs = numel(x_tmp);
                A = zeros(n_obs, degree + 1);
                A(:, 1) = ones(n_obs, 1);
                for d = 1 : degree
                    A(:, d + 1) = x_tmp .^ d;
                end
                
                if (nargin < 4) && numel(x_out) == numel(x_tmp)
                    A2 = A;
                else
                    n_out = numel(x_out);
                    A2 = zeros(n_out, degree + 1);
                    A2(:, 1) = ones(n_out, 1);
                    for d = 1 : degree
                        A2(:, d + 1) = x_out .^ d;
                    end
                end
                
                warning('off')
                if min(size(y_in)) == 1
                    y_out = A2 * ((A' * A) \ (A' * y_tmp(:)));
                    y_out = reshape(y_out, size(x_out, 1), size(x_out, 2));
                else
                    y_out(:,c) = A2 * ((A' * A + eye(size(A,2))) \ (A' * y_tmp(:)));
                end
                warning('on')                
            end
        end
        
        function val = linInterpLatLonTime(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % Interpolate values froma data on a gepgraphical grid with multiple epoch
            % data structure: 
            %        first dimension : dlat (+) south pole -> north pole
            %        second dimension : dlon (+) west -> east
            %        third dimension : dr (+) time usual direction
            %        NOTE: dlat, dlon,dt do not have to be positive
            % 
            % INPUT:
            %      data - the data to be interpolate
            %      fist_lat - value of first lat value (max lat)
            %      dlat - px size lat
            %      first_lon - value of first lon value
            %      dlon - px size lon
            %      first_t - value of first time
            %      dt - px size time
            %      lat - lat at what we want to interpolate
            %      lon - lon at what we ant to interpolate
            %      gps_time - time at what we want to interpolate
            % NOTES 1 - all lat values should have same unit of measure
            %       2 - all lon values should have same unit of measure
            %       3 - all time values should have same unit of measure
            %       4 - the method will interpolate first in the dimesnion with less time
            % IMPORTANT : no double values at the borders should coexist: e.g. -180 180 or 0 360
            [nlat , nlon, nt] = size(data);
            n_in_lat = length(lat);
            n_in_lon = length(lon);
            n_in_t = length(t);
            assert(n_in_lat == n_in_lon);
            [ it, st, ilons, ilone, slon, ilat, slat] = Core_Utils.getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon,t);
            if n_in_lat > n_in_t % time first
                
                it = it*ones(size(ilat));
                % interpolate along time
                % [ 1 2  <= index of the cell at the smae time
                %   3 4]
                idx1 = sub2ind([nlat nlon nt], ilat, ilons, it);
                idx2 = sub2ind([nlat nlon nt], ilat, ilons, it+1);
                vallu = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt], ilat   , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat   , ilone , it+1);
                valru = data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilons , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilons , it+1);
                valld =  data(idx1).*(1-st) + data(idx2).*st;
                idx1 = sub2ind([nlat nlon nt],ilat+1 , ilone , it);
                idx2 = sub2ind([nlat nlon nt],ilat+1 , ilone , it+1);
                valrd =  data(idx1).*(1-st) + data(idx2).*st;
                
                %interpolate along long
                valu = vallu.*(1-slon) + valru.*slon;
                vald = valld.*(1-slon) + valrd.*slon;
                
                %interpolate along lat
                val = valu.*(1-slat) + vald.*slat;
                
            else %space first % NOTE: consider speed up in case only one time is present, unnecessary operations done
                % interpolate along lon
                valbu = permute(data(ilat   , ilons , it  ).*(1-slon) + data(ilat   , ilone , it  ).*slon,[3 1 2]);
                valau = permute(data(ilat   , ilons , min(it+1,size(data,3))).*(1-slon) + data(ilat   , ilone , min(it+1,size(data,3))).*slon,[3 1 2]);
                valbd = permute(data(ilat+1 , ilons , it  ).*(1-slon) + data(ilat+1 , ilone , it  ).*slon,[3 1 2]);
                valad = permute(data(ilat+1 , ilons , min(it+1,size(data,3))).*(1-slon) + data(ilat+1 , ilone , min(it+1,size(data,3))).*slon,[3 1 2]);
                
                %interpolate along lat
                valb = valbd.*(1-slat) + valbu.*slat;
                vala = valad.*(1-slat) + valau.*slat;
                
                %interpolate along time
                val = valb.*(1-st) + vala.*st;
            end
            
        end

        function [ it, st, ilons, ilone, slon, ilat, slat] = getIntIdx(data, first_lat, dlat, first_lon, dlon, first_t, dt, lat, lon, t)
            % get  interpolating index
            [nlat , nlon, nt] = size(data);
            
           
            lon(lon < first_lon) = lon(lon < first_lon) + nlon * dlon; %% to account for earth circularity 
            % find indexes and interpolating length
            % time
            it = max(min(floor((t - first_t)/ dt)+1,nt-1),1);
            st = max(min(t - first_t - (it-1)*dt, dt), 0) / dt;
            st = serialize(st);
            
            % lat
            ilat = max(min(floor((lat - first_lat)/ dlat)+1,nlat-1),1);
            slat = min(max(lat - first_lat - (ilat-1)*dlat, dlat), 0) / dlat;
            
            % lon
            ilons = max(min(floor((lon - first_lon)/ dlon)+1,nlon),1);
            ilone = ilons +1;
            ilone(ilone > nlon) = 1;
            slon = max(min(lon - first_lon- (ilons-1)*dlon, dlon), 0) / dlon;        
        end
            
        function response = permutedEqual(str1, str2)
            % check if the two variables are permuted version of the same sequence
            %
            % SYNTAX:
            %     response =Core_Utils.pertutedEqual(var1, var2)
            if length(str1) ~= length(str2)
                response = false;
            else
                ll = length(str1);
                found_all = true;
                i = 1;
                j = 1;
                while i <= ll
                    found = false;
                    while j <= ll && ~found
                        found = found || str1(i) == str2(j);
                        j = j+1;
                    end
                    found_all = found_all && found;
                    i = i+1;
                end
                response = found_all;
            end
            
        end
        
        function A = remBFromA(A,B)
            % reomve the lement of B from A
            for i = 1 : length(B)
                A(A==B(i)) = [];
            end
        end
        
        function [fb, frac_b_mat]= estimateFracBias(obs_cy, cycle_slip)
            % estimate the common factional bias to all the obesravtions
            %
            % SYNTAX:
            %    fb = Network.estimateFracBias(obs_cy, cycle_slip)
            amb_idx = Core_Utils.getAmbIdx(cycle_slip, obs_cy);
            frac_cy = obs_cy;
            n_arcs = max(amb_idx(:,end));
            frac_b = zeros(n_arcs,1);
            num_ep = zeros(n_arcs,1);
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb),'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = frac_cy(idx_amb) - amb;
                frac_b(i) = median(frac_cy(idx_amb),'omitnan'); %strongMean(frac_cy(idx_amb),1, 0.95,2.5);
                num_ep(i) = sum(sum(idx_amb));
                frac_b_mat(idx_amb) = frac_b(i);
            end
            frac_b_mat = frac_b_mat+0*obs_cy;
            fb = Core_Utils.circularModedRobustMean(frac_b_mat(:));
            frac_cy = nan(size(obs_cy));
            frac_b_mat = nan(size(obs_cy));
            for i = 1 : n_arcs
                idx_amb = amb_idx == i;
                amb = round(median(obs_cy(idx_amb) - fb,'omitnan'));%floor(strongMean(obs_cy(idx_amb),1, 0.95,2.5));
                frac_cy(idx_amb) = obs_cy(idx_amb) - amb - fb;
                frac_b_mat(idx_amb) = amb;
            end
            fb = fb + strongMean(frac_cy(:),1, 0.90,2.5);
        end
        
        function fr_cy = circularMean(cycle, obs_weigth)
            % estimate the mean for data over 0 -1
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularMean(cycle)
            
            % elimnate the nan
            if nargin < 2
                obs_weigth = ones(size(cycle));
            end
            idx_nan = isnan(cycle);
            cycle(idx_nan) = [];
            obs_weigth(idx_nan) = [];
            cycle = cycle*2*pi;
            % compute the unit vectors, take the mean and recompute the
            % average
            
            unit = [cos(cycle) sin(cycle)];
            obs_weigth = obs_weigth./ sum(obs_weigth);
            mean_vec = [mean(unit(:,1) .* obs_weigth) mean(unit(:,2) .* obs_weigth)];
            fr_cy = atan2(mean_vec(1), mean_vec(2))/(2*pi);
            
            
        end
        
        function fr_cy = circularModedRobustMean(cycle)
            % estimate a roubust mean mean for data over 0 -1 
            %
            % SYNTAX:
            %     fr_cy = Core_Utils.circularModedRobustMean(cycle)
            
            % elimnate the nan
            cycle(isnan(cycle)) = [];
            mode_cy = mode(cycle);
             
            % center around the mode and then take the string(robust) mean
            idx_inf = (cycle - mode_cy) < -0.5;
            idx_sup = (cycle - mode_cy) > 0.5;
            cycle(idx_inf) = cycle(idx_inf) +0.5;
            cycle(idx_sup) = cycle(idx_sup) -0.5;
            
            fr_cy = strongMean(cycle,1,0.95,2.5);
            
        end
        
        function amb_idx = getAmbIdx(cycle_slip , obs)
            % get matrix of same dimesion of the observation showing the ambiguity index of the obsarvation
            %
            % SYNTAX:
            % this.getAmbIdx()
            
            amb_idx = ones(size(cycle_slip));
            n_epochs = size(amb_idx,1);
            n_stream = size(amb_idx,2);
            for s = 1:n_stream
                if s > 1
                    amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                end
                cs = find(cycle_slip(:, s) > 0)';
                for c = cs
                    amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                end
            end
            amb_idx = zero2nan(amb_idx .* (obs ~= 0));
            amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx);
        end
        
        function createEmptyProject(base_dir, prj_name)
            % create empty config file
            %
            % SYNTAX
            %    createEmptyProject(base_dir, prj_name)
            %    createEmptyProject(prj_name)
            
            fnp = File_Name_Processor();
            state = Main_Settings('');

            if nargin == 1
                prj_name = base_dir;
                base_dir = fnp.getFullDirPath([state.getHomeDir filesep '..']);
            end
            
            log = Core.getLogger();
            log.addMarkedMessage(sprintf('Creating a new project "%s" into %s', prj_name, [base_dir filesep prj_name]));
            
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name]));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'config']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'out']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'RINEX']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station/CRD']));
            [status, msg, msgID] = mkdir(fnp.checkPath([base_dir filesep prj_name filesep 'station/MET']));
            state.setPrjHome(fnp.checkPath([base_dir filesep prj_name]));
            state.prj_name = prj_name;
            config_path = fnp.checkPath([base_dir filesep prj_name filesep 'config' filesep 'config.ini']);
            state.save(config_path);
            Global_Configuration.getCurrentSettings.import(state);
        end
        
        function y = fillNan1D(y,x)
            % fill the nan into the y
            %
            % SYNTAX
            % y = Core_Utils.fillNan1D(y,<x>)
            
            if nargin < x
                x = 1: length(y);
            end
            idx_nan = isnan(y);
            int_data = interp1(x(~idx_nan),y(~idx_nan),x(idx_nan));
            y(idx_nan) = int_data;
        end
        
        function [Amp,Phase,f] = compute_spectrum(y,smpl_rate)
            % compute the spectrum with fft
            %
            % SYNTAX:
            %  [Amp,Phase,f] = Core_Utils.compute_spectrum(y,smpl_rate);
            Y = fft(y);
            
            L = length(y);
            Fs = 1 /smpl_rate;
            f = Fs*(0:(L/2))/L;
            Y = Y(1:(L/2 +1));
            Amp = abs(Y/L)*2;
            Phase = angle(Y);
            
        end
        
        function [val] = spline(t,order)
            % Compute matrix entry for spline 
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            switch order
                case 1
                    val = Core_Utils.linearSpline(t);
                case 2
                    %%% tBD
                case 3
                    val = Core_Utils.cubicSpline(t);
            end
            
        end
        
        function [val] = cubicSpline(t)
            % Compute matrix entry for cubic spline 
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),4);
            val(:,1) = (1 - t).^3/6;
            val(:,2) = ((2-t).^3 - 4*(1-t).^3)/6;
            val(:,3) = ((1+t).^3 - 4*(t).^3)/6;
            val(:,4) = (t).^3/6;
        end
        
        function [val] = linearSpline(t)
            % Compute matrix entry for linear spline
            %
            % INPUT
            %   t -> 0 : 1
            %   order -> 1,3
            %
            % SYNTAX:
            %  Core_Utils.cubicSplic(t)
            val = zeros(numel(t),2);
            val(:,1) = 1 -t;
            val(:,2) = t;
        end
         
        function [idx, val] = hemisphereCubicSpline(n_az, n_el, az, el)
            % give the index of the hemisphere spline idx
            % first the equator then the first parallel then the second
            % parallel
            %
            % SYNTAX:
            %  [idx] = hemisphereSpline(n_az,n_el,az,el)
            el_step = 90/n_el;
            idx_el = repmat(ceil(el / el_step),1,4);
            idx_el(:,2) = idx_el(:,2) + 1;
            idx_el(:,3) = idx_el(:,3) + 2;
            idx_el(:,4) = idx_el(:,4) + 3;
            az_step = 360/n_el;
            idx_az = repmat(ceil(az / az_step),1,4);
            idx_az(:,2) = idx_az(:,2) + 1;
            idx_az(:,3) = idx_az(:,3) + 2;
            idx_az(:,4) = idx_az(:,4) + 3;
            idx_az(idx_az > n_az) = idx_az(idx_az > n_az) - n_az; 
            idx = idx_az + (idx_el -1).*n_el;
            idx = [idx idx+1 idx+2 idx+3];
            
            t_el = rem(el/el_step);
            t_az = rem(az/az_step);
            val_el = Core_Utils.cubicSpline(t_el);
            val_az = Core_Utils.cubicSpline(t_az);
            
            val = [val_az.*repmat(val_el(:,1),1,4) val_az.*repmat(val_el(:,2),1,4) val_az.*repmat(val_el(:,3),1,4) val_az.*repmat(val_el(:,4),1,4)];
            
            
        end
    end
end
