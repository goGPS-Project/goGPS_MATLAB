     
classdef Core_Utils < handle
    properties
    end
    methods
    function [antenna_PCV] = read_antenna_pcv(this, filename, antmod, date_limits)
            % SYNTAX:
            %   [antPCV] = this.read_antenna_PCV(filename, antmod, date_start, date_stop);
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

            for m = numel(antmod) : -1 : 1
                antenna_PCV(m) = struct('name', antmod{m}, ...
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
                                                id_ant = strfind(ant_char,line(1:20));
                                            else
                                                id_ant = strfind(ant_char, line(21:23));
                                            end
                                            if ~isempty(id_ant)
                                                if (nargin == 2) % receiver
                                                    m = (id_ant - 1) / 20 + 1; % I'm reading the antenna
                                                else
                                                    m = (id_ant - 1) / 3 + 1; % I'm reading the antenna
                                                end

                                                if ~(antenna_PCV(m(1)).available)

                                                    for a = 1:length(m)
                                                        this.log.addMessage(sprintf('Reading antenna %d => %s', m(a), antmod{m(a)}),100);
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
                                                                this.log.addMessage(sprintf(' - out of range -> (%s : %s) not after %s', date_limits.first.toString(), date_limits.last.toString(), datestr(validity_start)), 100)
                                                            else
                                                                this.log.addMessage(sprintf(' - out of range -> (%s : %s) not intersecting (%s : %s)', date_limits.first.toString(), date_limits.last.toString(), datestr(validity_start), datestr(validity_end)), 100)
                                                            end
                                                        end
                                                    else  %receiver antenna
                                                    end

                                                    if ~(invalid_date) % continue parsing
                                                        for a = 1:length(m)
                                                            this.log.addMessage(sprintf('Found a valid antenna %s', antmod{m(a)}), 50);
                                                        end
                                                        l = l_start;

                                                        % get TYPE
                                                        antenna_PCV(m(1)).type = line(1:20);

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
                                                            antenna_PCV(m(1)).frequency_name(frequencies_found,:)=sscanf(line(4:6),'%s');
                                                            antenna_PCV(m(1)).frequency(frequencies_found)=sscanf(line(6),'%d');

                                                            switch sscanf(line(4),'%c')
                                                                case 'G'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 1;
                                                                case 'R'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 2;
                                                                case 'E'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 3;
                                                                case 'C'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 4;
                                                                case 'J'
                                                                    antenna_PCV(m(1)).sys(frequencies_found) = 5;
                                                            end
                                                            antenna_PCV(m(1)).sysfreq(frequencies_found)=antenna_PCV(m(1)).sys(frequencies_found)*10+antenna_PCV(m(1)).frequency(frequencies_found);

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
                                                    this.log.addMessage('There are no more antenna!!!',100);
                                                end
                                            end
                                        end
                                    end
                                case 2


                                case 3


                                case 0

                            end
                        else
                            this.log.addWarning('PCO/PCV file not loaded.\n');
                        end
                    else
                        this.log.addWarning('PCO/PCV file not loaded.\n');
                    end
                end
            end

            idx_not_found = find(~antenna_found);
            if ~isempty(idx_not_found)
                w_msg = sprintf('The PCO/PCV model for the following antennas has not been found,\nsome models are missing or not defined at the time of processing');
                for a = 1 : length(idx_not_found)
                    w_msg = sprintf('%s\n -  antenna model for "%s" is missing', w_msg, cell2mat(antmod(idx_not_found(a))));
                end
                this.log.addWarning(w_msg);
            end
        end
    end
end
