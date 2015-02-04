function [antenna_PCV] = read_antenna_PCV(filename, antmod)

% SYNTAX:
%   [antPCV] = read_antenna_PCV(filename, antmod);
%
% INPUT:
%   filename = antenna phase center offset/variation file
%   antmod = cell-array containing antenna model strings
%
% OUTPUT:
%   antenna_PCV (see description below)
%
% DESCRIPTION:
%   Extracts antenna phase center offset values from a PCO/PCV file in ATX format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
%
% Portions of code contributed by Stefano Caldera
%----------------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------------

% antenna_PCV struct definition
% antenna_PCV.name           : antenna name (with radome code)
% antenna_PCV.n_frequency    : number of available frequencies
% antenna_PCV.frequency_name : array with name of available frequencies ({'G01';'G02';'R01',...})
% antenna_PCV.frequency      : array with list of frequencies (carrier number) corresponding to the frequencies name ({'1';'2';'1',...})
% antenna_PCV.sys            : array with code id of the system constellation of each frequency (1: GPS, 2: GLONASS, ...)
% antenna_PCV.sysfreq        : array with codes of the system constellation and carrier of each frequency (11: GPS L1, 12: GPS L2, 21: GLONASS L1, ...)
% antenna_PCV.offset         : ENU offset (one array for each frequency)
% antenna_PCV.dazi           : increment of the azimuth (0.0 for non-azimuth-dependent phase center variations)
% antenna_PCV.zen1           : Definition of the grid in zenith angle: minimum zenith angle
% antenna_PCV.zen2           : Definition of the grid in zenith angle: maximum zenith angle
% antenna_PCV.dzen           : Definition of the grid in zenith angle: increment of the zenith angle
% antenna_PCV.tableNOAZI     : PCV values for NOAZI, in a cell array with a vector for each frequency [m]
% antenna_PCV.tablePCV       : PCV values elev/azim depentend, in a cell array with a matrix for each frequency [m]
% antenna_PCV.tablePCV_zen   : zenith angles corresponding to each column of antenna_PCV.tablePCV
% antenna_PCV.tablePCV_azi   : azimutal angles corresponding to each row of antenna_PCV.tablePCV


antenna_PCV=[];
for m = 1 : length(antmod)
    antenna_PCV(m).name=antmod{m};
    antenna_PCV(m).n_frequency=0;
end
antenna_found=zeros(length(antmod),1);


for file_pcv=1:size(filename,1)
    if sum(antenna_found)<length(antmod)
        if (~isempty(filename))
            fid = fopen(char(filename(file_pcv,:)),'r');
            if (fid ~= -1)
                found = 0;
                format = 0;
                % get format (1: ATX, 2: Bernese 5.0, 3: Bernese 5.2)
                line = fgetl(fid);
                if ~isempty(strfind(line,'ANTEX VERSION / SYST'))
                    format = 1;
                end
                if ~isempty(strfind(line,'MODEL NAME:'))
                    format = 2;
                end
                if ~isempty(strfind(line,'ANTENNA PHASE CENTER VARIATIONS DERIVED FROM ANTEX FILE'))
                    format = 3;
                end
                
                switch format
                    %% ATX
                    case 1
                        while (~feof(fid) && found < length(antmod))
                            %parse the next line
                            line = fgetl(fid);
                            
                            %check if any of the requested antenna models is in this line
                            for m = 1 : length(antmod)
                                %antenna_PCV(m).name=antmod{m};
                                %antenna_PCV(m).n_frequency=0;
                                %antenna_PCV(m).offset=[0 0 0];
                                
                                answer1 = strfind(line,antmod{m});
                                answer2 = strfind(line,'TYPE / SERIAL NO');
                                if (~isempty(answer1) && ~isempty(answer2) && ~isempty(find(antmod{m} ~= ' ', 1)))
                                    
                                    %get DAZI
                                    while (isempty(strfind(line,'DAZI')))
                                        line = fgetl(fid);
                                    end
                                    antenna_PCV(m).dazi=sscanf(line(1:8),'%f');
                                    
                                    %get ZEN1 / ZEN2 / DZEN
                                    while (isempty(strfind(line,'ZEN1 / ZEN2 / DZEN')))
                                        line = fgetl(fid);
                                    end
                                    antenna_PCV(m).zen1=sscanf(line(1:8),'%f');
                                    antenna_PCV(m).zen2=sscanf(line(9:14),'%f');
                                    antenna_PCV(m).dzen=sscanf(line(15:20),'%f');
                                    
                                    %get FREQUENCIES
                                    while (isempty(strfind(line,'# OF FREQUENCIES')))
                                        line = fgetl(fid);
                                    end
                                    antenna_PCV(m).n_frequency=sscanf(line(1:8),'%d');
                                    antenna_PCV(m).offset=zeros(1,3,antenna_PCV(m).n_frequency);
                                    
                                    %get information of each frequency
                                    frequencies_found=0;
                                    
                                    while frequencies_found<antenna_PCV(m).n_frequency
                                        while (isempty(strfind(line,'START OF FREQUENCY')))
                                            line = fgetl(fid);
                                        end
                                        frequencies_found=frequencies_found+1;
                                        antenna_PCV(m).frequency_name(frequencies_found,:)=sscanf(line(4:6),'%s');
                                        antenna_PCV(m).frequency(frequencies_found)=sscanf(line(6),'%d');
                                        
                                        switch sscanf(line(4),'%c')
                                            case 'G'
                                                antenna_PCV(m).sys(frequencies_found) = 1;
                                            case 'R'
                                                antenna_PCV(m).sys(frequencies_found) = 2;
                                            case 'E'
                                                antenna_PCV(m).sys(frequencies_found) = 3;
                                            case 'C'
                                                antenna_PCV(m).sys(frequencies_found) = 4;
                                            case 'J'
                                                antenna_PCV(m).sys(frequencies_found) = 5;
                                        end
                                        antenna_PCV(m).sysfreq(frequencies_found)=antenna_PCV(m).sys(frequencies_found)*10+antenna_PCV(m).frequency(frequencies_found);
                                        
                                        while (isempty(strfind(line,'NORTH / EAST / UP')))
                                            line = fgetl(fid);
                                        end
                                        
                                        antenna_PCV(m).offset(1,1:3,frequencies_found)=[sscanf(line(11:20),'%f'),sscanf(line(1:10),'%f'),sscanf(line(21:30),'%f')].*1e-3; %E,N,U
                                        
                                        number_of_zenith=(antenna_PCV(m).zen2-antenna_PCV(m).zen1)/antenna_PCV(m).dzen+1;
                                        if antenna_PCV(m).dazi~=0
                                            number_of_azimuth=(360-0)/antenna_PCV(m).dazi+1;
                                        else
                                            number_of_azimuth=0;
                                        end
                                            
                                        % NOAZI LINE
                                        line = fgetl(fid);
                                        antenna_PCV(m).tableNOAZI(1,:,frequencies_found)=sscanf(line(9:end),'%f')'.*1e-3;
                                        antenna_PCV(m).tablePCV_zen(1,1:number_of_zenith,1)=antenna_PCV(m).zen1:antenna_PCV(m).dzen:antenna_PCV(m).zen2;
                                        % TABLE AZI/ZEN DEPENDENT
                                        if number_of_azimuth ~= 0
                                            antenna_PCV(m).tablePCV_azi(1,1:number_of_azimuth,1)=NaN(number_of_azimuth,1);
                                            antenna_PCV(m).tablePCV(:,:,frequencies_found)=NaN(number_of_azimuth,number_of_zenith);
                                        else
                                            antenna_PCV(m).tablePCV_azi=NaN;
                                            antenna_PCV(m).tablePCV=NaN;                                          
                                        end
                                        
                                        line = fgetl(fid);
                                        if (isempty(strfind(line,'END OF FREQUENCY')))
                                            tablePCV=zeros(number_of_azimuth,number_of_zenith);
                                            for i=1:number_of_azimuth
                                                tablePCV(i,:)=sscanf(line(9:end),'%f')'.*1e-3;
                                                line = fgetl(fid);
                                            end
                                            antenna_PCV(m).tablePCV(:,:,frequencies_found)=tablePCV;
                                            antenna_PCV(m).tablePCV_azi(:,1:number_of_azimuth,1)=0:antenna_PCV(m).dazi:360;
                                        end
                                        if number_of_azimuth == 0
                                            antenna_PCV(m).tablePCV=NaN;
                                        end
                                    end
                                    found = found + 1;
                                end
                                for n = m+1 : length(antmod)
                                    if (strcmp(antmod{n}, antmod{m}))
                                        antenna_PCV(n) = antenna_PCV(m);
                                    end
                                end
                            end
                        end
                        
                    case 2
                        
                        
                    case 3
                        
                        
                    case 0
                        
                end
                
                fclose(fid);
            else
                fprintf('... WARNING: PCO/PCV file not loaded.\n');
            end
        else
            fprintf('... WARNING: PCO/PCV file not loaded.\n');
        end
    end
end



