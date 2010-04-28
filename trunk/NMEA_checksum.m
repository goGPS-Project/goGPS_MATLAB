function checksum = NMEA_checksum(nmeastring)

% SYNTAX:
%   checksum = NMEA_checksum(nmeastring);
%
% INPUT:
%   nmeastring = NMEA sentence without checksum
%
% OUTPUT:
%   checksum = checksum to be appended to NMEA sentence
%
% DESCRIPTION:
%   Checksum computation as required by NMEA 0183 format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% (code publicly available at
% http://www.mathworks.com/matlabcentral/fileexchange/15080)
%
%----------------------------------------------------------------------------------------------

checksum = 0;

% see if string contains the * which starts the checksum and keep string
% upto * for generating checksum
nmeastring = strtok(nmeastring,'*');

nmeastring_d = double(nmeastring);                    % convert characters in string to double values
for count = 2:length(nmeastring)                      % checksum computation ignores $ at start
    checksum = bitxor(checksum,nmeastring_d(count));  % checksum computation
    checksum = uint16(checksum);                      % make sure that checksum is unsigned int16
end

% convert checksum to hex value
checksum = double(checksum);
checksum = dec2hex(checksum);

% add leading zero to checksum if it is a single digit, e.g. 4 has a 0
% added so that the checksum is 04
[null, nchar] = size(checksum); %#ok<ASGLU>
if (nchar == 1)
    [checksum] = sprintf('0%s',checksum);
end
