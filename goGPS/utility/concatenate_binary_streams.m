function concatenate_binary_streams(filerootIN, filenameOUT, hour_idx)

% SYNTAX:
%   concatenate_binary_streams(filerootIN, filenameOUT, hour_idx);
%
% INPUT:
%   filerootIN  = prefix of input files (e.g. '../data/VRS_20140819_master')
%   filenameOUT = output filename (e.g. '../data/VRS_20140819.rtcm')
%   hour_idx    = starting hour index (e.g. for a file list starting with
%                 *_rover_036.bin --> hour_idx = 36)
%
% DESCRIPTION:
%   Function that concatenate hourly binary stream files produced by goGPS
%   monitors into one file.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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

if (nargin < 3)
    hour = 0;                                                       %hour index (integer)
else
    hour = hour_idx;
end
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([filerootIN '_' hour_str '.bin']);                          %first input file
fid_out = fopen(filenameOUT,'w+');                                  %output file opening
while ~isempty(d)
    fprintf('%s\n',['Reading: ' filerootIN '_' hour_str '.bin']);
    num_bytes = d.bytes;                                            %input file size (number of bytes)
    fid_in = fopen([filerootIN '_' hour_str '.bin'],'r+');          %input file opening
    buf = fread(fid_in,num_bytes,'uint8');                          %input file reading
    fclose(fid_in);                                                 %input file closing
    fwrite(fid_out,buf,'uint8');                                    %output file writing
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
    d = dir([filerootIN '_' hour_str '.bin']);                      %file to be read
end
fclose(fid_out);
