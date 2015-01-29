function [data] = decode_AID_HUI(msg)

% SYNTAX:
%   [data] = decode_AID_HUI(msg);
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
%
% OUTPUT:
%   data = cell-array that contains the AID-HUI packet information
%          1.1) message class-id (AID-HUI)
%          2.1) satellite health
%          3.1) UTC parameter A1
%          3.2) UTC parameter A0
%          3.3) UTC reference Time-Of-Week
%          3.4) UTC reference Week Number
%          3.5) UTC leap seconds before event
%          3.6) UTC Week Number when next leap second event occurs
%          3.7) UTC Day-Of-Week when next leap second event occurs
%          3.8) UTC leap seconds after event
%          3.9)  ionosphere parameter (alpha 0)
%          3.10) ionosphere parameter (alpha 1)
%          3.11) ionosphere parameter (alpha 2)
%          3.12) ionosphere parameter (alpha 3)
%          3.13) ionosphere parameter (beta 0)
%          3.14) ionosphere parameter (beta 1)
%          3.15) ionosphere parameter (beta 2)
%          3.16) ionosphere parameter (beta 3)
%
% DESCRIPTION:
%   AID-HUI binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

% first message initial index
pos = 1;

% output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(32,1);
data{3} = zeros(16,1);

% output data save
data{1} = 'AID-HUI';

% satellite health bitmask (4 bytes: bit set = sat healthy)
health = msg(pos:pos+31);  pos = pos + 32;

for i = 1 : 32
    data{2}(i) = str2double(health(i));
end

%------------------------------------------------

% UTC parameter A1 (8 bytes)
A1field = msg(pos:pos+63);
pos = pos + 64;

% byte order inversion (little endian)
A1field = fliplr(reshape(A1field,8,[]));
A1field = A1field(:)';

% floating point value decoding (double floating point)
sign = str2num(A1field(1));
esp  = fbin2dec(A1field(2:12));
mant = fbin2dec(A1field(13:64)) / 2^52;
A1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

%------------------------------------------------

% UTC parameter A0 (8 bytes)
A0field = msg(pos:pos+63);
pos = pos + 64;

% byte order inversion (little endian)
A0field = fliplr(reshape(A0field,8,[]));
A0field = A0field(:)';

% floating point value decoding (double floating point)
sign = str2num(A0field(1));
esp  = fbin2dec(A0field(2:12));
mant = fbin2dec(A0field(13:64)) / 2^52;
A0 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

%------------------------------------------------

% UTC reference Time-Of-Week (4 bytes)
TOW1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW3 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW4 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW = TOW1 + (TOW2 * 2^8) + (TOW3 * 2^16) + (TOW4 * 2^24);  % little endian
% TOW = TOW / 1000;

%------------------------------------------------

% UTC Week Number (2 bytes)
WEEK1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK = WEEK1 + (WEEK2 * 2^8);        % little endian

%------------------------------------------------

% UTC Leap seconds (2 bytes)
LS1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
LS2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
LS = LS1 + (LS2 * 2^8);        % little endian

%------------------------------------------------

% UTC Week Number when next leap seconds event occurs (2 bytes)
WNF1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WNF2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WNF = WNF1 + (WNF2 * 2^8);        % little endian

%------------------------------------------------

% UTC Day Of Week when next leap seconds event occurs (2 bytes)
DN1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
DN2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
DN = DN1 + (DN2 * 2^8);        % little endian

%------------------------------------------------

% UTC leap seconds after event (2 bytes)
LSF1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
LSF2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
LSF = LSF1 + (LSF2 * 2^8);        % little endian

%------------------------------------------------

% Spare slot (2 bytes)
pos = pos + 16;

%------------------------------------------------

% ionosphere parameter alpha 0 [s]
a0field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
a0field = fliplr(reshape(a0field,8,[]));
a0field = a0field(:)';

% floating point value decoding (single floating point)
sign = str2num(a0field(1));
esp  = fbin2dec(a0field(2:9));
mant = fbin2dec(a0field(10:32)) / 2^23;
alpha0 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter alpha 1 [s/semicircle]
a1field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
a1field = fliplr(reshape(a1field,8,[]));
a1field = a1field(:)';

% floating point value decoding (single floating point)
sign = str2num(a1field(1));
esp  = fbin2dec(a1field(2:9));
mant = fbin2dec(a1field(10:32)) / 2^23;
alpha1 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter alpha 2 [s/semicircle^2]
a2field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
a2field = fliplr(reshape(a2field,8,[]));
a2field = a2field(:)';

% floating point value decoding (single floating point)
sign = str2num(a2field(1));
esp  = fbin2dec(a2field(2:9));
mant = fbin2dec(a2field(10:32)) / 2^23;
alpha2 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter alpha 3 [s/semicircle^3]
a3field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
a3field = fliplr(reshape(a3field,8,[]));
a3field = a3field(:)';

% floating point value decoding (single floating point)
sign = str2num(a3field(1));
esp  = fbin2dec(a3field(2:9));
mant = fbin2dec(a3field(10:32)) / 2^23;
alpha3 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter beta 0 [s]
b0field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
b0field = fliplr(reshape(b0field,8,[]));
b0field = b0field(:)';

% floating point value decoding (single floating point)
sign = str2num(b0field(1));
esp  = fbin2dec(b0field(2:9));
mant = fbin2dec(b0field(10:32)) / 2^23;
beta0 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter beta 1 [s]
b1field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
b1field = fliplr(reshape(b1field,8,[]));
b1field = b1field(:)';

% floating point value decoding (single floating point)
sign = str2num(b1field(1));
esp  = fbin2dec(b1field(2:9));
mant = fbin2dec(b1field(10:32)) / 2^23;
beta1 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter beta 2 [s]
b2field = msg(pos:pos+31);
pos = pos + 32;

% byte order inversion (little endian)
b2field = fliplr(reshape(b2field,8,[]));
b2field = b2field(:)';

% floating point value decoding (single floating point)
sign = str2num(b2field(1));
esp  = fbin2dec(b2field(2:9));
mant = fbin2dec(b2field(10:32)) / 2^23;
beta2 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

% ionosphere parameter beta 3 [s]
b3field = msg(pos:pos+31);
pos = pos + 32; %#ok<NASGU>

% byte order inversion (little endian)
b3field = fliplr(reshape(b3field,8,[]));
b3field = b3field(:)';

% floating point value decoding (single floating point)
sign = str2num(b3field(1));
esp  = fbin2dec(b3field(2:9));
mant = fbin2dec(b3field(10:32)) / 2^23;
beta3 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

%------------------------------------------------

%data output save
data{3}(1) = A1;
data{3}(2) = A0;
data{3}(3) = TOW;
data{3}(4) = WEEK;
data{3}(5) = LS;
data{3}(6) = WNF;
data{3}(7) = DN;
data{3}(8) = LSF;
data{3}(9)  = alpha0;
data{3}(10) = alpha1;
data{3}(11) = alpha2;
data{3}(12) = alpha3;
data{3}(13) = beta0;
data{3}(14) = beta1;
data{3}(15) = beta2;
data{3}(16) = beta3;
