function [data] = decode_4Ah(msg)

% SYNTAX:
%   [data] = decode_4Ah(msg);
%
% INPUT:
%   msg = message transmitted by the NVS receiver
%
% OUTPUT:
%   data = cell-array that contains the 4Ah packet information
%          1.1) message id (4Ah)
%          2.1) ionosphere parameter (alpha 0)
%          2.2) ionosphere parameter (alpha 1)
%          2.3) ionosphere parameter (alpha 2)
%          2.4) ionosphere parameter (alpha 3)
%          2.5) ionosphere parameter (beta 0)
%          2.6) ionosphere parameter (beta 1)
%          2.7) ionosphere parameter (beta 2)
%          2.8) ionosphere parameter (beta 3)
%
% DESCRIPTION:
%   BINR 4Ah binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Daisuke Yoshida
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
data{2} = zeros(8,1);

% output data save
data{1} = '4Ah';

%------------------------------------------------

%check the minimum allowed length of the message
if (length(msg) < 256)
    return
end

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

param = [alpha0 alpha1 alpha2 alpha3 beta0 beta1 beta2 beta3];

if (any(param-param(1)))
    data{2}(1) = alpha0;
    data{2}(2) = alpha1;
    data{2}(3) = alpha2;
    data{2}(4) = alpha3;
    data{2}(5) = beta0;
    data{2}(6) = beta1;
    data{2}(7) = beta2;
    data{2}(8) = beta3;
end
