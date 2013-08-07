function [data] = decode_F7h(msg)

% SYNTAX:
%   [data] = decode_F7h(msg);
%
% INPUT:
%   msg = message transmitted by the NVS receiver
%
% OUTPUT:
%   data = cell-array that contains the AID-EPH packet information
%          1.1) message id (F7h)
%          2.1) GPS satellite id
%          2.2) GPS af2
%          2.3) GPS M0
%          2.4) GPS root A
%          2.5) GPS delta-N
%          2.6) GPS eccentricity 
%          2.7) GPS omega 
%          2.8) GPS Cuc
%          2.9) GPS Cus
%          2.10)GPS Crc
%          2.11)GPS Crs
%          2.12)GPS i0
%          2.13)GPS IDOT 
%          2.14)GPS Cic
%          2.15)GPS Cis
%          2.16)GPS omega0
%          2.17)GPS omegadot
%          2.18)GPS toe
%          2.19)GPS af0
%          2.20)GPS af1
%          2.21)GPS toc
%          2.22)GPS IODE
%          2.23)GPS codes; 
%          2.24)GPS weekno;
%          2.25)GPS L2flag;
%          2.26)GPS svaccur;
%          2.27)GPS svhealth; -- not available
%          2.28)GPS tgd;
%          2.29)GPS fit_int; -- not available
%          2.30)multi-constellation satellite index (here only GPS is assumed)
%
% DESCRIPTION:
%   BINR F7h binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(33,1);

%output data save
data{1} = 'F7h';

%To check whether GPS(1) or GLONASS(2)
type = msg(pos:pos+7); pos = pos + 8;
type = fliplr(reshape(type,8,[]));                % byte order inversion (little endian)
type = type(:)';
type = fbin2dec((type(1:8)));

%satellite PRN (4 bytes)
PRN = msg(pos:pos+7); pos = pos + 8;
PRN = fliplr(reshape(PRN,8,[]));                  % byte order inversion (little endian)
PRN = PRN(:)';
PRN = fbin2dec((PRN(1:8)));

if (type == 1 && length(msg(pos:end)) >= 1088)

    %------------------------------------------------
    % Crs   
    Crs = msg(pos:pos+31); pos = pos + 32;
    Crs = fliplr(reshape(Crs,8,[]));              % byte order inversion (little endian)
    Crs = Crs(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Crs(1));
    esp  = fbin2dec(Crs(2:9));
    mant = fbin2dec(Crs(10:32)) / 2^23;
    Crs = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % delta_n 
    delta_n = msg(pos:pos+31); pos = pos + 32;
    delta_n = fliplr(reshape(delta_n,8,[]));      % byte order inversion (little endian)
    delta_n = delta_n(:)';

    % floating point value decoding (single floating point)
    sign = str2num(delta_n(1));
    esp  = fbin2dec(delta_n(2:9));
    mant = fbin2dec(delta_n(10:32)) / 2^23;
    delta_n = (-1)^sign * (2^(esp - 127)) * (1 + mant);
    delta_n = delta_n*1e3; %from rad/ms to rad/s

    %------------------------------------------------
    % M0      
    M0 = msg(pos:pos+63); pos = pos + 64;
    M0 = fliplr(reshape(M0,8,[]));                % byte order inversion (little endian)
    M0 = M0(:)';

    % floating point value decoding (double floating point)
    sign = str2num(M0(1));
    esp  = fbin2dec(M0(2:12));
    mant = fbin2dec(M0(13:64)) / 2^52;
    M0 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % Cuc     
    Cuc = msg(pos:pos+31); pos = pos + 32;
    Cuc = fliplr(reshape(Cuc,8,[]));              % byte order inversion (little endian)
    Cuc = Cuc(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Cuc(1));
    esp  = fbin2dec(Cuc(2:9));
    mant = fbin2dec(Cuc(10:32)) / 2^23;
    Cuc = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % e       
    e = msg(pos:pos+63); pos = pos + 64;
    e = fliplr(reshape(e,8,[]));                  % byte order inversion (little endian)
    e = e(:)';

    % floating point value decoding (double floating point)
    sign = str2num(e(1));
    esp  = fbin2dec(e(2:12));
    mant = fbin2dec(e(13:64)) / 2^52;
    e = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % Cus    
    Cus = msg(pos:pos+31); pos = pos + 32;
    Cus = fliplr(reshape(Cus,8,[]));              % byte order inversion (little endian)
    Cus = Cus(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Cus(1));
    esp  = fbin2dec(Cus(2:9));
    mant = fbin2dec(Cus(10:32)) / 2^23;
    Cus = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % root_A  
    root_A = msg(pos:pos+63); pos = pos + 64;
    root_A = fliplr(reshape(root_A,8,[]));        % byte order inversion (little endian)
    root_A = root_A(:)';

    % floating point value decoding (double floating point)
    sign = str2num(root_A(1));
    esp  = fbin2dec(root_A(2:12));
    mant = fbin2dec(root_A(13:64)) / 2^52;
    root_A = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % toe     
    toe = msg(pos:pos+63); pos = pos + 64;
    toe = fliplr(reshape(toe,8,[]));              % byte order inversion (little endian)
    toe = toe(:)';

    % floating point value decoding (double floating point)
    sign = str2num(toe(1));
    esp  = fbin2dec(toe(2:12));
    mant = fbin2dec(toe(13:64)) / 2^52;
    toe = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    toe = toe*1e-3; %from ms to s
  
    %------------------------------------------------
    % Cic   
    Cic = msg(pos:pos+31); pos = pos + 32;
    Cic = fliplr(reshape(Cic,8,[]));              % byte order inversion (little endian)
    Cic = Cic(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Cic(1));
    esp  = fbin2dec(Cic(2:9));
    mant = fbin2dec(Cic(10:32)) / 2^23;
    Cic = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % omega0   
    omega0 = msg(pos:pos+63); pos = pos + 64;
    omega0 = fliplr(reshape(omega0,8,[]));        % byte order inversion (little endian)
    omega0 = omega0(:)';

    % floating point value decoding (double floating point)
    sign = str2num(omega0(1));
    esp  = fbin2dec(omega0(2:12));
    mant = fbin2dec(omega0(13:64)) / 2^52;
    omega0 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % Cis      
    Cis = msg(pos:pos+31); pos = pos + 32;
    Cis = fliplr(reshape(Cis,8,[]));              % byte order inversion (little endian)
    Cis = Cis(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Cis(1));
    esp  = fbin2dec(Cis(2:9));
    mant = fbin2dec(Cis(10:32)) / 2^23;
    Cis = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % i0     
    i0 = msg(pos:pos+63); pos = pos + 64;
    i0 = fliplr(reshape(i0,8,[]));                 % byte order inversion (little endian)
    i0 = i0(:)';

    % floating point value decoding (double floating point)
    sign = str2num(i0(1));
    esp  = fbin2dec(i0(2:12));
    mant = fbin2dec(i0(13:64)) / 2^52;
    i0 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % Crc     
    Crc = msg(pos:pos+31); pos = pos + 32;
    Crc = fliplr(reshape(Crc,8,[]));                % byte order inversion (little endian)
    Crc = Crc(:)';

    % floating point value decoding (single floating point)
    sign = str2num(Crc(1));
    esp  = fbin2dec(Crc(2:9));
    mant = fbin2dec(Crc(10:32)) / 2^23;
    Crc = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    % W  
    omega = msg(pos:pos+63); pos = pos + 64;
    omega = fliplr(reshape(omega,8,[]));            % byte order inversion (little endian)
    omega = omega(:)';
    
    % floating point value decoding (double floating point)
    sign = str2num(omega(1));
    esp  = fbin2dec(omega(2:12));
    mant = fbin2dec(omega(13:64)) / 2^52;
    omega = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------
    % omegadot 
    omegadot = msg(pos:pos+63); pos = pos + 64;
    omegadot = fliplr(reshape(omegadot,8,[]));       % byte order inversion (little endian)
    omegadot = omegadot(:)';
    
    % floating point value decoding (double floating point)
    sign = str2num(omegadot(1));
    esp  = fbin2dec(omegadot(2:12));
    mant = fbin2dec(omegadot(13:64)) / 2^52;
    omegadot = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    omegadot = omegadot*1e3; %from rad/ms to rad/s

    %------------------------------------------------
    % IDOT     
    IDOT = msg(pos:pos+63); pos = pos + 64;
    IDOT = fliplr(reshape(IDOT,8,[]));                % byte order inversion (little endian)
    IDOT = IDOT(:)';
    
    % floating point value decoding (double floating point)
    sign = str2num(IDOT(1));
    esp  = fbin2dec(IDOT(2:12));
    mant = fbin2dec(IDOT(13:64)) / 2^52;
    IDOT = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    IDOT = IDOT*1e3; %from rad/ms to rad/s

    %------------------------------------------------
    % tgd   
    tgd = msg(pos:pos+31); pos = pos + 32;
    tgd = fliplr(reshape(tgd,8,[]));                  % byte order inversion (little endian)
    tgd = tgd(:)';

    % floating point value decoding (single floating point)
    sign = str2num(tgd(1));
    esp  = fbin2dec(tgd(2:9));
    mant = fbin2dec(tgd(10:32)) / 2^23;
    tgd = (-1)^sign * (2^(esp - 127)) * (1 + mant);
    tgd = tgd*1e-3; %from ms to s
   
    %------------------------------------------------
    % toc
    toc = msg(pos:pos+63); pos = pos + 64;
    toc = fliplr(reshape(toc,8,[]));                  % byte order inversion (little endian)
    toc = toc(:)';
    
    % floating point value decoding (double floating point)
    sign = str2num(toc(1));
    esp  = fbin2dec(toc(2:12));
    mant = fbin2dec(toc(13:64)) / 2^52;
    toc = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    toc = toc*1e-3; %from ms to s
    
    %------------------------------------------------
    % af2
    af2 = msg(pos:pos+31); pos = pos + 32;
    af2 = fliplr(reshape(af2,8,[]));                  % byte order inversion (little endian)
    af2 = af2(:)';

    % floating point value decoding (single floating point)
    sign = str2num(af2(1));
    esp  = fbin2dec(af2(2:9));
    mant = fbin2dec(af2(10:32)) / 2^23;
    af2 = (-1)^sign * (2^(esp - 127)) * (1 + mant);
    af2 = af2*1e3; %from ms/ms^2 to s/s^2
    
    %------------------------------------------------
    % af1
    af1 = msg(pos:pos+31); pos = pos + 32;
    af1 = fliplr(reshape(af1,8,[]));                  % byte order inversion (little endian)
    af1 = af1(:)';

    % floating point value decoding (single floating point)
    sign = str2num(af1(1));
    esp  = fbin2dec(af1(2:9));
    mant = fbin2dec(af1(10:32)) / 2^23;
    af1 = (-1)^sign * (2^(esp - 127)) * (1 + mant);
    
    %------------------------------------------------
    % af0
    af0 = msg(pos:pos+31); pos = pos + 32;
    af0 = fliplr(reshape(af0,8,[]));                  % byte order inversion (little endian)
    af0 = af0(:)';

    % floating point value decoding (single floating point)
    sign = str2num(af0(1));
    esp  = fbin2dec(af0(2:9));
    mant = fbin2dec(af0(10:32)) / 2^23;
    af0 = (-1)^sign * (2^(esp - 127)) * (1 + mant);
    af0 = af0*1e-3; %from ms to s
    
    %------------------------------------------------
    % URA (svaccur) 
    svaccur_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    svaccur_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    svaccur = svaccur_1 + (svaccur_2 * 2^8);                 % little endian

    %------------------------------------------------
    %IODE   
    IODE_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    IODE_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    IODE = IODE_1 + (IODE_2 * 2^8);                          % little endian
    
    %------------------------------------------------
    % IODC  
    IODC_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    IODC_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    IODC = IODC_1 + (IODC_2 * 2^8);                          % little endian

    %------------------------------------------------ 
    % code_on_L2
    code_on_L2_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    code_on_L2_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    code_on_L2 = code_on_L2_1 + (code_on_L2_2 * 2^8);        % little endian
    
    %------------------------------------------------ 
    % L2flag
    L2flag_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    L2flag_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    L2flag = L2flag_1 + (L2flag_2 * 2^8);                    % little endian
    
    %------------------------------------------------ 
    % weekno
    weekno_1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    weekno_2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    weekno = weekno_1 + (weekno_2 * 2^8);                    % little endian
 
    %output and reorder ephemerides data (if IODC == IODE)
    if (IODC == IODE)
        data{2}(1) = PRN;
        data{2}(2) = af2;
        data{2}(3) = M0;
        data{2}(4) = root_A;
        data{2}(5) = delta_n;
        data{2}(6) = e;
        data{2}(7) = omega;
        data{2}(8) = Cuc;
        data{2}(9) = Cus;
        data{2}(10) = Crc;
        data{2}(11) = Crs;
        data{2}(12) = i0;
        data{2}(13) = IDOT;
        data{2}(14) = Cic;
        data{2}(15) = Cis;
        data{2}(16) = omega0;
        data{2}(17) = omegadot;
        data{2}(18) = toe;
        data{2}(19) = af0;
        data{2}(20) = af1;
        data{2}(21) = toc;
        data{2}(22) = IODE;
        data{2}(23) = code_on_L2;
        data{2}(24) = weekno;
        data{2}(25) = L2flag;
        data{2}(26) = svaccur;
        data{2}(28) = tgd;
        data{2}(30) = PRN;       %assume only GPS (not multi-constellation)
        data{2}(31) = int8('G'); %assume only GPS (not multi-constellation)
        data{2}(32) = weektow2time(weekno, toe, 'G'); %wrong: weekno is mod(1024)... taken care of by the caller
        data{2}(33) = weektow2time(weekno, toc, 'G'); %wrong: weekno is mod(1024)... taken care of by the caller
    end

elseif (type == 2 && length(msg(pos:end)) >= 728)
    
    %disp('GLONAS PRN');
    %disp(PRN); 
    
end
