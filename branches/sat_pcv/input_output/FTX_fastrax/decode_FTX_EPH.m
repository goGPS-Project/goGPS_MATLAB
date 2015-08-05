function [data] = decode_FTX_EPH(msg, constellations)

% SYNTAX:
%   [data] = decode_FTX_EPH(msg, constellations);
%
% INPUT:
%   msg = message transmitted by the fastrax receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the FTX_EPH packet information
%          1.1) message class-id (FTX_EPH)
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
%          2.27)GPS svhealth;
%          2.28)GPS tgd;
%          2.29)GPS fit_int;
%          2.30)multi-constellation satellite index (here only GPS is assumed)
%
% DESCRIPTION:
%   FTX_EPH binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Ivan Reguzzoni
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

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

% first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(33,1);

%output data save
data{1} = 'FTX-EPH';

% RAW_EPHEMERIS.wPrn 		WORD 	The PRN for which the data is intended.
[PRN, pos]          = FTX_TypeConv('WORD', msg, pos);

% RAW_EPHEMERIS.wFitPeriod 		WORD 	Orbit fit period, for downloaded data use X hours.
[FitPeriod, pos]	= FTX_TypeConv('WORD', msg, pos);
% Curve Fit Intervals. Bit 17 in word 10 of subframe 2 is a "fit interval" flag which indicates the curvefit
% interval used by the CS (Block II/IIA/IIR/IIR-M/IIF) and SS (Block IIIA) in determining the ephemeris
% parameters, as follows:
% 0 = 4 hours,
% 1 = greater than 4 hours.
if FitPeriod  == 4
    FitPeriod = 0;
else
    FitPeriod = 1;
end

% RAW_EPHEMERIS.wHealth 		WORD 	Satellite health status from subframes 4 and 5.
[Health, pos]       = FTX_TypeConv('WORD', msg, pos);

% RAW_EPHEMERIS.wFlags 		WORD 	Flags that indicate the contents of the structure (see EPH_ flags).
[Flags, pos]        = FTX_TypeConv('WORD', msg, pos);

% RAW_SV_CLOCK_SUB RAW_EPHEMERIS.ClockSub: The clock data which are common with the RAW_ALMANAC structure.		
% RAW_EPHEMERIS.ClockSub.lAf0 	INT32 	SV Clock Correction: (sec)
[Af0, pos]          = FTX_TypeConv('INT32', msg, pos);
Af0 = Af0 * (2^(-31));

% RAW_EPHEMERIS.ClockSub.iAf1 	INT16 	SV Clock Correction: (sec/sec)
[Af1, pos]          = FTX_TypeConv('INT16', msg, pos);
Af1 = Af1 * (2^(-43));

% GPS_TIME RAW_EPHEMERIS.ClockSub.Toc: Reference time for clock data (seconds). In case of almanacs, this field stores the TOA.		
% RAW_EPHEMERIS.ClockSub.Toc.wWeek 	WORD 	Full week number.
[TocWeek, pos]      = FTX_TypeConv('WORD', msg, pos);

% RAW_EPHEMERIS.ClockSub.Toc.dwTowMs 	DWORD 	Time of week in milliseconds.
[TowMs, pos]        = FTX_TypeConv('DWORD', msg, pos);
TowMs = TowMs/1000;

% RAW_SV_ORBIT_SUB RAW_EPHEMERIS.OrbitSub: The orbit data which are common with the RAW_ALMANAC structure.		
% RAW_EPHEMERIS.OrbitSub.lEcc 	DWORD 	Eccentricity
[Ecc, pos]          = FTX_TypeConv('DWORD', msg, pos);
Ecc = Ecc * (2^(-33));

% RAW_EPHEMERIS.OrbitSub.lI0 	INT32 	Inclination at reference time (rad). In case of almanacs, this contains the correction with respect to 0.3.
[I0, pos]           = FTX_TypeConv('INT32', msg, pos);
I0 = I0 * (2^(-31));
I0 = I0 * pi();

% RAW_EPHEMERIS.OrbitSub.lOmegaDot 	INT32 	Rate of right ascension (rad/s).
[OmegaDot, pos]     = FTX_TypeConv('INT32', msg, pos);  
OmegaDot = OmegaDot * (2^(-43));
OmegaDot = OmegaDot * pi();

% RAW_EPHEMERIS.OrbitSub.lSqrta 	DWORD 	Square root of the semimajor axis (m^1/2).
[Sqrta, pos]        = FTX_TypeConv('DWORD', msg, pos);
Sqrta = Sqrta * (2^(-19));

% RAW_EPHEMERIS.OrbitSub.lOmega0 	INT32 	Right ascension at reference time toe (rad).
[Omega0, pos]       = FTX_TypeConv('INT32', msg, pos);
Omega0 = Omega0 * (2^(-31));
Omega0 = Omega0 * pi();

% RAW_EPHEMERIS.OrbitSub.lOmega 	INT32 	Argument of perigee (rad).
[Omega, pos]        = FTX_TypeConv('INT32', msg, pos);
Omega = Omega * (2^(-31));
Omega = Omega * pi();

% RAW_EPHEMERIS.OrbitSub.lM0 	INT32 	Mean anomaly at reference time toe (rad).
[M0, pos]           = FTX_TypeConv('INT32', msg, pos);
M0 = M0 * (2^(-31));
M0 = M0 * pi();

% RAW_SV_CLOCK_EXT RAW_EPHEMERIS.ClockExt: Higher-order terms of the clock data.		
% RAW_EPHEMERIS.ClockExt.wIODC 	WORD 	Issue Of Data for Clock
[IODC, pos]         = FTX_TypeConv('WORD', msg, pos);

% RAW_EPHEMERIS.ClockExt.iGroupDelay 	INT16 	Group delay (seconds)
[GroupDelay, pos]	= FTX_TypeConv('INT16', msg, pos);
GroupDelay = GroupDelay * (2^(-31));

% RAW_EPHEMERIS.ClockExt.iAf2 	INT16 	SV Clock Correction: (sec/sec^2)
[Af2, pos]      = FTX_TypeConv('INT16', msg, pos);
Af2 = Af2 * (2^(-55));

% RAW_SV_ORBIT_EXT RAW_EPHEMERIS.OrbitExt: Higher-order terms of the orbit data.		
% RAW_EPHEMERIS.OrbitExt.dwToe 	DWORD 	Reference time for ephemeris data (seconds).
[Toe, pos]      = FTX_TypeConv('DWORD', msg, pos);

% RAW_EPHEMERIS.OrbitExt.wIODE 	WORD 	Issue Of Data for Ephemeris
[IODE, pos]     = FTX_TypeConv('WORD', msg, pos);

% RAW_EPHEMERIS.OrbitExt.iDeltan 	INT16 	Mean motion difference (rad/s)
[Deltan, pos]	= FTX_TypeConv('INT16', msg, pos);
Deltan = Deltan * (2^(-43));
Deltan = Deltan * pi();

% RAW_EPHEMERIS.OrbitExt.iCrs 	INT16 	sin harmonic correction to orbit radius (m)
[Crs, pos]      = FTX_TypeConv('INT16', msg, pos);
Crs = Crs * (2^(-5));

% RAW_EPHEMERIS.OrbitExt.iCrc 	INT16 	cos harmonic correction to orbit radius (m)
[Crc, pos]      = FTX_TypeConv('INT16', msg, pos);
Crc = Crc * (2^(-5));

% RAW_EPHEMERIS.OrbitExt.iCus 	INT16 	sin harmonic correction to argument of latitude (rad)
[Cus, pos]      = FTX_TypeConv('INT16', msg, pos);
Cus = Cus * (2^(-29));

% RAW_EPHEMERIS.OrbitExt.iCuc 	INT16 	cos harmonic correction to argument of latitude (rad)
[Cuc, pos]      = FTX_TypeConv('INT16', msg, pos);
Cuc = Cuc * (2^(-29));

% RAW_EPHEMERIS.OrbitExt.iCis 	INT16 	sin harmonic correction to inclination (rad)
[Cis, pos]      = FTX_TypeConv('INT16', msg, pos);
Cis = Cis * (2^(-29));

% RAW_EPHEMERIS.OrbitExt.iCic 	INT16 	cos harmonic correction to inclination (rad)
[Cic, pos]      = FTX_TypeConv('INT16', msg, pos);
Cic = Cic * (2^(-29));

% RAW_EPHEMERIS.OrbitExt.iIdot 	INT16 	Rate of inclination (rad/s)
[Idot, pos]     = FTX_TypeConv('INT16', msg, pos);
Idot = Idot * (2^(-43));
Idot = Idot * pi();

%output and reorder ephemerides data (if IODC == IODE)
if ((IODC == IODE) && (IODC == IODE) && constellations.GPS.enabled)
    data{2}(1) = PRN;
    data{2}(2) = Af2;
    data{2}(3) = M0;
    data{2}(4) = Sqrta;
    data{2}(5) = Deltan;
    data{2}(6) = Ecc;
    data{2}(7) = Omega;
    data{2}(8) = Cuc;
    data{2}(9) = Cus;
    data{2}(10) = Crc;
    data{2}(11) = Crs;
    data{2}(12) = I0;
    data{2}(13) = Idot;
    data{2}(14) = Cic;
    data{2}(15) = Cis;
    data{2}(16) = Omega0;
    data{2}(17) = OmegaDot;
    data{2}(18) = Toe;
    data{2}(19) = Af0;
    data{2}(20) = Af1;
    data{2}(21) = TowMs;
    data{2}(22) = IODE;                       
    data{2}(23) = 1;    % code_on_L2;   % Inconsistent with ublox - Default setting
    data{2}(24) = TocWeek;              % Inconsistent with ublox
    data{2}(25) = Flags;                % Inconsistent with ublox
    data{2}(26) = 0;                    % Inconsistent with ublox - Default setting
    data{2}(27) = Health;
    data{2}(28) = GroupDelay;
    data{2}(29) = FitPeriod;
    data{2}(30) = constellations.GPS.indexes(PRN);
    data{2}(31) = int8('G');
    data{2}(32) = weektow2time(TocWeek, Toe,   'G');
    data{2}(33) = weektow2time(TocWeek, TowMs, 'G');
end

% Check, no PRN --> delete header to improve performance
if sum(data{2,1}) == 0
    data{1} = '';
end
