function [data] = decode_FTX_PSEUDO(msg, constellations)

% SYNTAX:
%   [data] = decode_FTX_PSEUDO(msg, constellations)
%
% INPUT:
%   msg = message transmitted by the Fastrax_IT03 receiver (string)
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the PSEUDO packet information
%          1.1) message class-id (PSEUDO)
%          2.1) TowMs   = Time of week [s] (TowMs+TowFrac)
%          2.2) TowWeek = GPS week - The GPS Week Number count began at approximately midnight on the evening of 05 January 1980 / morning of 06 January 1980.
%          2.3) NumObs  = Number of observations in this message.
%          2.4) Tick         = System tick at the measurement time. (Not Used).
%          2.5) Flags        = Flags common to all observations. (Not Used).
%          2.6) MMP          = Common millisecond multiple of all pseudoranges. (Not Used).
%          2.7) MeasLength   = Duration in ticks of the measurement interval for integrated Doppler. (Not Used).
%          3.1) CarrierPhase = Carrier phase (in degree). ------------------------------------------------------- (in cycles)
%          3.2) PseudoRange  = Code pseudorange [m] (corrected and possibly smoothed). -------------------------- (C/A code in meters)
%          3.3) Doppler      = Doppler [m/s]. -------------------------------------------------------------------- (in Hertz)
%          3.4) PRN          = space vehicle number
%          3.5) ObsFlags     = Status bitfield. The values of this field are defined by the PSEUDO_OBS_ flags.
%          3.6) SNR          = S/N [dBHz] range [0,63]
%          3.7) Corrections        = Which corrections have been applied. See PSEUDO_OBS_CORRECTED_ flags. 
%          3.8) LoopDopplerOffset  = Difference between Doppler measurement and the frequency from the software carrier tracking loop [cm/s].   
%          3.9) RangeErrEstim      = Error estimate for range measurement.
%          3.10) RateErrEstim      = Error estimate rate measurement.
%          3.11) EpochCount        = Counter of PRN code epochs since channel initialisation. The least significant 15 bits contain the epoch count (which will wrap around after 32767 epochs). The most significant bit is set to 1 after the first wraparound has happened.
%
% DESCRIPTION:
%   PSEUDO binary message decoding.
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

%retrieve GPS L1 wavelength
lambda = goGNSS.getWavelength(goGNSS.ID_GPS, 1);

% first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(8,1);
data{3} = zeros(constellations.nEnabledSat,11);

%output data save
data{1} = 'PSEUDO';

% PSEUDO.dwTick 	DWORD 	System tick at the measurement time.
[Tick, pos]     = FTX_TypeConv('DWORD', msg, pos);

% PSEUDO.wNumObs 	WORD 	Number of observations in this message.
[NumObs, pos]   = FTX_TypeConv('WORD', msg, pos);

% PSEUDO.swFlags 	WORD 	Flags common to all observations. See PSEUDO_ flags.
[Flags, pos]    = FTX_TypeConv('WORD', msg, pos);

% -- INT_GPS_TIME PSEUDO.RefGpsTow: The reference TOW with respect to which the pseudoranges are computed. This time is not necessarily a good indication of the receiver or satellite TOW.
% -- GPS_TIME PSEUDO.RefGpsTow.t: GPS week & TOW in ms.
% PSEUDO.RefGpsTow.t.wWeek      WORD 	Full week number.
[TowWeek, pos]  = FTX_TypeConv('WORD', msg, pos);

% PSEUDO.RefGpsTow.t.dwTowMs 	DWORD  	Time of week in milliseconds.
[TowMs, pos]    = FTX_TypeConv('DWORD', msg, pos);
TowMs = TowMs/1000;

% PSEUDO.RefGpsTow.lTowFrac 	INT32  	Sub-millisecond fractional part [ps].
[TowFrac, pos]  = FTX_TypeConv('INT32', msg, pos);
TowFrac = TowFrac/100000000000;
TowMs = TowMs + TowFrac;

% PSEUDO.lPseudoRangeMultiple 	INT32  	Common millisecond multiple of all pseudoranges.
[MMP, pos]      = FTX_TypeConv('INT32', msg, pos);
MMP = MMP/1000;
%TowMs = TowMs + MMP;

% PSEUDO.wMeasLength 	WORD 	Duration in ticks of the measurement interval for integrated Doppler.
[MeasLength, pos] = FTX_TypeConv('WORD', msg, pos);

% NOT IMPLEMENTED - PSEUDO.wReserved4 	WORD 	Reserved for future use.
pos = pos + 16;

%output data save
data{2}(1) = TowMs;
data{2}(2) = TowWeek;
data{2}(3) = NumObs;
data{2}(4) = Tick;
data{2}(5) = Flags;
data{2}(6) = MMP;
data{2}(7) = MeasLength;

% PSEUDO_OBS PSEUDO.Obs[12]: Observations.
for j = 1 : NumObs
    % PSEUDO.Obs[n].wPrn 	WORD 	LS byte: SV PRN. MS byte: reserved.
    [PRN, pos,PRN_part] = FTX_TypeConv('WORD', msg, pos);
    PRN = PRN_part(1);
    clear PRN_part    
    
    % PSEUDO.Obs[n].wSnr 	WORD 	LS byte: S/N [dBHz] range [0,63]. MS byte: reserved.
    [SNR, pos,SNR_part] = FTX_TypeConv('WORD', msg, pos);
    SNR = SNR_part(1);
    clear SNR_part

    % PSEUDO.Obs[n].swFlags 	WORD 	Status bitfield. The values of this field are defined by the PSEUDO_OBS_ flags.
    [ObsFlags, pos]     = FTX_TypeConv('WORD', msg, pos);

    % PSEUDO.Obs[n].swCorrections 	WORD 	Which corrections have been applied. See PSEUDO_OBS_CORRECTED_ flags.
    [Corrections, pos]  = FTX_TypeConv('WORD', msg, pos);

    % PSEUDO.Obs[n].dPseudoRange 	DOUBLE 	Code pseudorange [m] (corrected and possibly smoothed).
    [PseudoRange, pos]  = FTX_TypeConv('DOUBLE', msg, pos);
    
    % PSEUDO.Obs[n].dDoppler 	DOUBLE 	Doppler [m/s].
    [Doppler, pos]      = FTX_TypeConv('DOUBLE', msg, pos);
    Doppler = Doppler / lambda;

    % PSEUDO.Obs[n].iLoopDopplerOffset 	INT16 	Difference between Doppler measurement and the frequency from the software carrier tracking loop [cm/s].   
    [LoopDopplerOffset, pos] = FTX_TypeConv('INT16', msg, pos);
    
    % NOT IMPLEMENTED - PSEUDO.Obs[n].iOffsetEstim 	INT16 	Reserved for future use.
    pos = pos +16;
    
    % PSEUDO.Obs[n].iCarrierPhase 	INT16 	Carrier phase. (degree)
    [CarrierPhase, pos]  = FTX_TypeConv('INT16', msg, pos);
    %CarrierPhase = CarrierPhase / 360;                      %Convert to cycles

    % PSEUDO.Obs[n].wRangeErrEstim 	WORD 	Error estimate for range measurement.
    [RangeErrEstim, pos] = FTX_TypeConv('WORD', msg, pos);

    % PSEUDO.Obs[n].wRateErrEstim 	WORD 	Error estimate rate measurement.
    [RateErrEstim, pos]  = FTX_TypeConv('WORD', msg, pos);

    % PSEUDO.Obs[n].wEpochCount 	WORD 	Counter of PRN code epochs since channel initialisation. The least significant 15 bits contain the epoch count (which will wrap around after 32767 epochs). The most significant bit is set to 1 after the first wraparound has happened.
    wraparound  = msg(pos);                  pos = pos + 1;     % wraparound
    EpochCount1 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;
    EpochCount2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    EpochCount  = EpochCount1 + (EpochCount1 * 2^8);
    clear EpochCount1 EpochCount2
        
    % NOT IMPLEMENTED - PSEUDO.Obs[n].dwReserved4 	DWORD 	Reserved for future use.
    pos = pos +32;

    % FLAG CONTROL - ObsFlags - view end of file
    ObsFlagsBin = dec2bin(ObsFlags);
    [temp,zeros_idx] = find(ObsFlagsBin == '0');
    zeros_idx = (length(ObsFlagsBin) - zeros_idx) + 1;
    if min(zeros_idx) <= 4
        % Data Error --> Impose SNR = 0
        SNR = 0;
    end
    if min(zeros_idx) == 6
        % Data Error --> Impose SNR = 0
        SNR = 0;
    end

    % PSEUDO_OBS_DOPPLER_OK 	0x0001 	The doppler is valid
    % PSEUDO_OBS_PSEUDORANGE_OK 	0x0002 	The pseudorange is valid
    % PSEUDO_OBS_TOW_OK 	0x0004 	The pseudorange has been computed from a signal where a valid TOW has been decoded.
    % PSEUDO_OBS_PRN_OK 	0x0008 	The PRN number is one of the supported PRN numbers.
    % PSEUDO_OBS_ELEV_OK 	0x0010 	The SV is higher than the elevation mask
    % PSEUDO_OBS_SNR_OK 	0x0020 	The signal is not too weak.
    % PSEUDO_OBS_SV_HEALTHY 	0x0040 	The SV is marked as healthy
    % PSEUDO_OBS_NO_CROSS_CORR 	0x0080 	No cross correlation has been detected.
    % PSEUDO_OBS_DATA_EXISTS 	0x0100 	There is valid data for this observation.
    % PSEUDO_OBS_DATA_GOOD 	0x0200 	The data for this observation is accurate enough for navigation.
    % PSEUDO_OBS_BIT_LOCK 	0x0400 	Bit edges are being detected on the channel.
    % PSEUDO_OBS_FIRST_MEAS 	0x0800 	This is the first measurement for this PRN.
    % PSEUDO_OBS_RAIM_P_OK 	0x1000 	Observation not considered as outlier for position by RAIM.
    % PSEUDO_OBS_RAIM_V_OK 	0x2000 	Observation not considered as outlier for velocity by RAIM.
    % PSEUDO_OBS_RAIM_T_OK 	0x4000 	Observation not considered as outlier for time by RAIM.
    % PSEUDO_OBS_PLL 	0x8000U 	The doppler measurement has been obtained with the phase locked loop (PLL).
    % PSEUDO_OBS_MEAS_OK 	( PSEUDO_OBS_ELEV_OK | PSEUDO_OBS_SNR_OK | PSEUDO_OBS_PRN_OK | PSEUDO_OBS_NO_CROSS_CORR | PSEUDO_OBS_SV_HEALTHY | PSEUDO_OBS_DATA_EXISTS | PSEUDO_OBS_DATA_GOOD | PSEUDO_OBS_PSEUDORANGE_OK ) 	A common mask that can be used to determine if the pseudorange measurement is valid for navigation.
    % PSEUDO_OBS_DOPPLER_MEAS_OK 	( PSEUDO_OBS_ELEV_OK | PSEUDO_OBS_SNR_OK | PSEUDO_OBS_PRN_OK | PSEUDO_OBS_NO_CROSS_CORR | PSEUDO_OBS_SV_HEALTHY | PSEUDO_OBS_DATA_EXISTS | PSEUDO_OBS_DATA_GOOD | PSEUDO_OBS_DOPPLER_OK ) 	A common mask that can be used to determine if the pseudorange measurement is valid for navigation.

    % assign constellation-specific indexes
    % exclude satelites with wrong data
    idx = [];
    if (SV <= 32 && (SNR > 0))
        idx = constellations.GPS.indexes(SV);
    end
    
    %data output save
    data{3}(idx, 1) = CarrierPhase;
    data{3}(idx, 2) = PseudoRange;
    data{3}(idx, 3) = Doppler;
    data{3}(idx, 4) = PRN;
    data{3}(idx, 5) = ObsFlags;
    data{3}(idx, 6) = SNR;
    data{3}(idx, 7) = Corrections;
    data{3}(idx, 8) = LoopDopplerOffset;
    data{3}(idx, 9) = RangeErrEstim;
    data{3}(idx,10) = RateErrEstim;
    data{3}(idx,11) = EpochCount;
end

end


% 
% Data flags in the PSEUDO structure
% 
% These flags indicate the quality of the data in the PSEUDO .
% Name 	Value 	Description
% PSEUDO_OBS_DOPPLER_OK 	0x0001 	The doppler is valid
% PSEUDO_OBS_PSEUDORANGE_OK 	0x0002 	The pseudorange is valid
% PSEUDO_OBS_TOW_OK 	0x0004 	The pseudorange has been computed from a signal where a valid TOW has been decoded.
% PSEUDO_OBS_PRN_OK 	0x0008 	The PRN number is one of the supported PRN numbers.
% PSEUDO_OBS_ELEV_OK 	0x0010 	The SV is higher than the elevation mask
% PSEUDO_OBS_SNR_OK 	0x0020 	The signal is not too weak.
% PSEUDO_OBS_SV_HEALTHY 	0x0040 	The SV is marked as healthy
% PSEUDO_OBS_NO_CROSS_CORR 	0x0080 	No cross correlation has been detected.
% PSEUDO_OBS_DATA_EXISTS 	0x0100 	There is valid data for this observation.
% PSEUDO_OBS_DATA_GOOD 	0x0200 	The data for this observation is accurate enough for navigation.
% PSEUDO_OBS_BIT_LOCK 	0x0400 	Bit edges are being detected on the channel.
% PSEUDO_OBS_FIRST_MEAS 	0x0800 	This is the first measurement for this PRN.
% PSEUDO_OBS_RAIM_P_OK 	0x1000 	Observation not considered as outlier for position by RAIM.
% PSEUDO_OBS_RAIM_V_OK 	0x2000 	Observation not considered as outlier for velocity by RAIM.
% PSEUDO_OBS_RAIM_T_OK 	0x4000 	Observation not considered as outlier for time by RAIM.
% PSEUDO_OBS_PLL 	0x8000U 	The doppler measurement has been obtained with the phase locked loop (PLL).
% PSEUDO_OBS_MEAS_OK 	( PSEUDO_OBS_ELEV_OK | PSEUDO_OBS_SNR_OK | PSEUDO_OBS_PRN_OK | PSEUDO_OBS_NO_CROSS_CORR | PSEUDO_OBS_SV_HEALTHY | PSEUDO_OBS_DATA_EXISTS | PSEUDO_OBS_DATA_GOOD | PSEUDO_OBS_PSEUDORANGE_OK ) 	A common mask that can be used to determine if the pseudorange measurement is valid for navigation.
% PSEUDO_OBS_DOPPLER_MEAS_OK 	( PSEUDO_OBS_ELEV_OK | PSEUDO_OBS_SNR_OK | PSEUDO_OBS_PRN_OK | PSEUDO_OBS_NO_CROSS_CORR | PSEUDO_OBS_SV_HEALTHY | PSEUDO_OBS_DATA_EXISTS | PSEUDO_OBS_DATA_GOOD | PSEUDO_OBS_DOPPLER_OK ) 	A common mask that can be used to determine if the pseudorange measurement is valid for navigation.
% Flags in the pseudodata header
% 
% These flags are common to all channels in the PSEUDO structure.
% Name 	Value 	Description
% PSEUDO_TOW_WEEK_OK 	0x0001 	This flag is set in the PSEUDO header if the reference TOW week is valid.
% PSEUDO_TOW_OK 	0x0002 	This flag is set in the PSEUDO header if the reference TOW is valid.
% PSEUDO_RESYNCH 	0x0004 	This flag is set in the PSEUDO header to indicate that the receiver tow estimate may be inconsistent with the previous one. This may happen, for exameple, after channels have been lost and reacquired.
% PSEUDO_FIRST_MEAS 	0x0008 	This flag is set in the PSEUDO header if this is the first PSEUDO message after the obs task has been started.
% PSEUDO_UNSCHEDULED 	0x0010 	This flag is set in the PSEUDO header if the corresponding NAV_FIX should not be output.
% Correction flags in the PSEUDO structure
% 
% These flags indicate the corrections that have been applied to the data in the PSEUDO structure.
% Name 	Value 	Description
% PSEUDO_OBS_CORRECTED_AMBIGUOUS 	0x0001 	The observation has been corrected for millisecond ambiguity.
% PSEUDO_OBS_CORRECTED_BY_SMOOTHING 	0x0002 	The observation has been smoothed.
% PSEUDO_OBS_CORRECTED_BY_IONO 	0x0008 	Ionosphere correction has been applied.
% PSEUDO_OBS_CORRECTED_BY_TROPO 	0x0010 	Troposphere correction has been applied.
% PSEUDO_OBS_CORRECTED_BY_FAST_CORR 	0x0020 	WAAS fast corrections and integrity info has been applied. If this bit is set when the PSEUDO_OBS_SV_HEALTHY is not set, it means that the satellite has been declared unhealthy by WAAS.
% PSEUDO_OBS_CORRECTED_BY_DGPS 	0x0040 	Differential GPS correction has been applied.
% PSEUDO_OBS_CORRECTED_BY_SLOW_CORR 	0x0080 	Waas long term corrections have been applied
% PSEUDO_OBS_CORRECTED_BY_WAAS_IONO 	0x0100 	Waas ionospheric corrections have been applied
% PSEUDO_OBS_CORR_FRAME_LOCK 	0x8000u 	The channel has frame lock.
% PSEUDO_OBS_CORR_POSSIBLE_XCORR 	0x4000u 	
% PSEUDO_OBS_CORRECTED_BY_WAAS 	( PSEUDO_OBS_CORRECTED_BY_WAAS_IONO | PSEUDO_OBS_CORRECTED_BY_FAST_CORR) 	TBA: for now only waas ionos corretions
