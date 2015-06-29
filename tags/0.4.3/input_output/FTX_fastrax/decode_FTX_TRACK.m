function [data] = decode_FTX_TRACK(msg, constellations)

% SYNTAX:
%   [data] = decode_FTX_TRACK(msg, constellations)
%
% INPUT:
%   msg = message transmitted by the Fastrax_IT03 receiver (string)
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the TRACK packet information
%          1.1) message class-id (TRACK)
%          2.1) Tick    = Common info for carrier and code.
%          2.2) NumChan = <15bit-reserved-8bit><7bit-Number of channel in Lock(track)-0bit>
%          2.3) Flags 	= Flags which are not channel-specific.
%          3.1) PRN     = space vehicle number
%          3.2) ChIdx   = channel index.
%          3.3) Lock 	= Channel lock status. See the lock bits definitions.
%          3.4) Power   = Power measurement value from Carrier tracking.	
%          3.5) CarrCount = Carrier NCO cycle count
%          3.6) CarrPhase = Carrier phase. The 8 LSB represent the fractional Carrier Count so that the CarrierPhase measurement can be obtained as: dwCarrCount + (wCarrPhase & 0x00ff)/256.
%          3.7) CarrFreq  = Carrier NCO replica frequency from loop filter.
%          3.8) PrnTow    = Last decoded TOW.
%          3.9) ChipPhase = Code chip phase, HW offset 1A. 
%          3.10) ChipCount   = Code chip count, HW offset 1B.
%          3.11) EpochCount  = Code epoch count, HW offset 1C.
%          3.12) CyclesSLTOW = Code cycles since last decoded TOW.
%          3.13) Reserved1   = Signal-to-noise ratio (dBHz).
%          3.14) CodeNCOFreq = Code NCO replica frequency from loop filter.
%
% DESCRIPTION:
%   TRACK binary message decoding.
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

debug = 0;

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

% first message initial index
pos = 1;

% output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(3,1);
data{3} = zeros(constellations.nEnabledSat,14);

% output header
data{1} = 'TRACK';

% TRACK_HDR TRACK.TrackHdr: Common info for carrier and code.		
% TRACK.TrackHdr.dwTick 	DWORD 	Measurement time (system ticker value)
[Tick, pos]     = FTX_TypeConv('DWORD', msg, pos);

% TRACK.TrackHdr.swNumChan 	WORD 	<15bit-reserved-8bit><7bit-Number of channel in Lock(track)-0bit>
[NumChan, pos, NumChan_par] = FTX_TypeConv('WORD', msg, pos);
NumChan = NumChan_par(1);           

% TRACK.TrackHdr.swFlags 	WORD 	Flags which are not channel-specific.
[Flags, pos]	= FTX_TypeConv('WORD', msg, pos);

%output data save
data{2}(1) = Tick;
data{2}(2) = NumChan;
data{2}(3) = Flags;

for j = 1 : NumChan
    % TRACK_CH TRACK.TrackCh[12]: Channel-specific data.
    % TRACK.TrackCh[n].swLockedChInfo 	WORD 	Channel info. Low byte: PRN number. High byte: channel index.
    [Temp, pos, LChInfo]    = FTX_TypeConv('WORD', msg, pos);
    PRN     = LChInfo(1);
    ChIdx   = LChInfo(2);
    
    % TRACK.TrackCh[n].swLock 	WORD 	Channel lock status. See the lock bits definitions.
    [Lock, pos]      = FTX_TypeConv('WORD', msg, pos);
    if (debug == 1)
        fprintf('Lock: \t%6d\t  ', Lock);
    end
    
    % TRACK_CARRIER TRACK.TrackCh[n].TrackCarrier: Carrier tracking data		
    % TRACK.TrackCh[n].TrackCarrier.dwDet 	DWORD 	Power measurement value from Carrier tracking.
    [Power, pos]     = FTX_TypeConv('DWORD', msg, pos);

    % TRACK.TrackCh[n].TrackCarrier.dwCarrCount 	DWORD 	Carrier NCO cycle count
    [CarrCount, pos] = FTX_TypeConv('DWORD', msg, pos);
    if (debug == 1)
        fprintf('CarrCount: \t%9d\t  ', CarrCount);
    end
    
    % TRACK.TrackCh[n].TrackCarrier.wCarrPhase 	WORD 	Carrier phase. The 8 LSB represent the fractional Carrier Count so that the CarrierPhase measurement can be obtained as: dwCarrCount + (wCarrPhase & 0x00ff)/256.
    [CarrPhase, pos, CarrPhase_par] = FTX_TypeConv('WORD', msg, pos);
    CarrPhase = CarrCount + (CarrPhase_par(1)/256);
    if (debug == 1)
        fprintf('CarrPhase: \t%15.4f\t ', CarrPhase);
    end
    
    % TRACK.TrackCh[n].TrackCarrier.dwCarrNCOFreq 	DWORD 	Carrier NCO replica frequency from loop filter.
    [CarrFreq, pos]  = FTX_TypeConv('DWORD', msg, pos);
    
    % TRACK_CODE TRACK.TrackCh[n].TrackCode: Code tracking data		
    % TRACK.TrackCh[n].TrackCode.dwPrnTow 	DWORD 	Last decoded TOW
    [PrnTow, pos]    = FTX_TypeConv('DWORD', msg, pos);
    if (debug == 1)
        fprintf('PrnTow: \t%11d\t ', PrnTow);
    end
    
    % CORR_CH_PRN_CNT TRACK.TrackCh[n].TrackCode.TrPrnCnt:		
    % TRACK.TrackCh[n].TrackCode.TrPrnCnt.wChipPhase 	WORD 	Code chip phase, HW offset 1A
    [ChipPhase, pos]      = FTX_TypeConv('WORD', msg, pos);     
    % TRACK.TrackCh[n].TrackCode.TrPrnCnt.wChipCount 	WORD 	Code chip count, HW offset 1B
    [ChipCount, pos]      = FTX_TypeConv('WORD', msg, pos);     
    % TRACK.TrackCh[n].TrackCode.TrPrnCnt.wEpochCount 	WORD 	Code epoch count, HW offset 1C
    [EpochCount, pos]     = FTX_TypeConv('WORD', msg, pos);     
    % TRACK.TrackCh[n].TrackCode.wCyclesSinceLastTOW 	WORD 	Code cycles since last decoded TOW
    [CyclesSLTOW, pos]    = FTX_TypeConv('WORD', msg, pos);     
    % TRACK.TrackCh[n].TrackCode.wReserved1 	WORD 	Signal-to-noise ratio (dBHz)
    [Reserved1, pos]      = FTX_TypeConv('WORD', msg, pos);
    % TRACK.TrackCh[n].TrackCode.dwCodeNCOFreq 	DWORD 	Code NCO replica frequency from loop filter
    [CodeNCOFreq, pos]    = FTX_TypeConv('DWORD', msg, pos);
    
    % assign constellation-specific indexes
    idx = [];
    if (SV <= 32)
        idx = constellations.GPS.indexes(SV);
    end
    
    data{3}(idx, 1) = PRN;
    data{3}(idx, 2) = ChIdx;
    data{3}(idx, 3) = Lock;
    data{3}(idx, 4) = Power;
    data{3}(idx, 5) = CarrCount;
    data{3}(idx, 6) = CarrPhase;
    data{3}(idx, 7) = CarrFreq;
    data{3}(idx, 8) = PrnTow;
    data{3}(idx, 9) = ChipPhase;
    data{3}(idx,10) = ChipCount;
    data{3}(idx,11) = EpochCount;
    data{3}(idx,12) = CyclesSLTOW;
    data{3}(idx,13) = Reserved1;
    data{3}(idx,14) = CodeNCOFreq;
end