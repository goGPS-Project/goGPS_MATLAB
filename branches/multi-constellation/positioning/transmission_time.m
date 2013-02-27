function [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR)

% SYNTAX:
%   [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR);
%
% INPUT:
%   time_rx   = reception time
%   range     = observed range
%   sat       = satellite index
%   Eph       = ephemeris
%   SP3       = structure containing precise ephemeris data
%   sbas      = SBAS corrections
%   err_tropo = tropospheric delays
%   err_iono  = ionospheric delays
%   dtR       = receiver clock offset
%
% OUTPUT:
%   time_tx = transmission time
%   dtS     = satellite clock error
%
% DESCRIPTION:
%   Compute the signal transmission time.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) Mirko Reguzzoni, Eugenio Realini, 2012
%
%----------------------------------------------------------------------------------------------

global v_light

%SBAS clock offsets
if (~isempty(sbas))
    dtsbas = sbas.doffset(sat);
else
    dtsbas = zeros(1,length(sat));
end

time_tx_RAW = time_rx - (range - err_tropo - err_iono) / v_light + dtR;

% tcorr0 = 0;
% tcorr = sat_clock_error_correction(time_tx_RAW, Eph);
% while (tcorr - tcorr0 > 0.001)
%     tcorr0 = tcorr;
%     tcorr = sat_clock_error_correction(time_tx_RAW - tcorr0, Eph);
% end

if (isempty(SP3))
    
    %if GLONASS
    if (strcmp(char(Eph(31)),'R'))
        
        dtS = sat_clock_error_correction(time_tx_RAW, Eph);
        dtS = sat_clock_error_correction(time_tx_RAW - dtS, Eph);
    else
        %relativistic correction term
        dtrel = relativistic_error_correction(time_tx_RAW, Eph);
        
        %group delay correction term
        %--- this correction term is only for the benefit of "single-frequency"
        %    (L1 P(Y) or L2 P(Y)) users
        tgd = Eph(28);
        
        dtS = sat_clock_error_correction(time_tx_RAW, Eph);
        dtS = dtS + dtrel - tgd + dtsbas;
        dtS = sat_clock_error_correction(time_tx_RAW - dtS, Eph);
        dtS = dtS + dtrel - tgd + dtsbas;
    end
else
    %interpolate SP3 satellite clock correction term
    dtS = interpolate_SP3_clock(time_tx_RAW, SP3, sat);
end

time_tx = time_tx_RAW - dtS;
