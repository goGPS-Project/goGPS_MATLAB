function [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3_time, SP3_clck, err_tropo, err_iono, dtR)
% SYNTAX:
%   [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3_time, SP3_clck, err_tropo, err_iono, dtR);
%
% INPUT:
%
% OUTPUT:
%
% DESCRIPTION:
%   Compute initial receiver and satellite position and clocks using
%   Bancroft and least-squares iterative correction. It requires at least
%   four satellites available.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) Mirko Reguzzoni, Eugenio Realini, 2012
%
%----------------------------------------------------------------------------------------------

global v_light

time_tx_RAW = time_rx - (range - err_tropo - err_iono) / v_light + dtR;

%     tcorr0 = 0;
%     tcorr = sat_clock_error_correction(tx_RAW, Eph);
%     while (tcorr - tcorr0 > 0.001)
%         tcorr0 = tcorr;
%         tcorr = sat_clock_error_correction(tx_RAW - tcorr0, Eph);
%     end

if (isempty(SP3_time))
    dtS = sat_clock_error_correction(time_tx_RAW, Eph);
    dtS = sat_clock_error_correction(time_tx_RAW - dtS, Eph);
else
    %interpolate SP3 satellite clock correction term
    dtS = interpolate_SP3_clock(time_tx_RAW, SP3_time, SP3_clck(sat,:));
end

time_tx = time_tx_RAW - dtS;
