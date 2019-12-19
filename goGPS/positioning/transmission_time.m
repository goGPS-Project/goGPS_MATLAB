function [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR, frequencies, obs_comb, lambda)

% SYNTAX:
%   [time_tx, dtS] = transmission_time(time_rx, range, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR, frequencies, obs_comb, lambda);
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
%   frequencies = L1 carrier (phase=1), L2 carrier (phase=2)
%   obs_comb    = observations combination (e.g. iono-free: obs_comb = 'IONO_FREE')
%   lambda      = matrix containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   time_tx = transmission time
%   dtS     = satellite clock error
%
% DESCRIPTION:
%   Compute the signal transmission time.

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

%SBAS clock offsets
if (~isempty(sbas))
    dtsbas = sbas.doffset(sat);
else
    dtsbas = zeros(1,length(sat));
end

time_tx_RAW = time_rx - (range - err_tropo - err_iono) / Core_Utils.V_LIGHT + dtR;

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
        dtrel = relativistic_clock_error_correction(time_tx_RAW, Eph, SP3);

        %group delay correction term
        if (nargin > 9 && ~strcmp(obs_comb,'IONO_FREE'))
            tgd = Eph(28);
            if (length(frequencies) == 2)
                %TO BE DONE!
            else
                if (frequencies(1) == 2)
                    gamma = (lambda(1,2)/lambda(1,1))^2;
                    tgd = gamma*tgd;
                end
            end
        else
            tgd = 0;
        end

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
