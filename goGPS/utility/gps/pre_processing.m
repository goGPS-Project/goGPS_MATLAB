function [time, pr1, ph1, pr2, ph2, XR, dtR, dtRdot, el, az, bad_sats, bad_epochs, var_dtR, var_SPP, status_obs, status_cs, eclipsed, ISBs, var_ISBs] = pre_processing(time_ref, time, XR0, pr1, ph1, pr2, ph2, dop1, dop2, snr1, Eph, SP3, iono, lambda, frequencies, obs_comb, nSatTot, waitbar_handle, flag_XR, sbas, cc, flag_full_prepro, order)

% SYNTAX:
%   [time, pr1, ph1, pr2, ph2, XR, dtR, dtRdot, el, az, bad_sats, bad_epochs, var_dtR, var_SPP, status_obs, status_cs, eclipsed, ISBs, var_ISBs] = pre_processing(time_ref, time, XR0, pr1, ph1, pr2, ph2, dop1, dop2, snr1, Eph, SP3, iono, lambda, frequencies, obs_comb, nSatTot, waitbar_handle, flag_XR, sbas, cc, flag_full_prepro, order);
%
% INPUT:
%   time_ref = GPS reference time
%   time = GPS nominal time (as read from RINEX file)
%   XR0 = receiver position (=[] if not available)
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio
%   Eph = matrix containing 33 ephemerides for each satellite
%   SP3 = structure with precise ephemeris and clock
%   iono = ionosphere parameters (Klobuchar)
%   frequencies = L1 carrier (phase=1) L2 carrier (phase=2)
%   obs_comb = observations combination (e.g. iono-free: obs_comb = 'IONO_FREE')
%   lambda  = wavelength matrix (depending on the enabled constellations)
%   nSatTot = maximum number of satellites (given the enabled constellations)
%   waitbar_handle = handle to the waitbar object
%   flag_XR = 2: coordinate fixed to XR0 values
%           = 1: approximate coordinates available
%           = 0: no apriori coordinates available
%   sbas = SBAS corrections
%   cc = struct with multi-constellation settings
%   flag_full_prepro = do  a full preprocessing
%   order = dynamic model order (1: static; >1 kinematic or epoch-by-epoch)

% OUTPUT:
%   time = GPS nominal time (as read from RINEX file) and then corrected to the reference time during the pre-processing
%   pr1 = processed code observation (L1 carrier)
%   ph1 = processed phase observation (L1 carrier)
%   pr2 = processed code observation (L2 carrier)
%   ph2 = processed phase observation (L2 carrier)
%   XR  = receiver position
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%   el  = elevation of the satellites
%   az  = azimuth of the satellites
%   bad_sats = vector for flagging "bad" satellites (e.g. too few observations, code without phase, etc)
%   bad_epochs  = vector with 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   var_SPP   = [code single point positioning a posteriori sigma, sum of
%                weighted squared residuals, redoundancy], one row per epoch
%   status_obs = for each satellite. NaN: not observed, 0: observed only, 1: used, -1: outlier
%   status_cs = [satellite_number, frequency, epoch, cs_correction_fix, cs_correction_float, cs_corrected(0:no, 1:yes)]: vector with cs information
%   eclipsed = satellites under eclipse condition (vector) (0: OK, 1: eclipsed)
%   ISBs        = estimated inter-system biases
%   var_ISBs    = variance of estimation errors (inter-system biases)
%
% DESCRIPTION:
%   Pre-processing of code and phase observations to correct them for
%    the receiver clock error.

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
%  Contributors:     Stefano Caldera
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

    global cutoff snr_threshold n_sys flag_doppler_cs

    state = Core.getCurrentSettings();
    log = Logger.getInstance();

    % Check doppler validity, if dop1 == dop2 [ !! ]
    if any(dop1(:)) && any(dop2(:)) && ~any(serialize(nan2zero(zero2nan(dop1)-zero2nan(dop2))))
        dop1 = dop1 * 0;
        dop2 = dop2 * 0;
    end
    p_rate = state.getProcessingRate();
    v_light = Core_Utils.V_LIGHT;

    %iono-free coefficients
    alpha1 = lambda(:,4);
    alpha2 = lambda(:,5);

    %number of epochs
    nEpochs = length(time);

    %receiver clock error
    dtR = zeros(nEpochs,1);
    XR = nan(3, nEpochs);

    %inter-system biases
    if (~isempty(SP3))
        nisbs = length(unique(SP3.sys(SP3.sys~=0)));
    else
        nisbs = length(unique(Eph(31,Eph(31,:)~=0)));
    end
    if (any(strfind(char(Eph(31,Eph(31,:)~=0)),'R')) || (~isempty(SP3) && any(strfind(char(SP3.sys(SP3.sys~=0))','R'))))
        if (cc.glo.isActive())
            if (nisbs == 1)
                nisbs = nisbs + 14; %if only GLONASS is present, just add IFBs
            else
                nisbs = nisbs - 1 + 14; %if GLONASS is present with other constellations, remove the GLONASS ISB and add IFBs
            end
        end
    end
    if (nisbs > 1)
        ISBs = zeros(nisbs-1, 1);
    else
        ISBs = [];
    end

    %receiver clock drift
    dtRdot = zeros(nEpochs-1,1);

    %vector for flagging "bad" satellites (e.g. too few observations, code without phase, etc)
    bad_sats = zeros(nSatTot,1);

    %vector with bad epochs definition
    bad_epochs=NaN(nEpochs,1);

    % vector with SPP a posteriori sigma
    var_SPP=NaN(nEpochs,3);

    %Lagrange interpolation order
    lagr_order = 10;

    %--------------------------------------------------------------------------------------------
    % APPROXIMATE POSITION
    %--------------------------------------------------------------------------------------------

    % if ((sum(abs(XR0)) == 0) || isempty(XR0))
    %     %approximate position not available
    %     flag_XR = 0;
    % else
    %     %approximate position available
    %     flag_XR = 1;
    % end

    err_iono = zeros(nSatTot,nEpochs);
    el = zeros(nSatTot,nEpochs);
    az = zeros(nSatTot,nEpochs);
    eclipsed = zeros(nSatTot,nEpochs);
    cond_num = zeros(nEpochs,1); %#ok<*NASGU>
    cov_XR = zeros(3,3,nEpochs);
    var_dtR = NaN(nEpochs,1);
    var_ISBs = NaN(nisbs-1,nEpochs);
    status_obs = NaN(nSatTot,nEpochs);
    status_cs=[];

    %----------------------------------------------------------------------------------------------
    % RECEIVER CLOCK DRIFT DISCONTINUITIES
    %----------------------------------------------------------------------------------------------
    % correct nominal time desynchronization and jumps ---------------------------------------------------
    
    % nominal time desynchronization (e.g. with some low-cost receivers)
    time_desync = round((time_ref - time) * 1e7) / 1e7; % computed here, used in init_positioning
    
    ph1_bk = ph1;
    ph2_bk = ph2;
    pr1_bk = pr1;
    pr2_bk = pr2;
    ph = zero2nan(bsxfun(@times, [ph1_bk; ph2_bk], [lambda(:, 1); lambda(:,2)])');
    pr = zero2nan([pr1_bk; pr2_bk]');
    %figure(1); plot(diff(zero2nan(ph),4));
    %figure(2); plot(diff(zero2nan(pr),4))
    [pr, ph, dt_pr, dt_ph] = Core_Pre_Processing.correctTimeDesync(time_ref, time, pr, ph);
    
    %flag_ph = false(size(ph));
    %flag_pr = false(size(pr));
    [ph, flag_ph] = Core_PP.flagRawObsD4(ph, time_ref - dt_ph, time_ref, 6, 5); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8
    [pr, flag_pr] = Core_PP.flagRawObsD4(pr, time_ref - dt_pr, time_ref, 6, 5); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8
    %[ph, flag_ph] = Core_PP.flagRawObsD4(ph, time_ref - dt_ph, time_ref, 6); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8
    %[pr, flag_pr] = Core_PP.flagRawObsD4(pr, time_ref - dt_pr, time_ref, 6); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8

    % flag by high deviation of the 4th derivate
    ph = nan2zero(bsxfun(@rdivide, zero2nan(ph), [lambda(:, 1); lambda(:,2)]'));
    ph1 = ph(:,1:size(ph1_bk,1))';
    ph2 = ph(:,(size(ph2_bk,1)+1):end)';
    
    % flag by high deviation of the 4th derivate
    pr1 = nan2zero(pr(:,1:size(pr1_bk,1))');
    pr2 = nan2zero(pr(:,(size(pr2_bk,1)+1):end)');
    
    %figure(3); plot(diff(zero2nan(ph),4));
    %figure(4); plot(diff(zero2nan(pr),4))

    % ----------------------------------------------------------------------------------------------------
        
    % [pr1, ph1] = correct_time_desync(time_ref, time, pr1, ph1, lambda(:,1));
    % [pr2, ph2] = correct_time_desync(time_ref, time, pr2, ph2, lambda(:,2));
    % time = time_ref;

    % remove short arcs
    min_arc = max([state.getMinArc() lagr_order]);
    % log.addMessage(sprintf('Trimming arcs shorter than %d epochs', min_arc));
    pr1 = remove_short_arcs(pr1, min_arc);
    pr2 = remove_short_arcs(pr2, min_arc);
    ph1 = remove_short_arcs(ph1, min_arc);
    ph2 = remove_short_arcs(ph2, min_arc);

    if not(flag_full_prepro)
        XR = repmat(XR0, 1, length(dtR));
        dtRdot(end+1) = dtRdot(end);
    else

        %------------------------------------------------------------------------------------------------------------------
        % RECEIVER CLOCK  AND INTER-SYSTEM BIAS (IF ANY) ESTIMATION BY MULTI-EPOCH LEAST-SQUARES ADJUSTMENT: INITIALIZATION
        %------------------------------------------------------------------------------------------------------------------

        %modulo time
        interval = median(diff(time)); %seconds
        mt = ceil(interval/3);
        nEpochs_reduced = floor(length(time)/mt);

        %number of position solutions to be estimated
        if (order == 1)
            npos = 1;
        else
            npos = nEpochs_reduced;
        end

        %vector to keep track of the number of observations available at each epoch
        n_obs_epoch = NaN(nEpochs_reduced, 1);

        %cumulative number of observations, to keep track of where each epoch begins
        epoch_index = NaN(nEpochs_reduced, 1);
        epoch_track = 0;

        %total number of observations (for matrix initialization)
        n_obs_tot = floor(sum(sum(pr1(:,:,1) ~= 0))/mt);

        if (nisbs > 1)
            y0_all = NaN(n_obs_tot,1);
            b_all  = NaN(n_obs_tot,1);
            A_all  = NaN(n_obs_tot,npos*3+nEpochs_reduced+(nisbs-1));
            Q_all  = NaN(n_obs_tot,1);
        end

        r = 1;
        for i = 1 : nEpochs

            %--------------------------------------------------------------------------------------------
            % SATELLITE AND EPHEMERIS SELECTION
            %--------------------------------------------------------------------------------------------

            if (frequencies(1) == 1)
                if (length(frequencies) < 2 || ~strcmp(obs_comb,'IONO_FREE'))
                    sat0 = find(pr1(:,i) ~= 0);
                else
                    sat0 = find(pr1(:,i) ~= 0 & pr2(:,i) ~= 0);
                end
            else
                sat0 = find(pr2(:,i) ~= 0);
            end

            status_obs(sat0,i) = 0; % satellite observed

            Eph_t = rt_find_eph (Eph, time(i), nSatTot);
            sbas_t = find_sbas(sbas, i);

            %----------------------------------------------------------------------------------------------
            % EPOCH-BY-EPOCH RECEIVER POSITION AND CLOCK ERROR
            %----------------------------------------------------------------------------------------------

            min_nsat_LS = 3 + n_sys;

            if (length(sat0) >= min_nsat_LS)
                if (frequencies(1) == 1)
                    if (length(frequencies) < 2 || ~strcmp(obs_comb,'IONO_FREE'))
                        [XR_tmp, dtR_tmp, ~, ~, ~, ~, ~, ~, err_iono_tmp, sat, el_tmp, az_tmp, ~, ~, cov_XR_tmp, var_dtR_tmp, ~, ~, ~, cond_num_tmp, bad_sat_i, bad_epochs(i), var_SPP(i,:), ~, eclipsed_tmp, ISBs_tmp, var_ISBs_tmp, y0, b, A, Q] = init_positioning(time_ref(i) - dt_pr(i) - time_desync(i), pr1(sat0,i), snr1(sat0,i), Eph_t, SP3, iono, sbas_t, XR0, [], [], sat0, [], lambda(sat0,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 0, nisbs > 1, 1); %#ok<ASGLU>
                    else
                        [XR_tmp, dtR_tmp, sat_pos_tmp, ~, ~, sat_vel_tmp, ~, ~, err_iono_tmp, sat, el_tmp, az_tmp, ~, ~, cov_XR_tmp, var_dtR_tmp, ~, ~, ~, cond_num_tmp, bad_sat_i, bad_epochs(i), var_SPP(i,:), ~, eclipsed_tmp, ISBs_tmp, var_ISBs_tmp, y0, b, A, Q] = init_positioning(time_ref(i) - dt_pr(i) - time_desync(i), alpha1(sat0).*pr1(sat0,i) - alpha2(sat0).*pr2(sat0,i), snr1(sat0,i), Eph_t, SP3, zeros(8,1), sbas_t, XR0, [], [], sat0, [], zeros(length(sat0),2), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 0, nisbs > 1, 1); %#ok<ASGLU>
                        %ids = find(sat == 2);
                        %if ~isempty(ids)
                        %    sat_pos(i,:) = sat_pos_tmp(ids, : );
                        %    sat_vel(i,:) = sat_vel_tmp(ids, : );
                        %end
                    end
                else
                    [XR_tmp, dtR_tmp, ~, ~, ~, ~, ~, ~, err_iono_tmp, sat, el_tmp, az_tmp, ~, ~, cov_XR_tmp, var_dtR_tmp, ~, ~, ~, cond_num_tmp, bad_sat_i, bad_epochs(i), var_SPP(i,:), ~, eclipsed_tmp, ISBs_tmp, var_ISBs_tmp, y0, b, A, Q] = init_positioning(time_ref(i) - dt_pr(i) - time_desync(i), pr2(sat0,i), snr1(sat0,i), Eph_t, SP3, iono, sbas_t, XR0, [], [], sat0, [], lambda(sat0,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 0, nisbs > 1, 1); %#ok<ASGLU>
                end

                if (~isempty(A) && (nisbs > 1) && (mod(i,mt) == 0))
                    n_obs_epoch(r) = length(y0);
                    y0_all(epoch_track+1:epoch_track+n_obs_epoch(r)) = y0;
                    b_all( epoch_track+1:epoch_track+n_obs_epoch(r)) =  b;
                    A_all( epoch_track+1:epoch_track+n_obs_epoch(r),(1:3)+any(npos-1)*(r-1)*3) = A(:,1:3);
                    if (nisbs > 1)
                        A_all( epoch_track+1:epoch_track+n_obs_epoch(r),end-(nisbs-2):end) = A(:,end-(nisbs-2):end);
                    end
                    Q_all( epoch_track+1:epoch_track+n_obs_epoch(r),1) = diag(Q);
                    epoch_index(r) = epoch_track + n_obs_epoch(r);
                    epoch_track = epoch_index(r);
                    r = r + 1;
                end

                if isempty(var_dtR_tmp)
                    var_dtR_tmp=NaN;
                end

                status_obs(sat,i) = 1; % satellite used
                status_obs(find(bad_sat_i==1),i)=-1; %#ok<FNDSB> % satellite outlier

                if (~isempty(dtR_tmp) && ~isempty(sat))
                    XR(:,i) = XR_tmp;
                    dtR(i) = dtR_tmp;
                    %             if (~isempty(ISBs_tmp))
                    %                 ISBs(:,i) = ISBs_tmp;
                    %             end
                    err_iono(sat,i) = err_iono_tmp;
                    el(sat,i) = el_tmp;
                    az(sat,i) = az_tmp;
                    eclipsed(sat,i) = eclipsed_tmp;
                    %             cond_num(i,1) = cond_num_tmp;
                    %             if (~isempty(var_dtR_tmp))
                    %                 cov_XR(:,:,i) = cov_XR_tmp;
                    %                 var_dtR(i,1) = var_dtR_tmp;
                    %             end
                end

                if (size(sat,1) >= min_nsat_LS)

                    if (i > 1)
                        %compute receiver clock drift
                        if (dtR(i) ~= 0 && dtR(i-1) ~= 0)
                            dtRdot(i-1) = (dtR(i) - dtR(i-1))/(time(i) - time(i-1));
                        end
                    end
                else
                    if (i > 2)
                        dtR(i) = 0;
                        dtRdot(i-1) = 0;
                        time(i) = 0;
                    end
                end
            else
                if (i > 2)
                    dtR(i) = 0;
                    dtRdot(i-1) = 0;
                    time(i) = 0;
                end
            end

            if (nargin > 15 && ~isempty(waitbar_handle))
                waitbar_handle.goTime(i);
            end
        end
        
        %-------------------------------------------------------------------------------------------------------------
        % RECEIVER CLOCK AND INTER-SYSTEM BIAS (IF ANY) ESTIMATION BY MULTI-EPOCH LEAST-SQUARES ADJUSTMENT: PROCESSING
        %-------------------------------------------------------------------------------------------------------------

        if (nisbs > 1)
            %unknown_index = find(~isnan(n_obs_epoch));
            n_obs_epoch(isnan(n_obs_epoch)) = [];
            epoch_index(isnan(epoch_index)) = [];

            index_nan = find(isnan(y0_all),1);
            y0_all(index_nan:end) = [];
            b_all(index_nan:end)  = [];

            n = length(y0_all);
            m = 3*npos + length(n_obs_epoch) + (nisbs-1);
            if (~isempty(index_nan))
                A_all = A_all(1:index_nan-1,1:size(A_all,2));
                Q_all = Q_all(1:index_nan-1,1);
            end
            A_all(isnan(A_all)) = 0;
            Q_all(isnan(Q_all)) = 0;

            %set the design matrix to estimate the receiver clock
            A_all(1:epoch_index(1),3*npos+1) = 1;
            for e = 2 : length(epoch_index)
                A_all(epoch_index(e-1)+1:epoch_index(e),3*npos+e) = 1;
            end

            %remove unavailable epochs from the unknowns (receiver clock)
            avail_index = any(A_all,1);
            avail_IFBs = [];
            if (cc.glo.isActive())
                %check if one of the GLONASS IFBs was removed from the unknowns
                avail_IFBs = avail_index(3*npos+nEpochs_reduced+1:3*npos+nEpochs_reduced+14);
            end
            avail_ISBs = [avail_IFBs avail_index(3*npos+nEpochs_reduced+length(avail_IFBs)+1:end)];
            A_all(:,~avail_index) = [];
            m = m - length(find(~avail_index));

            A  = A_all;
            y0 = y0_all;
            b  = b_all;
            Q  = diag(Q_all);

            %least-squares solution (broken down to improve computation speed)
            K = A';
            P = Q\A;
            N = K*P;
            Y = (y0-b);
            R = Q\Y;
            L = K*R;
            x = N\L;

            %variance of the estimation error
            y_hat = A*x + b;
            v_hat = y0 - y_hat;
            V = v_hat';
            T = Q\v_hat;
            sigma02_hat = (V*T)/(n-m);

            %estimated values
            % dtR = zeros(size(dtR));
            % dtR(avail_index(4:end-(nisbs-1))) = x(4:end-(nisbs-1))/v_light;

            ISBs(avail_ISBs) = x(end-(nisbs-2-sum(~avail_IFBs)):end)/v_light;
        end
        
        %----------------------------------------------------------------------------------------------
        % RECEIVER CLOCK DRIFT DISCONTINUITIES
        %----------------------------------------------------------------------------------------------

        %check if there is any discontinuity in the clock drift
        clock_thresh = 1e-5;
        dtRdot = zero2nan(dtRdot);
        disc = find(abs(dtRdot-mean(dtRdot(~isnan(dtRdot)))) > clock_thresh);

        %remove discontinuities from the clock drift
        for i = 1 : length(disc)
            if (disc(i) < 5)
                dtRdot(disc(i)) = median(dtRdot(1:10));
            elseif (disc(i) <= nEpochs-6)
                dtRdot(disc(i)) = median(dtRdot(disc(i)-4:disc(i)+5));
            elseif (disc(i) > nEpochs-6)
                dtRdot(disc(i)) = median(dtRdot(nEpochs-10:nEpochs-1));
            end
        end

        dtRdot(end+1) = dtRdot(end);

        % %check if it is needed to correct observations for receiver clocks offsets
        % % (some RINEX files contain clock-corrected observations, although they are
        % %  not respecting the specifications); clock offsets lower than 1
        % %  microsecond don't need to be corrected
        % if (max(abs(dtR)) < 1e-6)
        %     return
        % end

        %check which observations must be corrected for receiver clock offsets
        % (some receivers have inconsistent observations, e.g. code with clock
        %  jumps, phase without)

%         %jump detection threshold
%         j_thres = (clock_thresh*10)*v_light;
%
%         %flags
%         flag_jumps_pr1 = 0;
%         flag_jumps_pr2 = 0;
%         flag_jumps_ph1 = 0;
%         flag_jumps_ph2 = 0;
%
%         for i = 1 : length(disc)
%
%             for s = 1 : nSatTot
%
%                 %check code on L1
%                 if (pr1(s,disc(i):disc(i)+1) ~= 0)
%                     if (abs(diff(pr1(s,disc(i):disc(i)+1))) > j_thres)
%                         flag_jumps_pr1 = 1;
%                     end
%                 end
%
%                 %check code on L2
%                 if (pr2(s,disc(i):disc(i)+1) ~= 0)
%                     if (abs(diff(pr2(s,disc(i):disc(i)+1))) > j_thres)
%                         flag_jumps_pr2 = 1;
%                     end
%                 end
%
%                 %check phase on L1
%                 if (ph1(s,disc(i):disc(i)+1) ~= 0)
%                     if (abs(diff(ph1(s,disc(i):disc(i)+1)))*lambda(s,1) > j_thres)
%                         flag_jumps_ph1 = 1;
%                     end
%                 end
%
%                 %check phase on L2
%                 if (ph2(s,disc(i):disc(i)+1) ~= 0)
%                     if (abs(diff(ph2(s,disc(i):disc(i)+1)))*lambda(s,2) > j_thres)
%                         flag_jumps_ph2 = 1;
%                     end
%                 end
%
%                 %no need to go through all satellites
%                 if (any([flag_jumps_pr1 flag_jumps_pr2 flag_jumps_ph1 flag_jumps_ph2]))
%                     break
%                 end
%             end
%
%             %no need to go through all discontinuities
%             if (any([flag_jumps_pr1 flag_jumps_pr2 flag_jumps_ph1 flag_jumps_ph2]))
%                 break
%             end
%         end

        %-------------------------------------------------------------------------------------------------------------
        % OBSERVATION CORRECTION FOR CLOCK ERROR (v_light * dt)
        %-------------------------------------------------------------------------------------------------------------
        
        ph = zero2nan(bsxfun(@times, zero2nan([ph1; ph2]), [lambda(:, 1); lambda(:,2)])');
        % Filter dtR
        lim = getOutliers(dtR ~= 0);
        dtR = simpleFill1D(zero2nan(dtR), dtR == 0, 'spline');
        for i = 1 : size(lim, 1)
            dtR(lim(i,1) : lim(i,2)) = splinerMat([], dtR(lim(i,1) : lim(i,2)), 3);
        end
        ph = bsxfun(@minus, ph, v_light * dtR);
        
        % apply new dtR to pseudo-ranges
        pr = zero2nan([pr1; pr2]');
        pr = bsxfun(@minus, pr, v_light * dtR);

        % flag by high deviation of the 4th derivate
        sensor = Core_Utils.diffAndPred(zero2nan(ph), 4);
        sensor = abs(bsxfun(@minus, sensor, median(sensor, 2, 'omitnan')));
        flag = sensor > 20;
        ph1(flag(:,1:size(ph1,1))') = NaN;
        ph2(flag(:,(size(ph2,1)+1):end)') = NaN;

        ph = nan2zero(bsxfun(@rdivide, zero2nan(ph), [lambda(:, 1); lambda(:,2)]'));
        ph1 = zero2nan(ph(:,1:size(ph1,1))');
        ph2 = zero2nan(ph(:,(size(ph2,1)+1):end)');
        pr1 = (pr(:,1:size(pr1,1))');
        pr2 = (pr(:,(size(pr2,1)+1):end)');
        
        %----------------------------------------------------------------------------------------------
        % GEOMETRY FREE OBSERVABLES
        %----------------------------------------------------------------------------------------------

        ph_GF = compute_geometry_free(ph1, ph2, lambda, err_iono);
        
        %----------------------------------------------------------------------------------------------
        % SEARCH AND REJECT OUTLIERS ON THE GF
        %----------------------------------------------------------------------------------------------
        
        dGF = Core_Utils.diffAndPred(ph_GF', 3);
        flag = abs(dGF)' > 10 * perc((movstd(dGF(:), 30, 'omitnan')), 0.9);
        ph1(flag) = NaN;
        ph2(flag) = NaN;
        ph_GF(flag) = NaN;
        flag_ph = [flag' flag'] & flag_ph;
        
        %----------------------------------------------------------------------------------------------
        % WIDE LANE, NARROW LANE and MELBOURNE-WUBBENA OBSERVABLES
        %----------------------------------------------------------------------------------------------

        ph_MW = compute_melbourne_wubbena(ph1, ph2, pr1, pr2, lambda);

        ph1 = nan2zero(ph1);
        ph2 = nan2zero(ph2);
        pr1 = nan2zero(pr1);
        pr2 = nan2zero(pr2);
        ph_GF = nan2zero(ph_GF);
        ph_MW = nan2zero(ph_MW);
        
        %----------------------------------------------------------------------------------------------
        % OBSERVATION CORRECTION FOR CLOCK ERROR (time-shift)
        %----------------------------------------------------------------------------------------------

        %two types of corrections (as in http://www.navcen.uscg.gov/?pageName=RINEX):
        % 1. "frequency correction" c*dtR
        % 2. "receiver-satellite dynamics correction" by using Doppler if available,
        %     otherwise by interpolating observations on the time tag corrected by dtR

        % available epochs
        index_e = find(time ~= 0);

        % time "correction"
        time(index_e) = time(index_e) - dt_ph(index_e) - dtR(index_e) - time_desync(index_e);

        % variables to store interpolated observations
        pr1_interp = zeros(size(pr1));
        ph1_interp = zeros(size(ph1));
        pr2_interp = zeros(size(pr2));
        ph2_interp = zeros(size(ph2));

        freq1_required = (frequencies(1) == 1 || length(frequencies) > 1);
        freq2_required = (frequencies(1) == 2 || length(frequencies) > 1);

        for s = 1 : nSatTot

            if (any(pr1(s,:)) && freq1_required)

                index_s = find(pr1(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                pr1(s,index_x) = 0;

                if (length(index) > lagr_order)

    %                 if (any(dop1(s,index)))
    %                     corr = lambda(s,1).*dop1(s,index).*(time_desync(index) + dtR(index))';
    %                     pr1_interp(s,index) = pr1(s,index) - corr;
    %                 else
                        pr1_interp(s,index) = lagrange_interp1(time(index), pr1(s,index), time_ref(index), lagr_order);
    %                 end
                else
                    bad_sats(s,1) = 1;
                end
            end

            if (any(pr2(s,:)) && freq2_required)

                index_s = find(pr2(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                pr2(s,index_x) = 0;

                if (length(index) > lagr_order)

    %                 if (any(dop2(s,index)))
    %                     corr = lambda(s,2).*dop2(s,index).*(time_desync(index) + dtR(index))';
    %                     pr2_interp(s,index) = pr2(s,index) - corr;
    %                 else
                        pr2_interp(s,index) = lagrange_interp1(time(index), pr2(s,index), time_ref(index), lagr_order);
    %                 end
                else
                    bad_sats(s,1) = 1;
                end
            end

            if (any(ph1(s,:)) && freq1_required)
                index_s = find(ph1(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                ph1(s,index_x) = 0;

                if (length(index) > lagr_order)

                    %if (flag_jumps_ph1)
                        if (flag_doppler_cs && any(dop1(s,index)))
                            dop1(s,index) = dop1(s,index) + v_light*dtRdot(index)'/lambda(s,1);
                        end
                    %end
                end
            end

            if (any(ph2(s,:)) && freq2_required)
                index_s = find(ph2(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                ph2(s,index_x) = 0;

                if (length(index) > lagr_order)

                    %if (flag_jumps_ph2)
                        if (flag_doppler_cs && any(dop2(s,index)))
                            dop2(s,index) = dop2(s,index) + v_light*dtRdot(index)'/lambda(s,2);
                        end
                    %end
                end
            end

            if (any(ph1(s,:)) && freq1_required)

                index_s = find(ph1(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                ph1(s,index_x) = 0;

                if (length(index) > lagr_order)

                    [ph1(s,:), cs_found, cs_correction_i] = detect_and_fix_cycle_slips(time, pr1(s,:), ph1(s,:), pr2(s,:), ph2(s,:), ph_GF(s,:), ph_MW(s,:), dop1(s,:), el(s,:), err_iono(s,:), lambda(s,1), lambda(s,2));

                    if ~isempty(cs_correction_i)
                        cs_correction_i(:,1) = s;
                        cs_correction_i(:,2) = 1;
                        status_cs=[status_cs;cs_correction_i];
                    end

                    index_s = find(ph1(s,:) ~= 0);
                    index = intersect(index_e,index_s);

                    index_x = setdiff(index_s, index_e);
                    ph1(s,index_x) = 0;

    %                 if (any(dop1(s,index)))
    %                     corr = dop1(s,index).*(time_desync(index) + dtR(index))';
    %                     ph1_interp(s,index) = ph1(s,index) - corr;
    %                 else
                        ph1_interp(s,index) = lagrange_interp1(time(index), ph1(s,index), time_ref(index), lagr_order);
    %                 end

                    if (exist('cs_found', 'var') && cs_found)
                        fprintf('Pre-processing: %d cycle-slip(s) detected on L1 for satellite %02d\n', cs_found, s);
                    end
                else
                    bad_sats(s,1) = 1;
                end

            elseif (any(pr1(s,:)) && freq1_required)

                bad_sats(s,1) = 1;
            end

            if (any(ph2(s,:)) && freq2_required)

                index_s = find(ph2(s,:) ~= 0);
                index = intersect(index_e,index_s);

                index_x = setdiff(index_s, index_e);
                ph2(s,index_x) = 0;

                if (length(index) > lagr_order)

                    [ph2(s,:), cs_found, cs_correction_i] = detect_and_fix_cycle_slips(time, pr2(s,:), ph2(s,:), pr1(s,:), ph1(s,:), ph_GF(s,:), ph_MW(s,:), dop2(s,:), el(s,:), err_iono(s,:), lambda(s,2), lambda(s,1));

                    if ~isempty(cs_correction_i)
                        cs_correction_i(:,1) = s;
                        cs_correction_i(:,2) = 2;
                        status_cs=[status_cs;cs_correction_i];
                    end

                    index_s = find(ph2(s,:) ~= 0);
                    index = intersect(index_e,index_s);

                    index_x = setdiff(index_s, index_e);
                    ph2(s,index_x) = 0;

    %                 if (any(dop2(s,index)))
    %                     corr = dop2(s,index).*(time_desync(index) + dtR(index))';
    %                     ph2_interp(s,index) = ph2(s,index) - corr;
    %                 else
                        ph2_interp(s,index) = lagrange_interp1(time(index), ph2(s,index), time_ref(index), lagr_order);
    %                 end

                    if (exist('cs_found', 'var') && cs_found)
                        fprintf('Pre-processing: %d cycle-slip(s) detected on L2 for satellite %02d\n', cs_found, s);
                    end
                else
                    bad_sats(s,1) = 1;
                end

            elseif (any(pr2(s,:)) && freq2_required)

                bad_sats(s,1) = 1;
            end
        end

    %    for s = 1 : nSatTot
    %         if (freq1_required)
    %             if (any(ph1(s,:)))
    %                 [pr1(s,:)] = code_smoother(pr1(s,:), ph1(s,:), lambda(s,1), lagr_order);
    %             else
    %                 pr1(s,:) = 0;
    %             end
    %         end
    %         if (freq2_required)
    %             if (any(ph2(s,:)))
    %                 [pr2(s,:)] = code_smoother(pr2(s,:), ph2(s,:), lambda(s,2), lagr_order);
    %             else
    %                 pr2(s,:) = 0;
    %             end
    %         end
    %    end

        % flag by high deviation of the 4th derivate
        ph = zero2nan(bsxfun(@times, [ph1_interp; ph2_interp], [lambda(:, 1); lambda(:,2)])');
        pr = zero2nan([pr1_interp; pr2_interp]');
        sensor = Core_Utils.diffAndPred(ph, 4);
        sensor = abs(bsxfun(@minus, sensor, median(sensor, 2, 'omitnan')));
        flag = sensor > 3;
        flag_ph = flag | flag_ph;
        
        [ph, flag] = Core_PP.flagRawObsD4(ph, [], [], 6);
        flag_ph = flag | flag_ph;
        [pr, flag] = Core_PP.flagRawObsD4(pr, [], [], 6);
        flag_pr = flag | flag_pr;
        % flagging borders could decrease the performance of the final solution
        % skipping
        %[ph, flag] = Core_Pre_Processing.flagBorders(ph, 20);
        %flag_ph = flag | flag_ph;
        %[pr, flag] = Core_Pre_Processing.flagBorders(pr, 20);
        %flag_pr = flag | flag_pr;
        
        ph = nan2zero(bsxfun(@rdivide, zero2nan(ph), [lambda(:, 1); lambda(:,2)]'));
        ph1 = (ph(:,1:size(ph1,1))');
        ph2 = (ph(:,(size(ph2,1)+1):end)');
        pr1 = nan2zero(pr(:,1:size(pr1,1))');
        pr2 = nan2zero(pr(:,(size(pr2,1)+1):end)');
                
        for s = 1 : nSatTot
            %repeat remove short arcs after cycle slip detection
            % remove short arcs
            min_arc = max([state.getMinArc() lagr_order]);
            pr1(s,:) = remove_short_arcs(pr1(s,:), min_arc);
            pr2(s,:) = remove_short_arcs(pr2(s,:), min_arc);
            ph1(s,:) = remove_short_arcs(ph1(s,:), min_arc);
            ph2(s,:) = remove_short_arcs(ph2(s,:), min_arc);
        end
    end
    % %flag epochs with 4 or more slipped satellites as "bad"
    % [num_cs_occur, epoch] = hist(status_cs(:,3),unique(status_cs(:,3)));
    % idx_cs_occur = num_cs_occur >= 4;
    % bad_epochs(epoch(idx_cs_occur)) = 1;
    
    ph1 = jmpFix(ph1, lambda(:,1), state);
    ph2 = jmpFix(ph2, lambda(:,2), state);
    %%% DEBBUG data = ph1'; figure; subplot(2,2,1); plot(zero2nan(data)); subplot(2,2,2); plot(diff(zero2nan(data),1));subplot(2,2,3); plot(diff(zero2nan(data),2));subplot(2,2,4); plot(bsxfun(@minus,diff(zero2nan(data),4),median(diff(zero2nan(data),4),2,'omitnan'))); subplot(2,2,1); title('data'); subplot(2,2,2); title('first derivative'); subplot(2,2,3); title('second derivative'); subplot(2,2,4); title('forth derivative (reduced)');
    
    % At this point the data is syncronized to the reference time, and corrected for dtR and de-sync
    % resetting time, dtR, dtRdot
    time = time_ref;
    % for debug reason the dtR is not reset
    dtR = dtR + dt_ph;
    %dtRdot = dtRdot * 0;
end

function [ph_main, cs_correction_count, cs_correction_i] = detect_and_fix_cycle_slips(time, pr_main, ph_main, pr_sec, ph_sec, ph_GF, ph_MW, dop, el, err_iono, lambda_main, lambda_sec)

    global cutoff cs_threshold_preprocessing
    cs_resolution = 1;

    flag_plot = 0;
    flag_doppler_cs = 1;

    cs_correction_i=[];
    cs_correction_count = 0;

    cutoff_idx = find(el(1,:) > cutoff);
    avail_idx = find(ph_main(1,:) ~= 0);
    idx = intersect(cutoff_idx, avail_idx);

    if (length(idx) > 3)
        [p,~,mu] = polyfit(idx,err_iono(1,idx),3);
        err_iono_fit = polyval(p,idx,[],mu);
    else
        err_iono_fit = err_iono(1,idx);
    end
    % if (flag_plot)
    %     figure
    %     plot(idx,err_iono(1,idx),'r.')
    %     hold on
    %     plot(idx,err_iono_fit,'g.')
    % end

    N_mat = zeros(size(ph_main));
    N_mat(1,idx) = (pr_main(1,idx) - lambda_main(1,1).*ph_main(1,idx) - 2.*err_iono_fit)./lambda_main(1,1);

    if (~isempty(N_mat(1,N_mat(1,:)~=0)))

        interval = diff(time);
        interval = [interval; interval(end)]';

        %initialization
        jmp_code    = (1:(length(ph_main)-1))';
        jmp_doppler = (1:(length(ph_main)-1))';
        jmp_deriv   = (1:(length(ph_main)-1))';
        jmp_GF      = (1:(length(ph_main)-1))';
        jmp_MW      = (1:(length(ph_main)-1))';

        %detection (code)
        N_mat(1,N_mat(1,:)==0) = NaN;
    %     delta_N = diff(N_mat(1,:));
    %     delta_code = (diff(N_mat(1,:))./interval(1:end-1))';
        delta_code = diff(N_mat(1,:))';

        delta_test = delta_code;
        not_zero = find(delta_test ~= 0);
        not_nan  = find(~isnan(delta_test));
        avail_code = intersect(not_zero, not_nan);

        if (~isempty(delta_test(avail_code)))
            outliers = batch_outlier_detection(delta_test(avail_code),median(round(interval)));
            [~,jmp_code] = intersect(delta_test,outliers);
        end

        %detection (Doppler)
        delta_doppler = zeros(size(delta_code));
        if (flag_doppler_cs && any(dop) && (sum(~~dop) == sum(~~ph_main)))

    %         pred_phase = ph(1:end-1) + ((dop(2:end)+dop(1:end-1))/2).*interval(1:end-1);
    %         pred_phase = ph(1:end-1) - ((dop(2:end)+dop(1:end-1))/2).*interval(1:end-1);
            pred_phase = ph_main(1:end-1) - (dop(2:end)+dop(1:end-1))/2;

    %         delta_doppler_all = ((pred_phase - ph(2:end))./interval(1:end-1))';
            delta_doppler_all = (pred_phase - ph_main(2:end))';
            delta_doppler_all(ph_main(2:end)==0) = 0;
            delta_doppler_all(delta_doppler_all == 0) = NaN;
    %         delta_doppler_all = diff(delta_doppler_all)'./interval(1:end-2);
            delta_doppler_all = diff(delta_doppler_all)';
            pos1 = find(idx~=length(N_mat));
            pos2 = find(idx~=length(N_mat)-1);
            pos = intersect(pos1,pos2);
            delta_doppler(idx(pos)) = delta_doppler_all(idx(pos));
            delta_doppler(delta_doppler == 0) = NaN;
        end

        delta_test = delta_doppler;
        not_zero = (delta_test ~= 0);
        not_nan  = (~isnan(delta_test));
        avail_doppler = find(not_zero & not_nan);

        if (~isempty(delta_test(avail_doppler)))
            outliers = batch_outlier_detection(delta_test(avail_doppler),median(round(interval)));
            [~,jmp_doppler] = intersect(delta_test,outliers);
        end

        idx_interp = setdiff(avail_doppler,jmp_doppler);
        if (length(idx_interp) > 1)
            p = polyfit(idx_interp,delta_doppler(idx_interp),1);
            delta_doppler(idx_interp) = delta_doppler(idx_interp) - polyval(p,idx_interp);
        end

        jmp_doppler = sort(jmp_doppler);
        jmp_doppler(diff(jmp_doppler) == 1) = [];

        %detection (phase 3rd order derivative)
        delta_deriv = zeros(size(delta_code));
        ph_tmp = ph_main;
        ph_tmp(1,ph_tmp(1,:)==0) = NaN;
    %     delta_deriv_all = (diff(diff(diff(ph_tmp(1,:))))./interval(1:end-3).^3)';
        delta_deriv_all = diff(diff(diff(ph_tmp(1,:))))';
        pos1 = (idx~=length(N_mat));
        pos2 = (idx~=length(N_mat)-1);
        pos3 = (idx~=length(N_mat)-2);
        pos = pos1 & pos2 & pos3;
        delta_deriv(idx(pos)) = delta_deriv_all(idx(pos));
        delta_deriv(delta_deriv == 0) = NaN;

        delta_test = delta_deriv;
        not_zero = (delta_test ~= 0);
        not_nan  = ~isnan(delta_test);
        avail_deriv = not_zero & not_nan;

        if (~isempty(delta_test(avail_deriv)))
            outliers = batch_outlier_detection(delta_test(avail_deriv),median(round(interval)));
            [~,jmp_deriv] = intersect(delta_test,outliers);
        end
        jmp_deriv = sort(jmp_deriv);
        jmp_deriv(diff(jmp_deriv) == 1) = [];

        %detection (geometry free)
        delta_GF = zeros(size(delta_code));
        ph_GF(1,ph_GF(1,:)==0) = NaN;
    %     delta_GF_all = (diff(ph_GF)./interval(1:end-1))';
        delta_GF_all_ref = diff(ph_GF)';
        delta_GF_all = diff(diff(ph_GF))';
        pos1 = (idx~=length(N_mat));
        pos2 = (idx~=length(N_mat)-1);
        pos = pos1 & pos2;
        delta_GF_ref(idx(pos1)) = delta_GF_all_ref(idx(pos1));
        delta_GF(idx(pos)) = delta_GF_all(idx(pos));
        delta_GF(delta_GF == 0) = NaN;

        delta_test = delta_GF;
        not_zero = (delta_test ~= 0);
        not_nan  = (~isnan(delta_test));
        avail_GF = not_zero & not_nan;

        if (~isempty(delta_test(avail_GF)))
            outliers = batch_outlier_detection(delta_test(avail_GF),median(round(interval)));
            outliers(abs(outliers) < 0.04) = [];
            [~,jmp_GF] = intersect(delta_test,outliers);
        end
        jmp_GF = sort(jmp_GF);
        jmp_GF(diff(jmp_GF) == 1) = [];

        %detection (Melbourne-Wubbena)
        delta_MW = zeros(size(delta_code));
        ph_MW(1,ph_MW(1,:)==0) = NaN;
    %     delta_MW_all = (diff(ph_MW)./interval(1:end-1))';
        delta_MW_all = diff(ph_MW)';
        pos = (idx~=length(N_mat));
        delta_MW(idx(pos)) = delta_MW_all(idx(pos));

        delta_test = delta_MW;
        not_zero = find(delta_test ~= 0);
        not_nan  = find(~isnan(delta_test));
        avail_MW = intersect(not_zero, not_nan);

        if (~isempty(delta_test(avail_MW)))
            outliers = batch_outlier_detection(delta_test(avail_MW),median(round(interval)));
            outliers(abs(outliers) < 0.1) = [];
            [~,jmp_MW] = intersect(delta_test,outliers);
        end

        %select two observables with low standard deviation
        [min_std_code]    = detect_minimum_std(delta_code(avail_code));
        [min_std_doppler] = detect_minimum_std(delta_doppler(avail_doppler));
        [min_std_deriv]   = detect_minimum_std(delta_deriv(avail_deriv));
        [min_std_GF]      = detect_minimum_std(delta_GF(avail_GF));
        [min_std_MW]      = detect_minimum_std(delta_MW(avail_MW));

        min_stds = [min_std_code min_std_doppler min_std_deriv min_std_GF min_std_MW];
        jmps = {jmp_code; jmp_doppler; jmp_deriv; jmp_GF; jmp_MW};

        [~, pos1] = min(min_stds); min_stds(pos1) = 1e30;
        [~, pos2] = min(min_stds);

        %consider cycle slips detected by either geometry-free or Melbourne-Wubbena
        if ((pos1 == 4 && pos2 == 5) || (pos1 == 5 && pos2 == 4))
            jmp = sort(union(jmps{pos1},jmps{pos2}));
        else
            jmp = sort(intersect(jmps{pos1},jmps{pos2}));
        end

        %compute the dataset from which to extract the potential cs correction
        if (any(delta_GF_ref))
            delta = delta_GF_ref/lambda_main;
        elseif (any(delta_MW))
            freq_main = Core_Utils.V_LIGHT ./ lambda_main;
            freq_sec = Core_Utils.V_LIGHT ./ lambda_sec;
            delta = delta_MW*abs((freq_main - freq_sec)/(freq_main*lambda_main));
        elseif (any(delta_doppler))
            delta = -delta_doppler;
        else
            delta = -delta_deriv;
        end

        %ignore cycle slips smaller than cs_threshold_preprocessing
        jmp(roundmod(abs(delta(jmp)), cs_resolution) < cs_threshold_preprocessing) = [];

        %ignore cycle slips that cannot be fixed
        jmp(isnan(delta(jmp))) = [];

        %exclude observation epochs with subsequent cycle slips
        idx_bad_obs1 = find(diff(jmp) == 1);
        idx_bad_obs1 = unique([idx_bad_obs1; idx_bad_obs1+1]);

        %exclude observation epochs with computed cycle slip corrections "far" from integer values
        %thres = 0.1;   %ENABLED
        thres = 1e-10; %DISABLED (i.e. just exclude all)
        idx_bad_obs2 = find(abs(delta(jmp) - roundmod(delta(jmp), cs_resolution)) > thres);

        idx_bad_obs = union(idx_bad_obs1, idx_bad_obs2);
        if (~isempty(idx_bad_obs))
            ph_main(jmp(idx_bad_obs)) = 0;
            jmp_bad_obs = jmp(idx_bad_obs);
            for j = 1 : length(jmp_bad_obs)
                cs_correction_count = cs_correction_count + 1;
                cs_correction_i(cs_correction_count,3) = jmp_bad_obs(j); %#ok<*AGROW>
                cs_correction_i(cs_correction_count,4) = 0;
                cs_correction_i(cs_correction_count,5) = delta(jmp_bad_obs(j));
                cs_correction_i(cs_correction_count,6) = 0;
            end
            jmp(idx_bad_obs) = [];
        end

        %cycle slips detected by code and doppler observables must be of the same sign
        if ((pos1 == 1 && pos2 == 2) || (pos1 == 2 && pos2 == 1))
            sign_OK = (sign(delta_code(jmp)) .* sign(delta_doppler(jmp)) > 0);
            jmp = jmp(sign_OK);
        end

        %cycle slips detected by code and derivative observables must be of the same sign
        if ((pos1 == 1 && pos2 == 3) || (pos1 == 3 && pos2 == 1))
            sign_OK = (sign(delta_code(jmp)) .* sign(delta_deriv(jmp)) > 0);
            jmp = jmp(sign_OK);
        end

        %cycle slips detected by doppler and derivative observables must be of the same sign
        if ((pos1 == 2 && pos2 == 3) || (pos1 == 3 && pos2 == 2))
            sign_OK = (sign(delta_doppler(jmp)) .* sign(delta_deriv(jmp)) > 0);
            jmp = jmp(sign_OK);
        end

        if (isempty(jmp))
            return
        end

        if (flag_plot)
            figure; hold on %#ok<UNRCH>
            if (pos1 == 1 || pos2 == 1)
                plot(delta_code)
                plot(jmp,delta_code(jmp),'mx')
            end
            if (pos1 == 2 || pos2 == 2)
                plot(delta_doppler,'g')
                plot(jmp,delta_doppler(jmp),'rx')
            end
            if (pos1 == 3 || pos2 == 3)
                plot(delta_deriv,'c')
                plot(jmp,delta_deriv(jmp),'yx')
            end
            if (pos1 == 4 || pos2 == 4)
                plot(delta_GF,'k--')
                plot(jmp,delta_GF(jmp),'ko')
            end
            if (pos1 == 5 || pos2 == 5)
                plot(delta_MW,'k--')
                plot(jmp,delta_MW(jmp),'ko')
            end

            figure
            idx_plot = find(N_mat~=0);
            plot(idx_plot, N_mat(idx_plot));

            coltab = colorcube(2*length(jmp));
            c = 1;
        end

        %fixing
        for j = 1 : length(jmp)
            if (jmp(j) <= 1 || jmp(j) >= length(ph_main)-1)
                check = 0;
                pos_zeros = [];
            else
                pos_zeros = find(N_mat(1,jmp(j)-1:jmp(j)+1) == 0,1,'last');
                check = isempty(pos_zeros);
            end
            if (ismember(jmp(j),cutoff_idx))% || N_before_zero ~= 0)
                if (check)
                    idx_zeros = (ph_main(1,:)==0);
                    cs_correction = roundmod(delta(jmp(j)), cs_resolution);
                    cs_correction_count = cs_correction_count + 1;
                    cs_correction_i(cs_correction_count,3)=jmp(j); %#ok<*AGROW>
                    cs_correction_i(cs_correction_count,4)=cs_correction;
                    cs_correction_i(cs_correction_count,5)=delta(jmp(j));
                    cs_correction_i(cs_correction_count,6)=1;

                    if ((pos1 == 4 && pos2 == 5) || (pos1 == 5 && pos2 == 4)) %if detected only by GF and MW
                        ph_temp = ph_main;
                        ph_temp(1,jmp(j)+1:end) = ph_main(1,jmp(j)+1:end) + cs_correction;
                        e = err_iono;
                        e(idx) = err_iono_fit;
                        ph_GF_new = compute_geometry_free(ph_temp, ph_sec, [lambda_main, lambda_sec], e);
                        ph_MW_new = compute_melbourne_wubbena(ph_temp, ph_sec, pr_main, pr_sec, [lambda_main, lambda_sec]);
                        %if (abs(ph_MW_new(1,jmp(j)+1)-ph_MW_new(1,jmp(j))) > abs(ph_MW(1,jmp(j)+1)-ph_MW(1,jmp(j))) && ...
                        if (abs(ph_GF_new(1,jmp(j)+1)-ph_GF_new(1,jmp(j))) > abs(ph_GF(1,jmp(j)+1)-ph_GF(1,jmp(j)))) %correction applied to the wrong frequency
                            cs_correction_i = [];
                            cs_correction_count = cs_correction_count - 1;
                            continue
                        end
                    end

                    %apply correction
                    N_mat(1,jmp(j)+1:end) = N_mat(1,jmp(j)+1:end) - cs_correction;
                    ph_main(1,jmp(j)+1:end) = ph_main(1,jmp(j)+1:end) + cs_correction;
                    N_mat(1,idx_zeros) = 0;
                    ph_main(1,idx_zeros) = 0;

                    if (flag_plot)
                        hold on %#ok<UNRCH>
                        idx_plot = find(N_mat~=0);
                        plot(idx_plot, N_mat(idx_plot),'Color',coltab(2*c-1,:));

                        c = c + 1;
                    end

                elseif (pos_zeros == 3)

    %                 N_before_zero = N_mat(1,jmp(j));

                elseif (pos_zeros == 2)

    %                 idx = (ph(1,:)==0);
    %                 cs_correction = round(N_mat(1,jmp(j)+1) - N_before_zero);
    %                 N_mat(1,jmp(j)+1:end) = N_mat(1,jmp(j)+1:end) - cs_correction;
    %                 ph   (1,jmp(j)+1:end) = ph   (1,jmp(j)+1:end) + cs_correction;
    %                 N_mat(1,idx) = 0;
    %                 ph   (1,idx) = 0;
    %
    %                 if (flag_plot)
    %                     hold on
    %                     idx_plot = find(N_mat~=0);
    %                     plot(idx_plot, N_mat(idx_plot),'Color',coltab(2*c-1,:));
    %
    %                     c = c + 1;
    %                 end
    %
    %                 N_before_zero = 0;
    %                 N_after_zero = N_mat(1,jmp(j)+1);

                elseif (pos_zeros == 1)

    %                 idx = (ph(1,:)==0);
    %                 cs_correction = round(N_mat(1,jmp(j)+1) - N_after_zero);
    %                 N_mat(1,jmp(j)+1:end) = N_mat(1,jmp(j)+1:end) - cs_correction;
    %                 ph   (1,jmp(j)+1:end) = ph   (1,jmp(j)+1:end) + cs_correction;
    %                 N_mat(1,idx) = 0;
    %                 ph   (1,idx) = 0;
    %
    %                 if (flag_plot)
    %                     hold on
    %                     idx_plot = find(N_mat~=0);
    %                     plot(idx_plot, N_mat(idx_plot),'Color',coltab(2*c-1,:));
    %
    %                     c = c + 1;
    %                 end
    %
    %                 N_after_zero = 0;
                end
            end
        end
    end
end

function [yi] = lagrange_interp1(x,y,xi,n)
    d = n/2;
    yi = zeros(size(y));
    for t = 1 : length(xi)
        if (t<=d)
            yi(t) = LagrangeInter(x(1:n)', y(1:n), xi(t));
        elseif (t>(length(x)-d))
            yi(t) = LagrangeInter(x(end-d-1:end)', y(end-d-1:end), xi(t));
        else
            yi(t) = LagrangeInter(x(t-d:t+d)', y(t-d:t+d), xi(t));
        end
    end
end

function [min_std] = detect_minimum_std(time_series)
    min_std = 1e30;
    d = 5; %half window size
    if (length(time_series) < 2*d)
        return
    end
    if verLessThan('matlab', '9.0.0')
        %  explicit implementation of movstd
        mov_std = zeros(size(time_series));
        for t = 1 : length(time_series)
            if (t<=d)
                mov_std(t) = std(time_series(1:d));
            elseif (t>(length(time_series)-d))
                if (length(time_series)-d-1 == 0)
                    mov_std(t) = mov_std(t-1);
                else
                    mov_std(t) = std(time_series(end-d-1:end));
                end
            else
                mov_std(t) = std(time_series(t-d:t+d));
            end
        end
    else
        mov_std = movstd(time_series, 11);
    end
    if (~isempty(mov_std))
        min_std = min(mov_std);
    end
end

function [pr] = code_smoother(pr, ph, lambda, order)
    pr(pr == 0) = NaN;
    ph(ph == 0) = NaN;
    N = (pr - lambda*ph) /  lambda;
    % N = smoothing(N, order, 1);
    N_smar = smartFilter(N, order);
    % figure; plot(N,'.-'); hold on; plot(N_smar,'g');
    %N_smar = N;
    pr = lambda*(ph + N_smar);
    pr(isnan(pr)) = 0;
    ph(isnan(ph)) = 0;
end

function [ph] = jmpFix(ph, lambda, state)
    % after a missing data try to fix the jmp in the observations with an N * lambda correction
    cs_factor = state.cs_thr;
    ph = ph';
    for a = 1 : size(ph,2)
        arc = zero2nan(ph(:,a));
        x_arc = (1 : numel(arc))';
        arc_fit = zeros(size(arc));
        
        int_factor = zeros(size(arc));
        int_factor(~isnan(arc)) = lambda(a,1) * cs_factor;
        arc = arc ./ int_factor;
        if any(arc)
            %PLOT:figure;
            %PLOT:plot(x_arc, arc, 'LineWidth', 1); hold on;
            lim = getOutliers(~isnan(arc));
            for l = 1 : size(lim, 1) - 1
                id = lim(1, 1) : lim(l, 2);
                id = id(~isnan(arc(id)));
                id = id((max(1, numel(id) - state.getMinArc+1) : end));
                
                p = polyfit(x_arc(id), arc(id), 3);
                if l < size(lim, 1) && (lim(l+1, 1) - lim(l, 2)) < 100
                    arc_fit1 = polyval(p, x_arc([lim(l, 2) lim(l + 1, 1)]));
                    
                    id2 = lim(l + 1, 1) : lim(l + 1, 2);
                    id2 = id2(1 : min(numel(id2), state.getMinArc));
                    p2 = polyfit(x_arc(id2), arc(id2), 3);
                    arc_fit2 = polyval(p2, x_arc([lim(l, 2) lim(l + 1, 1)]));
                    
                    arc(lim(l + 1, 1) : end) = arc(lim(l + 1, 1) : end) + round(mean(arc_fit1 - arc_fit2));
                end
                %PLOT: plot(x_arc, arc);
            end
            %PLOT: plot(x_arc, arc, 'LineWidth', 2); hold on;
        end
        ph(:,a) = nan2zero(arc) .* int_factor;
    end
    ph = ph';

    %PLOT: figure; plot(zero2nan(ph));
    %PLOT: d_ph = diff(zero2nan([ph(1,:); ph])); figure; plot(d_ph - movmedian(d_ph, 3, 'omitnan'));
end
