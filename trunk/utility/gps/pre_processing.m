function [pr1, ph1, pr2, ph2, dtR, dtRdot, bad_sats, bad_epochs, var_dtR, var_SPP, status_obs, status_cs] = pre_processing(time_ref, time, XR0, pr1, ph1, pr2, ph2, dop1, dop2, snr1, Eph, SP3, iono, lambda, nSatTot, waitbar_handle, flag_XR, sbas)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dtR, dtRdot, bad_sats, bad_epochs] = pre_processing(time_ref, time, XR0, pr1, ph1, pr2, ph2, dop1, dop2, snr1, Eph, SP3, iono, lambda, nSatTot, waitbar_handle);
%
% INPUT:
%   time_ref = GPS reference time
%   time     = GPS nominal time (as read from RINEX file)
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
%   lambda  = wavelength matrix (depending on the enabled constellations)
%   nSatTot = maximum number of satellites (given the enabled constellations)
%   waitbar_handle = handle to the waitbar object
%   flag_XR = 2: coordinate fixed to XR0 values
%           = 1: approximate coordinates available
%           = 0: no apriori coordinates available

% OUTPUT:
%   pr1 = processed code observation (L1 carrier)
%   ph1 = processed phase observation (L1 carrier)
%   pr2 = processed code observation (L2 carrier)
%   ph2 = processed phase observation (L2 carrier)
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%   bad_sats = vector for flagging "bad" satellites (e.g. too few observations, code without phase, etc)
%   bad_epochs  = vector with 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   var_SPP   = [code single point positioning a posteriori sigma, sum of
%                weighted squared residuals, redoundancy], one row per epoch
%   status_obs = for each satellite. NaN: not observed, 0: observed only, 1: used, -1: outlier 
%   status_cs = [satellite_number, frequency, epoch, cs_correction_fix, cs_correction_float, cs_corrected(0:no, 1:yes)]: vector with cs information
%
% DESCRIPTION:
%   Pre-processing of code and phase observations to correct them for
%    the receiver clock error.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Stefano Caldera
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

global cutoff snr_threshold n_sys flag_doppler_cs

v_light = goGNSS.V_LIGHT;

%number of epochs
nEpochs = length(time);

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

%vector for flagging "bad" satellites (e.g. too few observations, code without phase, etc)
bad_sats = zeros(nSatTot,1);

%vector with bad epochs definition
bad_epochs=NaN(nEpochs,1);

% vector with SPP a posteriori sigma
var_SPP=NaN(nEpochs,3);



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
cond_num = zeros(nEpochs,1);
cov_XR = zeros(3,3,nEpochs);
var_dtR = NaN(nEpochs,1);
status_obs = NaN(nSatTot,nEpochs);
status_cs=[];


for i = 1 : nEpochs
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE AND EPHEMERIS SELECTION
    %--------------------------------------------------------------------------------------------
    
    sat0 = find(pr1(:,i) ~= 0);
    status_obs(sat0,i) = 0; % satellite observed

    Eph_t = rt_find_eph (Eph, time(i), nSatTot);
    sbas_t = find_sbas(sbas, i);
    
    %----------------------------------------------------------------------------------------------
    % RECEIVER POSITION AND CLOCK ERROR
    %----------------------------------------------------------------------------------------------
    
    min_nsat_LS = 3 + n_sys;
    
    if (length(sat0) >= min_nsat_LS)
        [~, dtR_tmp, ~, ~, ~, ~, ~, ~, err_iono_tmp, sat, el_tmp, ~, ~, ~, cov_XR_tmp, var_dtR_tmp, ~, ~, ~, cond_num_tmp, bad_sat_i, bad_epochs(i), var_SPP(i,:)] = init_positioning(time(i), pr1(sat0,i), snr1(sat0,i), Eph_t, SP3, iono, sbas_t, XR0, [], [], sat0, [], lambda(sat0,:), cutoff, snr_threshold, 1, flag_XR, 0, 1);
        
        if isempty(var_dtR_tmp)
            var_dtR_tmp=NaN;
        end
        
        status_obs(sat,i) = 1; % satellite used
        status_obs(find(bad_sat_i==1),i)=-1; % satellite outlier
        
        if (~isempty(dtR_tmp) && ~isempty(sat))
            dtR(i) = dtR_tmp;
            err_iono(sat,i) = err_iono_tmp;
            el(sat,i) = el_tmp;
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
                dtR(i) = dtR(i-1) + (dtR(i-1) - dtR(i-2));
                dtRdot(i-1) = (dtR(i) - dtR(i-1))/(time(i) - time(i-1));
            end
        end
    else
        if (i > 2)
            dtR(i) = dtR(i-1) + (dtR(i-1) - dtR(i-2));
            dtRdot(i-1) = (dtR(i) - dtR(i-1))/(time(i) - time(i-1));
        end
    end
    
    if (nargin > 15 && ~isempty(waitbar_handle))
        waitbar_handle.goTime(i);
    end
end

%bad_epochs = (var_dtR(:,1) > 5e-6);

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK DRIFT DISCONTINUITIES
%----------------------------------------------------------------------------------------------

%check if there is any discontinuity in the clock drift
clock_thresh = 1e-5;
disc = find(abs(dtRdot-mean(dtRdot)) > clock_thresh);

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

%jump detection threshold
j_thres = (clock_thresh*10)*v_light;

%flags
flag_jumps_pr1 = 0;
flag_jumps_pr2 = 0;
flag_jumps_ph1 = 0;
flag_jumps_ph2 = 0;

for i = 1 : length(disc)
    
    for s = 1 : nSatTot
        
        %check code on L1
        if (pr1(s,disc(i):disc(i)+1) ~= 0)
            if (abs(diff(pr1(s,disc(i):disc(i)+1))) > j_thres)
                flag_jumps_pr1 = 1;
            end
        end
        
        %check code on L2
        if (pr2(s,disc(i):disc(i)+1) ~= 0)
            if (abs(diff(pr2(s,disc(i):disc(i)+1))) > j_thres)
                flag_jumps_pr2 = 1;
            end
        end
        
        %check phase on L1
        if (ph1(s,disc(i):disc(i)+1) ~= 0)
            if (abs(diff(ph1(s,disc(i):disc(i)+1)))*lambda(s,1) > j_thres)
                flag_jumps_ph1 = 1;
            end
        end
        
        %check phase on L2
        if (ph2(s,disc(i):disc(i)+1) ~= 0)
            if (abs(diff(ph2(s,disc(i):disc(i)+1)))*lambda(s,2) > j_thres)
                flag_jumps_ph2 = 1;
            end
        end
        
        %no need to go through all satellites
        if (any([flag_jumps_pr1 flag_jumps_pr2 flag_jumps_ph1 flag_jumps_ph2]))
            break
        end
    end
    
    %no need to go through all discontinuities
    if (any([flag_jumps_pr1 flag_jumps_pr2 flag_jumps_ph1 flag_jumps_ph2]))
        break
    end
end

%----------------------------------------------------------------------------------------------
% GEOMETRY FREE OBSERVABLE
%----------------------------------------------------------------------------------------------
ph_GF = zeros(size(ph1));
for s = 1 : nSatTot
    if (any(ph1(s,:)) && any(ph2(s,:)))
        
        index_1 = find(ph1(s,:) ~= 0);
        index_2 = find(ph2(s,:) ~= 0);
        index_3 = find(err_iono(s,:) ~= 0);
        index = intersect(index_1, index_2);
        index = intersect(index,   index_3);

        ph_GF(s,index) = (lambda(s,1)*ph1(s,index) - lambda(s,2)*ph2(s,index)) - ((goGNSS.F1^2-goGNSS.F2^2)/goGNSS.F2^2)*err_iono(s,index);
    end
end

%----------------------------------------------------------------------------------------------
% OBSERVATION CORRECTION FOR CLOCK ERROR
%----------------------------------------------------------------------------------------------

%two types of corrections (as in http://www.navcen.uscg.gov/?pageName=RINEX):
% 1. "frequency correction" c*dtR
% 2. "receiver-satellite dynamics correction" by using Doppler if available,
%     otherwise by interpolating observations on the time tag corrected by dtR

%available epochs
index_e = find(time ~= 0);

%nominal time desynchronization (e.g. with some low-cost receivers)
time_desync = time_ref - time;

%reference time "correction"
time_ref(index_e) = time(index_e) + dtR(index_e) + time_desync(index_e);

%variables to store interpolated observations
pr1_interp = zeros(size(pr1));
ph1_interp = zeros(size(ph1));
pr2_interp = zeros(size(pr2));
ph2_interp = zeros(size(ph2));

%Lagrange interpolation order
lagr_order = 10;

for s = 1 : nSatTot

    if (any(pr1(s,:)))

        index_s = find(pr1(s,:) ~= 0);
        index = intersect(index_e,index_s);
        
        if (length(index) > lagr_order)
            
            if (flag_jumps_ph1)
                pr1(s,index) = pr1(s,index) - v_light*dtR(index)';
            end

            pr1_interp(s,index) = lagrange_interp1(time(index), pr1(s,index), time_ref(index), lagr_order);
        else
            bad_sats(s) = 1;
        end
    end
    
    if (any(pr2(s,:)))
        
        index_s = find(pr2(s,:) ~= 0);
        index = intersect(index_e,index_s);
        
        if (length(index) > lagr_order)
            
            if (flag_jumps_ph2)
                pr2(s,index) = pr2(s,index) - v_light*dtR(index)';
            end
            
            pr2_interp(s,index) = lagrange_interp1(time(index), pr2(s,index), time_ref(index), lagr_order);
        else
            bad_sats(s) = 1;
        end
    end
    
    if (any(ph1(s,:)))
        
        index_s = find(ph1(s,:) ~= 0);
        index = intersect(index_e,index_s);
        
        if (length(index) > lagr_order)
            
            if (flag_jumps_ph1)
                ph1(s,index) = ph1(s,index) - v_light*dtR(index)'/lambda(s,1);
                if (flag_doppler_cs && any(dop1(s,index)))
                    dop1(s,index) = dop1(s,index) + v_light*dtRdot(index)'/lambda(s,1);
                end
            end

            [ph1(s,:), cs_found, cs_correction_i] = detect_and_fix_cycle_slips(time, pr1(s,:), ph1(s,:), ph_GF(s,:), dop1(s,:), el(s,:), err_iono(s,:), lambda(s,1));
            
            if ~isempty(cs_correction_i)
                cs_correction_i(:,1) = s;
                cs_correction_i(:,2) = 1;
                status_cs=[status_cs;cs_correction_i];
            end
            
            index_s = find(ph1(s,:) ~= 0);
            index = intersect(index_e,index_s);

            ph1_interp(s,index) = lagrange_interp1(time(index), ph1(s,index), time_ref(index), lagr_order);

%             pr1_interp(s,:) = code_range_to_phase_range(pr1_interp(s,:), ph1_interp(s,:), el(s,:), err_iono(s,:), lambda(s,1));

            if (exist('cs_found', 'var') && cs_found)
                fprintf('Pre-processing: %d cycle-slip(s) detected and fixed on L1 for satellite %02d\n', cs_found, s);
            end
        else
            bad_sats(s) = 1;
        end
        
    elseif (any(pr1(s,:)))
        
        bad_sats(s) = 1;
    end
    
    if (any(ph2(s,:)))
        
        index_s = find(ph2(s,:) ~= 0);
        index = intersect(index_e,index_s);
        
        if (length(index) > lagr_order)
            
            if (flag_jumps_ph2)
                ph2(s,index) = ph2(s,index) - v_light*dtR(index)'/lambda(s,2);
                if (flag_doppler_cs && any(dop2(s,index)))
                    dop2(s,index) = dop2(s,index) + v_light*dtRdot(index)'/lambda(s,2);
                end
            end
            
            [ph2(s,:), cs_found, cs_correction_i] = detect_and_fix_cycle_slips(time, pr2(s,:), ph2(s,:), ph_GF(s,:), dop2(s,:), el(s,:), err_iono(s,:), lambda(s,2));
            
            if ~isempty(cs_correction_i)
                cs_correction_i(:,1) = s;
                cs_correction_i(:,2) = 2;
                status_cs=[status_cs;cs_correction_i];
            end
            
            index_s = find(ph2(s,:) ~= 0);
            index = intersect(index_e,index_s);
            
            ph2_interp(s,index) = lagrange_interp1(time(index), ph2(s,index), time_ref(index), lagr_order);

%             pr2_interp(s,:) = code_range_to_phase_range(pr2_interp(s,:), ph2_interp(s,:), el(s,:), err_iono(s,:), lambda(s,2));
            
            if (exist('cs_found', 'var') && cs_found)
                fprintf('Pre-processing: %d cycle-slip(s) detected and fixed on L2 for satellite %02d\n', cs_found, s);
            end
        else
            bad_sats(s) = 1;
        end
        
    elseif (any(pr2(s,:)))
        
        bad_sats(s) = 1;
    end
end
pr1 = pr1_interp;
pr2 = pr2_interp;
ph1 = ph1_interp;
ph2 = ph2_interp;

end


function [ph, cs_found, cs_correction_i] = detect_and_fix_cycle_slips(time, pr, ph, ph_GF, dop, el, err_iono, lambda)

global cutoff cs_threshold
cs_resolution = 1;
% if cs_threshold==0.5
%     cs_resolution=0.5;
% end

flag_plot = 0;
flag_doppler_cs = 1;

cs_found = 0;
cs_correction_i=[];
cs_correction_count = 0;

cutoff_idx = find(el(1,:) > cutoff);
avail_idx = find(ph(1,:) ~= 0);
idx = intersect(cutoff_idx, avail_idx);

p = polyfit(idx,err_iono(1,idx),3);
err_iono_fit = polyval(p,idx);
% if (flag_plot)
%     figure
%     plot(idx,err_iono(1,idx),'r.')
%     hold on
%     plot(idx,err_iono_fit,'g.')
% end

N_mat = zeros(size(ph));
N_mat(1,idx) = (pr(1,idx) - lambda(1,1).*ph(1,idx) - 2.*err_iono_fit)./lambda(1,1);

if (~isempty(N_mat(1,N_mat(1,:)~=0)))

    delta_thres = 1e30; %cycles
    
    interval = diff(time);
    interval = [interval; interval(end)]';
    
    %initialization
    jmp_code    = (1:(length(ph)-1))';
    jmp_doppler = (1:(length(ph)-1))';
    jmp_deriv   = (1:(length(ph)-1))';
    jmp_GF      = (1:(length(ph)-1))';
    
    %detection (code)
    N_mat(1,N_mat(1,:)==0) = NaN;
%     delta_N = diff(N_mat(1,:));
    delta_code = (diff(N_mat(1,:))./interval(1:end-1))';
    
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
    if (flag_doppler_cs && any(dop) && (sum(~~dop) == sum(~~ph)))
        
%         pred_phase = ph(1:end-1) + ((dop(2:end)+dop(1:end-1))/2).*interval(1:end-1);
        pred_phase = ph(1:end-1) - ((dop(2:end)+dop(1:end-1))/2).*interval(1:end-1);
        
        delta_doppler_all = ((pred_phase - ph(2:end))./interval(1:end-1))';
        delta_doppler_all(ph(2:end)==0) = 0;
        delta_doppler_all(delta_doppler_all == 0) = NaN;
        delta_doppler_all = diff(delta_doppler_all)'./interval(1:end-2);
        pos1 = find(idx~=length(N_mat));
        pos2 = find(idx~=length(N_mat)-1);
        pos = intersect(pos1,pos2);
        delta_doppler(idx(pos)) = delta_doppler_all(idx(pos));
        delta_doppler(delta_doppler == 0) = NaN;
    end
    
    delta_test = delta_doppler;
    not_zero = find(delta_test ~= 0);
    not_nan  = find(~isnan(delta_test));
    avail_doppler = intersect(not_zero, not_nan);

    if (~isempty(delta_test(avail_doppler)))
        outliers = batch_outlier_detection(delta_test(avail_doppler),median(round(interval)));
        [~,jmp_doppler] = intersect(delta_test,outliers);
    end

    idx_interp = setdiff(avail_doppler,jmp_doppler);
    p = polyfit(idx_interp,delta_doppler(idx_interp),1);
    delta_doppler(idx_interp) = delta_doppler(idx_interp) - polyval(p,idx_interp);

    jmp_doppler = sort(jmp_doppler);
    jmp_doppler(diff(jmp_doppler) == 1) = [];
    
    %detection (phase 3rd order derivative)
    delta_deriv = zeros(size(delta_code));
    ph_tmp = ph;
    ph_tmp(1,ph_tmp(1,:)==0) = NaN;
    delta_deriv_all = (diff(diff(diff(ph_tmp(1,:))))./interval(1:end-3).^3)';
    pos1 = find(idx~=length(N_mat));
    pos2 = find(idx~=length(N_mat)-1);
    pos3 = find(idx~=length(N_mat)-2);
    pos = intersect(intersect(pos1,pos2),pos3);
    delta_deriv(idx(pos)) = delta_deriv_all(idx(pos));
    delta_deriv(delta_deriv == 0) = NaN;
    
    delta_test = delta_deriv;
    not_zero = find(delta_test ~= 0);
    not_nan  = find(~isnan(delta_test));
    avail_deriv = intersect(not_zero, not_nan);

    if (~isempty(delta_test(avail_deriv)))
        outliers = batch_outlier_detection(delta_test(avail_deriv),median(round(interval)));
        [~,jmp_deriv] = intersect(delta_test,outliers);
    end
    jmp_deriv = sort(jmp_deriv);
    jmp_deriv(diff(jmp_deriv) == 1) = [];

    %detection (geometry free)
    delta_GF = zeros(size(delta_code));
    ph_GF(1,ph_GF(1,:)==0) = NaN;
    delta_GF_all = (diff(ph_GF)./interval(1:end-1))';
    pos = find(idx~=length(N_mat));
    delta_GF(idx(pos)) = delta_GF_all(idx(pos));
    
    delta_test = delta_GF;
    not_zero = find(delta_test ~= 0);
    not_nan  = find(~isnan(delta_test));
    avail_GF = intersect(not_zero, not_nan);

    if (~isempty(delta_test(avail_GF)))
        outliers = batch_outlier_detection(delta_test(avail_GF),median(round(interval)));
        outliers(outliers < 0.1) = [];
        [~,jmp_GF] = intersect(delta_test,outliers);
    end
    
    %select two observables with low standard deviation
    [min_std_code]    = detect_minimum_std(delta_code(avail_code));
    [min_std_doppler] = detect_minimum_std(delta_doppler(avail_doppler));
    [min_std_deriv]   = detect_minimum_std(delta_deriv(avail_deriv));
    [min_std_GF]      = detect_minimum_std(delta_GF(avail_GF));
    
    min_stds = [min_std_code min_std_doppler min_std_deriv min_std_GF];
    jmps = {jmp_code; jmp_doppler; jmp_deriv; jmp_GF};
    
    [~, pos1] = min(min_stds); min_stds(pos1) = 1e30;
    [~, pos2] = min(min_stds);
    
    jmp = sort(intersect(jmps{pos1},jmps{pos2}));
    
    if (any(delta_doppler))
        delta = -delta_doppler;
    else
        delta = -delta_deriv;
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

    %ignore cycle slips smaller than cs_threshold
    jmp(abs(delta(jmp)) < cs_threshold) = [];

%     jmp1 = intersect(jmp_deriv,jmp_doppler);
%     jmp2 = intersect(jmp1, jmp_code);
%     jmp3 = intersect(jmp2, jmp_GF);
%     jmp = sort(jmp3);

    if (isempty(jmp))
        return
    else
        cs_found = length(jmp);
    end

    if (flag_plot)
        figure; hold on
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
        
        figure
        idx_plot = find(N_mat~=0);
        plot(idx_plot, N_mat(idx_plot));
        
        coltab = colorcube(2*length(jmp));
        c = 1;
    end

%     N_before_zero = 0;
%     N_after_zero = 0;
    
%     if (isempty(jmp))
%         return
%     end

    %fixing
    for j = jmp'
        if (j <= 1 || j >= length(ph)-1)
            check = 0;
            pos_zeros = [];
        else
            pos_zeros = find(N_mat(1,j-1:j+1) == 0,1,'last');
            check = isempty(pos_zeros);
        end
        if (ismember(j,cutoff_idx))% || N_before_zero ~= 0)
            if (check)
                ph_propos = ph(1,j+1) + delta(j);
                if (j == 1)
                    ph_propag = ph_propos;
                else
                    ph_propag = interp1(time(j-1:j),ph(1,j-1:j),time(j+1),'linear','extrap');
                end
                if (abs(ph_propos - ph_propag) < delta_thres)
                    idx = (ph(1,:)==0); 
                    cs_correction = roundmod(delta(j), cs_resolution);
                    cs_correction_count = cs_correction_count + 1;
                    cs_correction_i(cs_correction_count,3)=j;
                    cs_correction_i(cs_correction_count,4)=cs_correction;
                    cs_correction_i(cs_correction_count,5)=delta(j);
                    cs_correction_i(cs_correction_count,6)=0;
%                     if abs(delta(j)-cs_correction) < 0.05
                        cs_correction_i(cs_correction_count,6)=1;
                        N_mat(1,j+1:end) = N_mat(1,j+1:end) - cs_correction;
                        ph   (1,j+1:end) = ph   (1,j+1:end) + cs_correction;
                        N_mat(1,idx) = 0;
                        ph   (1,idx) = 0;
%                     end
                    
                    if (flag_plot)
                        hold on
                        idx_plot = find(N_mat~=0);
                        plot(idx_plot, N_mat(idx_plot),'Color',coltab(2*c-1,:));
                        
                        c = c + 1;
                    end
                end
                
            elseif (pos_zeros == 3)
                
%                 N_before_zero = N_mat(1,j);
                
            elseif (pos_zeros == 2)
                
%                 idx = (ph(1,:)==0);
%                 cs_correction = round(N_mat(1,j+1) - N_before_zero);
%                 N_mat(1,j+1:end) = N_mat(1,j+1:end) - cs_correction;
%                 ph   (1,j+1:end) = ph   (1,j+1:end) + cs_correction;
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
%                 N_after_zero = N_mat(1,j+1);
                
            elseif (pos_zeros == 1)
                
%                 idx = (ph(1,:)==0);
%                 cs_correction = round(N_mat(1,j+1) - N_after_zero);
%                 N_mat(1,j+1:end) = N_mat(1,j+1:end) - cs_correction;
%                 ph   (1,j+1:end) = ph   (1,j+1:end) + cs_correction;
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


function [pr] = code_range_to_phase_range(pr, ph, el, err_iono, lambda)

global cutoff

N_mat = zeros(size(ph));
cutoff_idx = find(el(1,:) > cutoff);
N_mat(1,cutoff_idx) = (pr(1,cutoff_idx) - lambda(1,1).*ph(1,cutoff_idx) - 2.*err_iono(1,cutoff_idx))./lambda(1,1);
if (~isempty(N_mat(1,N_mat(1,:)~=0)))
    N = round(mean(N_mat(1,N_mat(1,:)~=0)));
    idx = find(ph(1,:) ~= 0);
    pr = zeros(size(ph));
    pr(1,idx) = lambda(1,1).*(ph(1,idx)+N);
end
end

function [yi] = lagrange_interp1(x,y,xi,n)
d = n/2;
yi = zeros(size(y));
for t = 1 : length(x)
    if (t<=d)
        yi(t) = LagrangeInter(x(1:n), y(1:n), xi(t));
    elseif (t>(length(x)-d))
        yi(t) = LagrangeInter(x(end-d-1:end), y(end-d-1:end), xi(t));
    else
        yi(t) = LagrangeInter(x(t-d:t+d), y(t-d:t+d), xi(t));
    end
end
end

function [min_std] = detect_minimum_std(time_series)
min_std = 1e30;
d = 5; %half window size
if (length(time_series) < 2*d)
    return
end
mov_std = zeros(size(time_series));
for t = 1 : length(time_series)
    if (t<=d)
        mov_std(t) = std(time_series(1:d));
    elseif (t>(length(time_series)-d))
        mov_std(t) = std(time_series(end-d-1:end));
    else
        mov_std(t) = std(time_series(t-d:t+d));
    end
end
if (~isempty(mov_std))
    min_std = min(mov_std);
end
end

function [outliers] = batch_outlier_detection(time_series, interval)

outlier_thres = 1e-3;
batch_size = 3600; %seconds

num_batches = floor(length(time_series)*interval/batch_size);
batch_idx = 1 : batch_size/interval : batch_size*num_batches/interval;

if (num_batches == 0 && ~isempty(time_series))
    num_batches = 1;
    batch_idx = 1;
end

outliers = [];
for b = 1 : num_batches
    start_idx = batch_idx(b);
    if (b == num_batches)
        end_idx = length(time_series);
    else
        end_idx = batch_idx(b+1)-1;
    end
    [~,~,outliers_batch] = deleteoutliers(time_series(start_idx:end_idx), outlier_thres);
    outliers = [outliers; outliers_batch];
end
end
