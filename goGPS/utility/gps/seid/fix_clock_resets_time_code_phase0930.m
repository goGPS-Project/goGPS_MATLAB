function [ pr1, ph1, pr2, ph2, dtR, dtRdot,time_RdtR ] = fix_clock_resets_time_code_phase0930( pr1, ph1, pr2, ph2, Eph, iono, snr1, time_rx, lambda )
%FIX_CLOCK_RESETS_TIME_CODE_PHASE ‚±‚ÌŠÖ?”‚ÌŠT—v‚ð‚±‚±‚É‹L?q
%   ?Ú?×?à–¾‚ð‚±‚±‚É‹L?q
% SYNTAX:
%   [pr1, ph1, pr2, ph2, dtR, dtRdot] = fix_clock_resets_code_phase(pr1, ph1, pr2, ph2, Eph, iono, snr1, time_rx);
%
% INPUT:
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   Eph = matrix containing 31 ephemerides for each satellite
%   iono = matrix containing ionosphere parameters
%   snr1 = signal-to-noise ratio
%   time_rx = GPS reception time
%
% OUTPUT:
%   pr1 = processed code observation (L1 carrier)
%   ph1 = processed phase observation (L1 carrier)
%   pr2 = processed code observation (L2 carrier)
%   ph2 = processed phase observation (L2 carrier)
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%
% DESCRIPTION:
%   Script that checks and in case fixes observation discontinuities due to
%    receiver clock resets. Estimated receiver clock error and drift for
%    each epoch can also be returned.


global cutoff snr_threshold

v_light = goGNSS.V_LIGHT;

nSatTot = 32;

%detect epochs without data (GPS)
delepochs = find(time_rx == 0);
if ~isempty(delepochs)
    fprintf('%d epochs without data\n', length(delepochs));
end

%number of epochs
nEpochs = length(time_rx);

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK ERROR AND DRIFT
%----------------------------------------------------------------------------------------------

fprintf('Computing receiver clock error and drift...\n');

for i = 1 : nEpochs
    
    Eph_t = rt_find_eph(Eph, time_rx(i), nSatTot);

    sat0 = find(pr1(:,i) ~= 0);

    sat0 = sat0(ismember(sat0, Eph_t(30,:)));
    
    if (size(sat0,1) >= 4)

        [~, dtR_tmp, ~, ~, ~, ~, ~, ~, ~, sat, ~, ~, ~, ~, ~, ~, ~, ~, ~] = init_positioning(time_rx(i), pr1(sat0,i), snr1(sat0,i), Eph_t, [], iono, [], [], [], [], sat0, [], lambda(sat0,:), cutoff, snr_threshold, 1, 0, 0);

%     [~, dtR(i), ~, ~, ~, ~, ~, ~, ~, sat, ~, ~, ~, ~, ~, ~, ~, ~, ~] = init_positioning(time_rx(i), pr1(sat0,i), snr1(sat0,i), Eph_t, [], [], [], iono, [], [], [], [], sat0, cutoff, snr_threshold, 0, 0);
% for reference;
% [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pseudorange, snr, Eph, SP3_time, SP3_coor, SP3_clck, iono, sbas, XR0, XS0, dtS0, sat0, cutoff_el, cutoff_snr, flag_XR, flag_XS)
%     [~, dtR(i), ~, ~, ~, ~, ~, ~, ~, sat] = init_positioning(time_rx(i), pr1(sat0,i), snr1(sat0,i), Eph_t, [], iono, [], [], [], [], sat0, cutoff, snr_threshold, 0, 0);
        if (size(sat,1) >= 4)
            
            dtR(i) = dtR_tmp;
            
            if (i > 1)
                %receiver clock drift
                if (dtR(i) ~= 0 && dtR(i-1) ~= 0)
                    dtRdot(i-1) = (dtR(i) - dtR(i-1))/(time_rx(i) - time_rx(i-1));
                end
            end
        else
            if (i > 2)
                dtR(i) = dtR(i-1) + (dtR(i-1) - dtR(i-2));
                dtRdot(i-1) = (dtR(i) - dtR(i-1));
            end
        end
    end
end

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK DISCONTINUITIES (DETECTION)
%----------------------------------------------------------------------------------------------

fprintf('Detecting clock resets... ');

%check if there is any discontinuity in the clock drift
clock_thresh = 1e-5;
disc = find(abs(dtRdot-mean(dtRdot)) > clock_thresh);

%remove unavailable epochs from discontinuities, because otherwise
%Doppler-predicted observations are going to be used on the wrong epochs;
%if an epoch is missing the discontinuity cannot be fixed
for i = 1 : length(delepochs)
	disc = setdiff(disc,(delepochs(i)-i));
end

% %display the clock error
% figure
% plot(dtR,'.')
% %display the absolute value of the clock drift
% figure
% plot(abs(dtRdot),'.')
% hold on
% plot([1, nEpochs], [clock_thresh clock_thresh],'r');

%----------------------------------------------------------------------------------------------
% OBSERVATION CORRECTION
%----------------------------------------------------------------------------------------------

%correct pseudorange for receiver clock error
% if (~isempty(disc))
    
    fprintf('%d detected.\n', length(disc));
    
    fprintf('Fixing %d clock resets...\n', length(disc));
    %index_e = find(time_rx ~= 0);
    time_RdtR = time_rx - dtR;
    for s = 1 : 32
%         for i = 1 : nEpochs
            if (any(ph1(s,:)))
                index = find(pr1(s,:)~=0);
                %index= intersect(index_e,index_s);
                pr1(s,index) = pr1(s,index) - v_light*dtR(index)';
                ph1(s,index) = ph1(s,index) - v_light*dtR(index)'/lambda(s,1);
               % figure
% plot(time_RdtR(index), pr1(s,index), '-*')
% hold on
                pr1(s,index) = interp1(time_RdtR(index), pr1(s,index), time_rx(index), 'spline','extrap');
                ph1(s,index) = interp1(time_RdtR(index), ph1(s,index), time_rx(index), 'spline','extrap');
% plot(time_rx(index), pr1(s,index), 'og')
% pause
            end
            if (any(ph2(s,:)))
                index = find(pr2(s,:)~=0);
                %index= intersect(index_e,index_s);
                pr2(s,index) = pr2(s,index) - v_light*dtR(index)';
                ph2(s,index) = ph2(s,index) - v_light*dtR(index)'/lambda(s,2);
                
                pr2(s,index) = interp1(time_RdtR(index), pr2(s,index), time_rx(index), 'spline','extrap');
                ph2(s,index) = interp1(time_RdtR(index), ph2(s,index), time_rx(index), 'spline','extrap');
            end
%         end
    end
    
    
     
% else
%     fprintf('none detected.\n');
% end

end

