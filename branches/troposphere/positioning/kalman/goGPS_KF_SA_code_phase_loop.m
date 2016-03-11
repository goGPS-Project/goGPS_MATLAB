function [check_on, check_off, check_pivot, check_cs] = goGPS_KF_SA_code_phase_loop(time_rx, pr1, ph1, dop1, pr2, ph2, dop2, snr, Eph, SP3, iono, sbas, lambda, frequencies, obs_comb, flag_IAR, flag_tropo, antenna_PCV) %#ok<INUSL>

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = goGPS_KF_SA_code_phase_loop(time_rx, pr1, ph1, dop1, pr2, ph2, dop2, snr, Eph, SP3, iono, sbas, lambda, frequencies, obs_comb, flag_IAR, flag_tropo, antenna_PCV);
%
% INPUT:
%   time_rx = GPS time
%   pr1  = ROVER-SATELLITE code pseudorange (L1 carrier)
%   ph1  = ROVER-SATELLITE phase observation (L1 carrier)
%   dop1 = ROVER_SATELLITE Doppler observation (L1 carrier)
%   pr2  = ROVER-SATELLITE code pseudorange (L2 carrier)
%   ph2  = ROVER-SATELLITE phase observation (L2 carrier)
%   dop2 = ROVER_SATELLITE Doppler observation (L2 carrier)
%   snr = signal-to-noise ratio for ROVER observations
%   Eph = satellite ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono =  ionospheric parameters (vector of zeroes if not available)
%   sbas = SBAS corrections
%   lambda = wavelength matrix (depending on the enabled constellations)
%   frequencies = L1 carrier (frequencies=1), L2 carrier (frequencies=2)
%   obs_comb = observations combination (e.g. iono-free: obs_comb = 'IONO_FREE')
%   flag_IAR = boolean variable to enable/disable integer ambiguity resolution
%   flag_tropo = boolean variable to enable/disable tropospheric delay estimation
%   antenna_PCV = antenna phase center variation
%
% OUTPUT:
%   check_on = boolean variable for satellite addition
%   check_off = boolean variable for satellite loss
%   check_pivot = boolean variable for pivot change
%   check_cs = boolean variable for cycle-slip
%
% DESCRIPTION:
%   Kalman filter for the ROVER trajectory computation.
%   Standalone positioning using code and phase.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Andrea Nardo, Stefano Caldera
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

global sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_tropo sigmaq_rclock sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_codIF sigmaq_ph sigmaq_phIF sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold o1 o2 o3 nN nT nC
global tile_header tile_georef dtm_dir
global h_antenna zero_time

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global doppler_pred_range1_R doppler_pred_range2_R
global ratiotest mutest succ_rate fixed_solution

global t residuals_fixed residuals_float outliers s02_ls

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables to point out events (satellite addition, losses, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1,1);

%azimuth, elevation and ROVER-satellite and MASTER-satellite distances
azR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

%iono-free coefficients
alpha1 = (goGNSS.F1^2/(goGNSS.F1^2 - goGNSS.F2^2));
alpha2 = (goGNSS.F2^2/(goGNSS.F1^2 - goGNSS.F2^2));

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
% (if a static model is used, no noise is added)
Cvv = zeros(o3+nN+nT+nC);
if (o1 > 1)
    Cvv(o1,o1) = sigmaq_vE;
    Cvv(o2,o2) = sigmaq_vN;
    Cvv(o3,o3) = sigmaq_vU;
    
    %propagate diagonal local cov matrix to global cov matrix
    Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), X_t1_t([1 o1+1 o2+1]));
end
if (flag_tropo)
    Cvv(o3+nN+1:o3+nN+nT,o3+nN+1:o3+nN+nT) = sigmaq_tropo * eye(nT);
end
Cvv(o3+nN+nT+1:o3+nN+nT+nC,o3+nN+nT+1:o3+nN+nT+nC) = sigmaq_rclock * eye(nC);

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
if (length(frequencies) == 2)
    sat_pr = find( (pr1(:,1) ~= 0) & (pr2(:,1) ~= 0) );
    sat = find( (pr1(:,1) ~= 0) & (ph1(:,1) ~= 0) & ...
                (pr2(:,1) ~= 0) & (ph2(:,1) ~= 0) );
else
    if (frequencies == 1)
        sat_pr = find(pr1(:,1) ~= 0);
        sat = find( (pr1(:,1) ~= 0) & (ph1(:,1) ~= 0) );
    else
        sat_pr = find(pr2(:,1) ~= 0);
        sat = find( (pr2(:,1) ~= 0) & (ph2(:,1) ~= 0) );
    end
end
if (isempty(SP3))
    eph_avail = Eph(30,:);
else
    eph_avail = SP3.avail;
end
sat_pr = sat_pr(ismember(sat_pr, eph_avail));
sat = sat(ismember(sat, eph_avail));

%only satellites with code and phase
%sat_pr = sat;

%previous satellite configuration (with phase measurements)
sat_old = find(conf_sat == 1);

%number of visible satellites
nsat = size(sat_pr,1);

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1(sat_pr) = pr1(sat_pr) + sbas.prc(sat_pr)';
end

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

%approximate position
XR0 = X_t1_t([1,o1+1,o2+1]);
flag_XR = 2;

%approximated coordinates X Y Z
X_app = XR0(1);
Y_app = XR0(2);
Z_app = XR0(3);

%----------------------------------------------------------------------------------------
% CONVERSION FROM CARTESIAN TO GEODETIC COORDINATES
%----------------------------------------------------------------------------------------
[phiR_app, lamR_app, hR_app] = cart2geod(X_app, Y_app, Z_app);

%----------------------------------------------------------------------------------------
% EXTRACTION OF THE HEIGHT PSEUDO-OBSERVATION FROM THE DTM
%----------------------------------------------------------------------------------------

%conversion from cartesian to geodetic coordinates
[phi_app, lam_app, h_app] = cart2geod(X_app, Y_app, Z_app);

%projection to UTM coordinates
[E_app, N_app] = geod2plan(phi_app, lam_app);

%dtm tile detection (in which the approximated position lies)
[tile_row,tile_col] = find ( (E_app > tile_georef(:,:,1)) & (E_app <= tile_georef(:,:,4)) & (N_app >= tile_georef(:,:,3)) & (N_app < tile_georef(:,:,2)));

%tile buffer dimension
tile_buffer_size = 3;

%check if the approximated position lies within one of the available tiles, otherwise set nodata value
if ( ~isempty(tile_row) && ~isempty(tile_col) )
    tile_buffer = cell(tile_buffer_size,tile_buffer_size);
    for i = -1 : 1
        for j = -1 : 1
            %definition of the path and the filename of the selected tile
            tile_path = strcat(dtm_dir,'/tiles/tile_',num2str(tile_row+i),'_',num2str(tile_col+j),'.mat');

            %check the existence of the file associated to the selected tile
            fid = fopen(tile_path,'r');
            if (fid ~= -1)
                fclose(fid);
                %load the selected tile
                load(tile_path, 'tile');
            else
                %load of a null tile
                tile(1:tile_header.nrows, 1:tile_header.ncols) = tile_header.nodata;
            end
            %buffer creation around the selected tile
            tile_buffer{i+2,j+2} =  tile;
        end
    end

    %buffer conversion from cell to matrix
    tile_buffer = cell2mat(tile_buffer);

    %computation of the tile buffer dimension (cell number)
    [tile_height, tile_width] = size(tile_buffer);

    %tile buffer lower left center coordinates extraction
    Ell = tile_georef(tile_row,tile_col,1) - tile_width/tile_buffer_size*tile_header.cellsize + tile_header.cellsize/2;
    Nll = tile_georef(tile_row,tile_col,3) - tile_height/tile_buffer_size*tile_header.cellsize + tile_header.cellsize/2;

    %extraction from the dtm of the height correspondent to the approximated position
    [h_dtm] = grid_bilin_interp(E_app, N_app, tile_buffer, tile_header.ncols*3, tile_header.nrows*3, tile_header.cellsize, Ell, Nll, tile_header.nodata);

    %antenna height addition
    h_dtm = h_dtm + h_antenna;
else
    h_dtm = tile_header.nodata;
end

%------------------------------------------------------------------------------------
% OBSERVATION EQUATIONS
%------------------------------------------------------------------------------------

%if the number of visible satellites is equal or greater than min_nsat
if (nsat >= min_nsat)
    
    %----------------------------------------------------------------------------------------
    % SATELLITE POSITION AND RECEIVER CLOCK ERROR COMPUTATION
    %----------------------------------------------------------------------------------------
    
    sat_pr_old = sat_pr;
    
    if (frequencies(1) == 1)
        if (length(frequencies) < 2 || ~strcmp(obs_comb,'IONO_FREE'))
            [~, ~, XS, dtS, ~, ~, ~, err_tropo, err_iono1, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys] = init_positioning(time_rx, pr1(sat_pr), snr(sat_pr), Eph, SP3, iono, sbas, XR0, [], [], sat_pr, [], lambda(sat_pr,:), cutoff, snr_threshold, frequencies, flag_XR, 0);
        else
            [~, ~, XS, dtS, ~, ~, ~, err_tropo, err_iono1, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys] = init_positioning(time_rx, alpha1*pr1(sat_pr) - alpha2*pr2(sat_pr), snr(sat_pr), Eph, SP3, zeros(8,1), sbas, XR0, [], [], sat_pr, [], zeros(length(sat_pr),2), cutoff, snr_threshold, frequencies, flag_XR, 0);
        end
        
        if (~isempty(sat_pr))
            err_iono2 = err_iono1 .* ionoFactor(sat_pr,2);
        end
    else
        [~, ~, XS, dtS, ~, ~, ~, err_tropo, err_iono2, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys] = init_positioning(time_rx, pr2(sat_pr), snr(sat_pr), Eph, SP3, iono, sbas, XR0, [], [], sat_pr, [], lambda(sat_pr,:), cutoff, snr_threshold, frequencies, flag_XR, 0);
        
        if (~isempty(sat_pr))
            err_iono1 = err_iono2 ./ ionoFactor(sat_pr,2);
        end
    end
    
    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_old, sat_pr);
    sat(ismember(sat,sat_removed)) = [];
    
    %index of satellite with phase among those with code
    [~, phase_index] = intersect(sat_pr,sat);
    
    %----------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING
    %----------------------------------------------------------------------------------------
    
    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr) = -1;
    conf_sat(sat) = +1;
    
    %cycle-slip configuration
    conf_cs = zeros(nSatTot,1);
    
    %number of visible satellites
    nsat = size(sat_pr,1);
    n = nsat;
    
    %previous pivot (not used)
    if (pivot ~= 0)
        pivot_old = pivot;
    end
    
    %current pivot (not used)
    if ~isempty(sat)
        [max_elR, i] = max(elR(sat)); %#ok<ASGLU>
        pivot = sat(i);
    else
        [max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
        pivot = sat_pr(i);
    end
    %pivot = find(elR == max(elR));
    
    %if the number of available satellites after the cutoffs is equal or greater than min_nsat
    if (nsat >= min_nsat)
        
        %------------------------------------------------------------------------------------
        % SATELLITE ADDITION/LOSS
        %------------------------------------------------------------------------------------

        %search for a lost satellite
        sat_dead = setdiff(sat_old,sat);
        
        if (~isempty(sat_dead))
            
            check_off = 1;
            
            %for lost satellites it is fundamental to set their ambiguities to 0.
            %Furthermore it could be convenient to raise their uncertainty
            %(not necessary - done when a new satellite is added)
            N1 = 0;
            N2 = 0;
            
            if (length(frequencies) == 2)
                if (strcmp(obs_comb,'NONE'))
                    X_t1_t(o3+sat_dead,1) = N1;
                    X_t1_t(o3+nSatTot+sat_dead,1) = N2;
                elseif (strcmp(obs_comb,'IONO_FREE'))
                    X_t1_t(o3+sat_dead,1) = N1;
                end
            else
                if (frequencies == 1)
                    X_t1_t(o3+sat_dead,1) = N1;
                else
                    X_t1_t(o3+sat_dead,1) = N2;
                end
            end
        end
        
        %search for a new satellite
        sat_born = setdiff(sat,sat_old);
        
        if (~isempty(sat_born))
            
            check_on = 1;
        end
        
        %------------------------------------------------------------------------------------
        % PIVOT CHANGE (NOT REQUIRED FOR STANDALONE)
        %------------------------------------------------------------------------------------
        
        %search for a possible PIVOT change
        if (pivot ~= pivot_old && pivot_old ~= 0)
            
            check_pivot = 1;
        end
        
        %------------------------------------------------------------------------------------
        % CYCLE-SLIP
        %------------------------------------------------------------------------------------
        
        if ~isempty(sat)
            
            %Test presence/absence of a cycle-slip at the current epoch.
            %The state of the system is not changed yet
            if (length(frequencies) == 2)
                if (strcmp(obs_comb,'NONE'))
                    [check_cs1, N_slip1, sat_slip1] = cycle_slip_detection_SA(X_t1_t(o3+1:o3+nSatTot),           pr1(sat), ph1(sat), err_iono1(phase_index), doppler_pred_range1_R(sat), sat, sat_born, cs_threshold, lambda(sat,1)); %#ok<ASGLU>
                    [check_cs2, N_slip2, sat_slip2] = cycle_slip_detection_SA(X_t1_t(o3+nSatTot+1:o3+nSatTot*2), pr2(sat), ph2(sat), err_iono2(phase_index), doppler_pred_range2_R(sat), sat, sat_born, cs_threshold, lambda(sat,2)); %#ok<ASGLU>
                    
                    if (check_cs1 || check_cs2)
                        check_cs = 1;
                    end
                elseif (strcmp(obs_comb,'IONO_FREE'))
                    [check_cs, N_slip, sat_slip] = cycle_slip_detection_SA(X_t1_t(o3+1:o3+nSatTot), alpha1*pr1(sat) - alpha2*pr2(sat), alpha1*lambda(sat,1).*ph1(sat) - alpha2*lambda(sat,2).*ph2(sat), zeros(size(sat)), alpha1*doppler_pred_range1_R(sat) - alpha2*doppler_pred_range2_R(sat), sat, sat_born, cs_threshold, ones(size(sat_pr))); %#ok<ASGLU>
                end
            else
                if (frequencies == 1)
                    [check_cs, N_slip, sat_slip] = cycle_slip_detection_SA(X_t1_t(o3+1:o3+nSatTot), pr1(sat), ph1(sat), err_iono1(phase_index), doppler_pred_range1_R(sat), sat, sat_born, cs_threshold, lambda(sat,1)); %#ok<ASGLU>
                else
                    [check_cs, N_slip, sat_slip] = cycle_slip_detection_SA(X_t1_t(o3+1:o3+nSatTot), pr2(sat), ph2(sat), err_iono2(phase_index), doppler_pred_range2_R(sat), sat, sat_born, cs_threshold, lambda(sat,2)); %#ok<ASGLU>
                end
            end
        else
            sat_slip1 = [];
            sat_slip2 = [];
            sat_slip = [];
            check_cs1 = 0;
            check_cs2 = 0;
            check_cs = 0;
        end
        
        %------------------------------------------------------------------------------------
        % PHASE AMBIGUITY ESTIMATION
        %------------------------------------------------------------------------------------
        
        if (length(frequencies) == 2)
            if (strcmp(obs_comb,'NONE'))
                
                [N1_slip, N1_born] = ambiguity_init_SA(XR0, XS, dtS, pr1(sat_pr), ph1(sat_pr), snr(sat_pr), elR(sat_pr), sat_pr, sat, sat_slip1, sat_born, distR(sat_pr), err_tropo, err_iono1, sys, lambda(sat_pr,1), X_t1_t(o3+        sat), Cee(o3+        sat, o3+        sat), X_t1_t(o3+nN+(1:nT)), Cee(o3+nN+(1:nT), o3+nN+(1:nT)), X_t1_t(o3+nN+nT+(1:nC)), Cee(o3+nN+nT+(1:nC), o3+nN+nT+(1:nC)));
                [N2_slip, N2_born] = ambiguity_init_SA(XR0, XS, dtS, pr2(sat_pr), ph2(sat_pr), snr(sat_pr), elR(sat_pr), sat_pr, sat, sat_slip2, sat_born, distR(sat_pr), err_tropo, err_iono2, sys, lambda(sat_pr,2), X_t1_t(o3+nSatTot+sat), Cee(o3+nSatTot+sat, o3+nSatTot+sat), X_t1_t(o3+nN+(1:nT)), Cee(o3+nN+(1:nT), o3+nN+(1:nT)), X_t1_t(o3+nN+nT+(1:nC)), Cee(o3+nN+nT+(1:nC), o3+nN+nT+(1:nC)));
                
                if (check_on)
                    X_t1_t(o3+sat_born,1) = N1_born;
                    X_t1_t(o3+nSatTot+sat_born,1) = N2_born;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1_born * eye(size(sat_born,1));
                    %Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = sigmaq_N2_born * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                end
                
                if (check_cs1)
                    conf_cs(sat_slip1) = 1;
                    X_t1_t(o3+sat_slip1) = N1_slip;
                    Cvv(o3+sat_slip1,o3+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
                end
                
                if (check_cs2)
                    conf_cs(sat_slip2) = 1;
                    X_t1_t(o3+nSatTot+sat_slip2) = N2_slip;
                    Cvv(o3+nSatTot+sat_slip2,o3+nSatTot+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
                end
            elseif (strcmp(obs_comb,'IONO_FREE'))
                
                [N_slip, N_born] = ambiguity_init_SA(XR0, XS, dtS, alpha1*pr1(sat_pr) - alpha2*pr2(sat_pr), alpha1*lambda(sat_pr,1).*ph1(sat_pr) - alpha2*lambda(sat_pr,2).*ph2(sat_pr), snr(sat_pr), elR(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), err_tropo, zeros(size(sat_pr)), sys, ones(size(sat_pr)), X_t1_t(o3+sat), Cee(o3+sat, o3+sat), X_t1_t(o3+nN+(1:nT)), Cee(o3+nN+(1:nT), o3+nN+(1:nT)), X_t1_t(o3+nN+nT+(1:nC)), Cee(o3+nN+nT+(1:nC), o3+nN+nT+(1:nC)));
                
                if (check_on)
                    X_t1_t(o3+sat_born,1) = N_born;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                end
                
                if (check_cs)
                    conf_cs(sat_slip) = 1;
                    X_t1_t(o3+sat_slip) = N_slip;
                    Cvv(o3+sat_slip,o3+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
                end
            end
        else
            if (frequencies == 1)
                [N_slip, N_born] = ambiguity_init_SA(XR0, XS, dtS, pr1(sat_pr), ph1(sat_pr), snr(sat_pr), elR(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), err_tropo, err_iono1, sys, lambda(sat_pr,1), X_t1_t(o3+sat), Cee(o3+sat, o3+sat), X_t1_t(o3+nN+(1:nT)), Cee(o3+nN+(1:nT), o3+nN+(1:nT)), X_t1_t(o3+nN+nT+(1:nC)), Cee(o3+nN+nT+(1:nC), o3+nN+nT+(1:nC)));
            else
                [N_slip, N_born] = ambiguity_init_SA(XR0, XS, dtS, pr2(sat_pr), ph2(sat_pr), snr(sat_pr), elR(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), err_tropo, err_iono2, sys, lambda(sat_pr,2), X_t1_t(o3+sat), Cee(o3+sat, o3+sat), X_t1_t(o3+nN+(1:nT)), Cee(o3+nN+(1:nT), o3+nN+(1:nT)), X_t1_t(o3+nN+nT+(1:nC)), Cee(o3+nN+nT+(1:nC), o3+nN+nT+(1:nC)));
            end
            
            if (check_on)
                X_t1_t(o3+sat_born,1) = N_born;
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
                Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
            
            if (check_cs)
                conf_cs(sat_slip) = 1;
                X_t1_t(o3+sat_slip) = N_slip;
                Cvv(o3+sat_slip,o3+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
            end
            
        end
        
        %------------------------------------------------------------------------------------
        % OBSERVATION EQUATIONS
        %------------------------------------------------------------------------------------
        
        %rows in which the phase observation is available
        p = find(ismember(sat_pr,sat)==1);
        
        %compute PCV: phase and code 1 
        [~, index_ph]=intersect(sat_pr,sat);
        
        if (~isempty(antenna_PCV) && antenna_PCV(1).n_frequency ~= 0) % rover
            index_rover=1; 
            PCV1=PCV_interp(antenna_PCV(index_rover), 90-elR(sat_pr), azR(sat_pr), sys, 1);
            pr1(sat_pr)=pr1(sat_pr)-PCV1;
            ph1(sat)=ph1(sat)-PCV1(index_ph)./lambda(sat,1);
            
            if (length(frequencies) == 2 || frequencies(1) == 2)
                PCV2=PCV_interp(antenna_PCV(index_rover), 90-elR(sat_pr), azR(sat_pr), sys, 2);
                pr2(sat_pr)=pr2(sat_pr)-PCV2;
                ph2(sat)=ph2(sat)-PCV2(index_ph)./lambda(sat,2);
            end
        end
        
        %when the tropospheric delay is estimated, only its hydrostatic part is removed from the observations
        if (flag_tropo)
            
            [week, sow] = time2weektow(time_rx + zero_time);
            date = gps2date(week, sow);
            [~, mjd] = date2jd(date);
            
            %pressure = goGNSS.STD_PRES;
            %temperature = goGNSS.STD_TEMP;
            %humidity = goGNSS.STD_HUMI;
            
            [pres_R, temp_R, undu_R] = gpt(mjd, phiR_app, lamR_app, hR_app);
            ZHD_R = saast_dry(pres_R, hR_app - undu_R, phiR_app*180/pi);
%             ZWD_R = saast_wet(temp_R, goGNSS.STD_HUMI, hR_app - undu_R);

%             ZHD_R = 2.3 * exp(-0.116e-3 * (hR_app - undu_R));
            ZWD_R = 0.1;
            
            gmfh_R = zeros(size(err_tropo));
            gmfw_R = zeros(size(err_tropo));
            err_tropo0 = zeros(size(err_tropo));
            beta_R = zeros(n,nT);
            for s = 1 : n
                [gmfh_R(s,1), gmfw_R(s,1)] = gmf_f_hu(mjd, phiR_app, lamR_app, hR_app, (90-elR(sat_pr(s),1))*pi/180);
                err_tropo0(s,1) = gmfh_R(s,1)*ZHD_R + gmfw_R(s,1)*ZWD_R;
            end
            
            beta_R(:,nT) = gmfw_R(:,ones(nT,1));
        else
            err_tropo0 = err_tropo;
            beta_R = zeros(n,nT);
        end
        
        %function that calculates the Kalman filter parameters
        [alpha, prapp_pr1, prapp_ph1, prapp_pr2, prapp_ph2, probs_prIF, probs_phIF, prapp_prIF, prapp_phIF] = input_kalman_SA(XR0, XS, pr1(sat_pr), ph1(sat_pr), pr2(sat_pr), ph2(sat_pr), distR(sat_pr), dtS, err_tropo0, err_iono1, err_iono2, lambda(sat_pr,:));
        
        %zeroes vector useful in matrix definitions
        Z_1_nN = zeros(1,nN);
        Z_n_nN = zeros(n,nN);
        Z_n_om = zeros(n,o1-1);
        Z_1_om = zeros(1,o1-1);
        Z_1_nT = zeros(1,nT);
        Z_1_nC = zeros(1,nC);
        gamma  = ones(n,nC);
        
        %H matrix computation for the code
        H_cod1  = [alpha(:,1) Z_n_om alpha(:,2) Z_n_om alpha(:,3) Z_n_om Z_n_nN beta_R gamma];
        H_cod2  = [alpha(:,1) Z_n_om alpha(:,2) Z_n_om alpha(:,3) Z_n_om Z_n_nN beta_R gamma];
        H_codIF = [alpha(:,1) Z_n_om alpha(:,2) Z_n_om alpha(:,3) Z_n_om Z_n_nN beta_R gamma];
        if (length(frequencies) == 2)
            if (strcmp(obs_comb,'NONE'))
                H_cod = [H_cod1; H_cod2];
            elseif (strcmp(obs_comb,'IONO_FREE'))
                H_cod = H_codIF;
            end
        else
            if (frequencies == 1)
                H_cod = H_cod1;
            else
                H_cod = H_cod2;
            end
        end
        
        %lambda positions computation
        L_pha1  = zeros(n,nSatTot);
        L_pha2  = zeros(n,nSatTot);
        L_phaIF = zeros(n,nSatTot);
        for u = 1 : n
            L_pha1(u,sat_pr(u))  = -(lambda(sat_pr(u),1));
            L_pha2(u,sat_pr(u))  = -(lambda(sat_pr(u),2));
            L_phaIF(u,sat_pr(u)) = -1;
        end
        
        %H matrix computation for the phase
        if ~isempty(p)
            H_pha1  = [alpha(p,1) Z_n_om(p,:) alpha(p,2) Z_n_om(p,:) alpha(p,3) Z_n_om(p,:) Z_n_nN(p,:) beta_R(p,:) gamma(p,:)];
            H_pha2  = [alpha(p,1) Z_n_om(p,:) alpha(p,2) Z_n_om(p,:) alpha(p,3) Z_n_om(p,:) Z_n_nN(p,:) beta_R(p,:) gamma(p,:)];
            H_phaIF = [alpha(p,1) Z_n_om(p,:) alpha(p,2) Z_n_om(p,:) alpha(p,3) Z_n_om(p,:) Z_n_nN(p,:) beta_R(p,:) gamma(p,:)];
            if (length(frequencies) == 2)
                if (strcmp(obs_comb,'NONE'))
                    H_pha1(:,o3+1:o3+nSatTot) = L_pha1(p,:);
                    H_pha2(:,o3+nSatTot+1:o3+nSatTot*2) = L_pha2(p,:);
                    H_pha = [H_pha1; H_pha2];
                elseif (strcmp(obs_comb,'IONO_FREE'))
                    H_phaIF(:,o3+1:o3+nSatTot) = L_phaIF(p,:);
                    H_pha = H_phaIF;
                end
            else
                if (frequencies == 1)
                    H_pha1(:,o3+1:o3+nSatTot) = L_pha1(p,:);
                    H_pha = H_pha1;
                else
                    H_pha2(:,o3+1:o3+nSatTot) = L_pha2(p,:);
                    H_pha = H_pha2;
                end
            end
        else
            H_pha = [];
        end
        
        %H matrix computation for the DTM pseudo-observation
        H_dtm = [];
        if (h_dtm ~= tile_header.nodata)
            H_dtm = [cos(phi_app)*cos(lam_app) Z_1_om cos(phi_app)*sin(lam_app) Z_1_om sin(phi_app) Z_1_om Z_1_nN Z_1_nT Z_1_nC];
        end
        
        %construction of the complete H matrix
        H = [H_cod; H_pha; H_dtm];
        
        %Y0 vector computation for the code
        y0_cod1  = pr1(sat_pr) - prapp_pr1  + alpha(:,1)*X_app + alpha(:,2)*Y_app + alpha(:,3)*Z_app;
        y0_cod2  = pr2(sat_pr) - prapp_pr2  + alpha(:,1)*X_app + alpha(:,2)*Y_app + alpha(:,3)*Z_app;
        y0_codIF = probs_prIF  - prapp_prIF + alpha(:,1)*X_app + alpha(:,2)*Y_app + alpha(:,3)*Z_app;
        
        %Y0 vector computation for the phase
        if ~isempty(p)
            y0_pha1  = ph1(sat).*lambda(sat,1) - prapp_ph1(p)  + alpha(p,1)*X_app + alpha(p,2)*Y_app + alpha(p,3)*Z_app;
            y0_pha2  = ph2(sat).*lambda(sat,2) - prapp_ph2(p)  + alpha(p,1)*X_app + alpha(p,2)*Y_app + alpha(p,3)*Z_app;
            y0_phaIF = probs_phIF(p)           - prapp_phIF(p) + alpha(p,1)*X_app + alpha(p,2)*Y_app + alpha(p,3)*Z_app;
        else
            y0_pha1 = [];
            y0_pha2 = [];
            y0_phaIF = [];
        end
        
        %Y0 vector computation for DTM constrain
        y0_dtm = [];
        if (h_dtm ~= tile_header.nodata)
            y0_dtm = h_dtm - h_app + cos(phi_app)*cos(lam_app)*X_app + cos(phi_app)*sin(lam_app)*Y_app + sin(phi_app)*Z_app;
        end
        
        %construction of the total Y0 vector
        if (length(frequencies) == 2)
            if (strcmp(obs_comb,'NONE'))
                y0_cod = [y0_cod1; y0_cod2];
                y0_pha = [y0_pha1; y0_pha2];
            elseif (strcmp(obs_comb,'IONO_FREE'))
                y0_cod = y0_codIF;
                y0_pha = y0_phaIF;
            end
        else
            if (frequencies == 1)
                y0_cod = y0_cod1;
                y0_pha = y0_pha1;
            else
                y0_cod = y0_cod2;
                y0_pha = y0_pha2;
            end
        end
        y0 = [y0_cod; y0_pha; y0_dtm];
        
        %------------------------------------------------------------------------------------
        % OBSERVATION COVARIANCE MATRIX
        %------------------------------------------------------------------------------------
        
        %construction of the cofactor matrix
        Q = cofactor_matrix_SA(elR(sat_pr), snr(sat_pr));
        
        %zeroes vector useful in matrix definitions
        Z_n_n = zeros(n,n);
        
        %multiplication by the code variance and the phase variance to build the matrix
        if ~isempty(p)
            if (length(frequencies) == 2)
                if (strcmp(obs_comb,'NONE'))
                    Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,:) Z_n_n(:,p) Z_n_n(:,p); Z_n_n(:,:) sigmaq_cod2*Q(:,:) Z_n_n(:,p) Z_n_n(:,p);
                        Z_n_n(p,:) Z_n_n(p,:) sigmaq_ph*Q(p,p) Z_n_n(p,p); Z_n_n(p,:) Z_n_n(p,:) Z_n_n(p,p) sigmaq_ph*Q(p,p)];
                elseif (strcmp(obs_comb,'IONO_FREE'))
                    Cnn = [sigmaq_codIF*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_phIF*Q(p,p)];
                end
            else
                if (frequencies == 1)
                    Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
                else
                    Cnn = [sigmaq_cod2*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
                end
            end
        else
            if (length(frequencies) == 2)
                Cnn = [sigmaq_cod1*Q Z_n_n; Z_n_n sigmaq_cod2*Q];
            else
                if (frequencies == 1)
                    Cnn = sigmaq_cod1*Q;
                else
                    Cnn = sigmaq_cod2*Q;
                end
            end
        end
        if (h_dtm ~= tile_header.nodata)
            Cnn(end+1,end+1) = sigmaq_dtm;
        end
        
        %------------------------------------------------------------------------------------
        % OUTLIER DETECTION (OPTIMIZED LEAVE ONE OUT)
        %------------------------------------------------------------------------------------
        
        search_for_outlier = 1;

        if (length(frequencies) == 2)
            if (strcmp(obs_comb,'NONE'))
                index_residuals_outlier=[sat_pr;nSatTot+sat_pr;nSatTot*2+sat;nSatTot*3+sat];  %[code;phase]
            elseif (strcmp(obs_comb,'IONO_FREE'))
                index_residuals_outlier=[sat_pr;nSatTot*2+sat];  %[code;phase]
            end
        else
            if (frequencies == 1)
                index_residuals_outlier=[sat_pr;nSatTot*2+sat];  %[code;phase]
            else
                index_residuals_outlier=[sat_pr;nSatTot*2+sat];  %[code;phase]
            end
        end

        if (h_dtm ~= tile_header.nodata)
            y0_residuals=y0(1:end-1);
            H1_residuals=H(1:end-1,:);
            y0_noamb=y0(1:end-1);
            H1=H(1:end-1,[1 o1+1 o2+1]);
        else
            y0_residuals=y0;
            H1_residuals=H;
            y0_noamb=y0;
            H1=H(:,[1 o1+1 o2+1]);
        end

        if (~isempty(sat))
            if (length(frequencies) == 2)
                if (strcmp(obs_comb,'NONE'))
                    y0_noamb(length(sat_pr)*2+            (1:length(sat)))=y0_noamb(length(sat_pr)*2+            (1:length(sat)))+lambda(sat,1).*X_t1_t(o3+        sat); %add predicted ambiguity to y0 (L1)
                    y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat)))=y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat)))+lambda(sat,2).*X_t1_t(o3+nSatTot+sat); %add predicted ambiguity to y0 (L2)
                elseif (strcmp(obs_comb,'IONO_FREE'))
                    y0_noamb(length(sat_pr)+1:end)=y0_noamb(length(sat_pr)+1:end)+X_t1_t(o3+sat); %add predicted ambiguity to y0 (IF)
                end
            else
                if (frequencies == 1)
                    y0_noamb(length(sat_pr)+1:end)=y0_noamb(length(sat_pr)+1:end)+lambda(sat,1).*X_t1_t(o3+sat); %add predicted ambiguity to y0 (L1)
                else
                    y0_noamb(length(sat_pr)+1:end)=y0_noamb(length(sat_pr)+1:end)+lambda(sat,2).*X_t1_t(o3+sat); %add predicted ambiguity to y0 (L2)
                end
            end
        end
        
        if (~isempty(sat_pr))
            if (length(frequencies) == 2 && strcmp(obs_comb,'NONE'))
                if (flag_tropo)
                    y0_noamb(1:length(sat_pr))                  = y0_noamb(1:length(sat_pr))+gmfw_R.*X_t1_t(o3+nN+1);
                    y0_noamb(length(sat_pr)+(1:length(sat_pr))) = y0_noamb(1:length(sat_pr))+gmfw_R.*X_t1_t(o3+nN+1);
                    y0_noamb(length(sat_pr)*2+            (1:length(sat))) = y0_noamb(length(sat_pr)*2+            (1:length(sat)))+gmfw_R.*X_t1_t(o3+nN+1);
                    y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat))) = y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat)))+gmfw_R.*X_t1_t(o3+nN+1);
                end
                y0_noamb(1:length(sat_pr))                  = y0_noamb(1:length(sat_pr))+X_t1_t(o3+nN+2);
                y0_noamb(length(sat_pr)+(1:length(sat_pr))) = y0_noamb(1:length(sat_pr))+X_t1_t(o3+nN+2);
                y0_noamb(length(sat_pr)*2+            (1:length(sat))) = y0_noamb(length(sat_pr)*2+            (1:length(sat)))+X_t1_t(o3+nN+2);
                y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat))) = y0_noamb(length(sat_pr)*2+length(sat)+(1:length(sat)))+X_t1_t(o3+nN+2);
            else
                if (flag_tropo)
                    y0_noamb(1:length(sat_pr))     = y0_noamb(1:length(sat_pr))    +gmfw_R.*X_t1_t(o3+nN+1);
                    y0_noamb(length(sat_pr)+1:end) = y0_noamb(length(sat_pr)+1:end)+gmfw_R.*X_t1_t(o3+nN+1);
                end
                y0_noamb(1:length(sat_pr))     = y0_noamb(1:length(sat_pr))    +X_t1_t(o3+nN+2);
                y0_noamb(length(sat_pr)+1:end) = y0_noamb(length(sat_pr)+1:end)+X_t1_t(o3+nN+2);
            end
        end

        % decomment to use only phase
%         y0_noamb=y0_noamb(length(sat_pr_np)+1:end);
%         H1=H(length(sat_pr_np)+1:end,[1 o1+1 o2+1]);
%         Cnn = Cnn(length(sat_pr_np)+1:end,length(sat_pr_np)+1:end);
%         H=H(length(sat_pr_np)+1:end,:);
%         y0=y0(length(sat_pr_np)+1:end);
%         index_residuals_outlier=[nSatTot+sat_np];

        % decomment to use only code
%         y0_noamb=y0_noamb(1:length(sat_pr_np));
%         H1=H(1:length(sat_pr_np),[1 o1+1 o2+1]);
%         Cnn = Cnn(1:length(sat_pr_np),1:length(sat_pr_np));
%         H=H(1:length(sat_pr_np),:);
%         y0=y0(1:length(sat_pr_np));
%         index_residuals_outlier=sat_pr_np;
        
        index_outlier_i=1:length(y0_noamb);

        while (search_for_outlier == 1)
            
            [index_outlier, ~, s02_ls(t)] = OLOO(H1, y0_noamb, Cnn);
            if (index_outlier ~= 0)
                %fprintf('\nOUTLIER FOUND! obs %d/%d\n',index_outlier,length(y0));
                H(index_outlier_i(index_outlier),:)   = [];
                y0(index_outlier_i(index_outlier),:)  = [];
                Cnn(index_outlier_i(index_outlier),:) = [];
                Cnn(:,index_outlier_i(index_outlier)) = [];
                y0_noamb(index_outlier,:)  = [];
                H1(index_outlier,:)  = [];
                outliers(index_residuals_outlier(index_outlier_i(index_outlier)))=1;
            else
                search_for_outlier = 0;
            end

        end
        
        %------------------------------------------------------------------------------------
        % DILUTION OF PRECISION
        %------------------------------------------------------------------------------------
        
        cov_XYZ = (alpha'*alpha)^-1;
        cov_ENU = global2localCov(cov_XYZ, XR0);
        
        PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
        HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
        VDOP = sqrt(cov_ENU(3,3));
        
        %--------------------------------------------------------------------------------------------
        % DOPPLER-BASED PREDICTION OF PHASE RANGES
        %--------------------------------------------------------------------------------------------
        doppler_pred_range1_R = zeros(nSatTot,1);
        doppler_pred_range2_R = zeros(nSatTot,1);
        if (dop1(sat))
            doppler_pred_range1_R(sat,1) = ph1(sat) - dop1(sat);
        end
        if (dop2(sat))
            doppler_pred_range2_R(sat,1) = ph2(sat) - dop2(sat);
        end
        
    else
        %to point out that notwithstanding the satellite configuration,
        %data were not analysed (motion by dynamics only).
        pivot = 0;
    end
else
    %to point out that notwithstanding the satellite configuration,
    %data were not analysed (motion by dynamics only).
    pivot = 0;
end

%----------------------------------------------------------------------------------------
% KALMAN FILTER
%----------------------------------------------------------------------------------------

%Kalman filter equations
if (nsat >= min_nsat)

    K = T*Cee*T' + Cvv;

    G = K*H' * (H*K*H' + Cnn)^(-1);

    Xhat_t_t = (I-G*H)*X_t1_t + G*y0;

    X_t1_t = T*Xhat_t_t;

    Cee = (I-G*H)*K;

else
    %positioning done only by system dynamics

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';

end

if exist('y0_residuals','var')
    if (length(frequencies) == 2 && strcmp(obs_comb,'NONE'))
        X_est = Xhat_t_t([[1 o1+1 o2+1]';o3+sat;o3+nSatTot+sat]);
        residuals_float([sat_pr;nSatTot+sat_pr;nSatTot*2+sat;nSatTot*3+sat]) = y0_residuals - H1_residuals(:,[[1 o1+1 o2+1]';o3+sat;o3+nSatTot+sat])*X_est;
    else
        X_est = Xhat_t_t([[1 o1+1 o2+1]';o3+sat;o3+nN+nT;o3+nN+nT+nC]);
        residuals_float([sat_pr;nSatTot*2+sat]) = y0_residuals - H1_residuals(:,[[1 o1+1 o2+1]';o3+sat;o3+nN+nT;o3+nN+nT+nC])*X_est;
    end
end

%--------------------------------------------------------------------------------------------
% INTEGER AMBIGUITY SOLVING BY LAMBDA METHOD
%--------------------------------------------------------------------------------------------

% if (flag_IAR && ~isempty(sat) && nsat >= min_nsat)
%     %try to solve integer ambiguities
%     [Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat)] = lambdafix(Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), Cee(o3+sat,o3+sat), Cee([1 o1+1 o2+1],o3+sat));
% else
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    succ_rate = [succ_rate NaN];
    fixed_solution = [fixed_solution 0];
% end

residuals_fixed = residuals_float;

%--------------------------------------------------------------------------------------------
% KALMAN FILTER DOP
%--------------------------------------------------------------------------------------------

%covariance propagation
Cee_XYZ = Cee([1 o1+1 o2+1],[1 o1+1 o2+1]);
Cee_ENU = global2localCov(Cee_XYZ, Xhat_t_t([1 o1+1 o2+1]));

%KF DOP computation
KPDOP = sqrt(Cee_XYZ(1,1) + Cee_XYZ(2,2) + Cee_XYZ(3,3));
KHDOP = sqrt(Cee_ENU(1,1) + Cee_ENU(2,2));
KVDOP = sqrt(Cee_ENU(3,3));

%positioning error
%sigma_rho = sqrt(Cee(1,1,end) + Cee(o1+1,o1+1,end) + Cee(o2+1,o2+1,end));
