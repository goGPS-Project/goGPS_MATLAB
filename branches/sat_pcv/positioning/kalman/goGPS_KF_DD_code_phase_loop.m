function [check_on, check_off, check_pivot, check_cs] = goGPS_KF_DD_code_phase_loop ...
         (XM, time_rx, pr1_R, pr1_M, ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
         dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, dtMdot, flag_IAR, antenna_PCV, sbas)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = goGPS_DD_code_phase_loop ...
%        (XM, time_rx, pr1_R, pr1_M, ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%        dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, dtMdot, flag_IAR, antenna_PCV, sbas);
%
% INPUT:
%   XM = master position (X,Y,Z)
%   time_rx = GPS time
%   pr1_R  = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr1_M  = MASTER-SATELLITE code pseudorange (L1 carrier)
%   ph1_R  = ROVER-SATELLITE phase observation (L1 carrier)
%   ph1_M  = MASTER-SATELLITE phase observation (L1 carrier)
%   dop1_R = ROVER-SATELLITE Doppler observation (L1 carrier)
%   dop1_M = MASTER-SATELLITE Doppler observation (L1 carrier)
%   pr2_R  = ROVER-SATELLITE code pseudorange (L2 carrier)
%   pr2_M  = MASTER-SATELLITE code pseudorange (L2 carrier)
%   ph2_R  = ROVER-SATELLITE phase observation (L2 carrier)
%   ph2_M  = MASTER-SATELLITE phase observation (L2 carrier)
%   dop2_R = ROVER-SATELLITE Doppler observation (L2 carrier)
%   dop2_M = MASTER-SATELLITE Doppler observation (L2 carrier)
%   snr_R  = signal-to-noise ratio for ROVER observations
%   snr_M  = signal-to-noise ratio for MASTER observations
%   Eph = satellite ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono =  ionospheric parameters (vector of zeroes if not available)
%   lambda = wavelength matrix (depending on the enabled constellations)
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%   dtMdot = master receiver clock drift
%   flag_IAR = boolean variable to enable/disable integer ambiguity resolution
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
%   Addition and loss of satellites, cycle slips e pivot changes are considered.

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

%--------------------------------------------------------------------------------------------
% KALMAN FILTER PARAMETERS
%--------------------------------------------------------------------------------------------

global sigmaq_vE sigmaq_vN sigmaq_vU %sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold o1 o2 o3 nN
global tile_header tile_georef dtm_dir
global h_antenna

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global doppler_pred_range1_R doppler_pred_range2_R
global doppler_pred_range1_M doppler_pred_range2_M
global ratiotest mutest succ_rate fixed_solution

global t residuals_fixed residuals_float outliers s02_ls
%global min_ambfixRMS min_ambfloatRMS

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables to point out events (satellite addition, losses, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

%azimuth, elevation and ROVER-satellite and MASTER-satellite distances
azR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
% (if a static model is used, no noise is added)
Cvv = zeros(o3+nN);
if (o1 > 1)
    Cvv(o1,o1) = sigmaq_vE;
    Cvv(o2,o2) = sigmaq_vN;
    Cvv(o3,o3) = sigmaq_vU;

    %propagate diagonal local cov matrix to global cov matrix
    Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), X_t1_t([1 o1+1 o2+1]));
end

%----------------------------------------------------------------------------------------
% POSITION ESTIMATION ERROR VARIANCE
%----------------------------------------------------------------------------------------

sigmaq_pos_R = diag(T*Cee*T');
sigmaq_pos_R = sigmaq_pos_R([1,o1+1,o2+1]);

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
if (length(phase) == 2)
    sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (pr2_R ~= 0) & (pr2_M ~= 0) );
    sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (ph1_R ~= 0) & (ph1_M ~= 0) & ...
                (pr2_R ~= 0) & (pr2_M ~= 0) & (ph2_R ~= 0) & (ph2_M ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
        sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                    (ph1_R ~= 0) & (ph1_M ~= 0) );
    else
        sat_pr = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
        sat = find( (pr2_R ~= 0) & (pr2_M ~= 0) & ...
                    (ph2_R ~= 0) & (ph2_M ~= 0) );
    end
end
sat_pr = sat_pr(ismember(sat_pr, Eph(30,:)));
sat = sat(ismember(sat, Eph(30,:)));

%only satellites with code and phase
%sat_pr = sat;

%previous satellite configuration (with phase measurements)
sat_old = find(conf_sat == 1);

%number of visible satellites
nsat = size(sat_pr,1);

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

%approximate position
XR0 = X_t1_t([1,o1+1,o2+1]);
%VR0 = X_t1_t([2,o1+2,o2+2]);
flag_XR = 2;

%approximated coordinates X Y Z
X_app = XR0(1);
Y_app = XR0(2);
Z_app = XR0(3);

%----------------------------------------------------------------------------------------
% CONVERSION FROM CARTESIAN TO GEODETIC COORDINATES
%----------------------------------------------------------------------------------------
[phiR_app, lamR_app, hR_app] = cart2geod(X_app, Y_app, Z_app);
%[phiM, lamM, hM] = cart2geod(XM(1), XM(2), XM(3));

%----------------------------------------------------------------------------------------
% EXTRACTION OF THE HEIGHT PSEUDO-OBSERVATION FROM THE DTM
%----------------------------------------------------------------------------------------

%projection to UTM coordinates
[E_app, N_app] = geod2plan(phiR_app, lamR_app);

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

%if the number of visible satellites is equal or greater than min_nsat
if (nsat >= min_nsat)
    
    %----------------------------------------------------------------------------------------
    % SATELLITE POSITION COMPUTATION
    %----------------------------------------------------------------------------------------
    
    sat_pr_nocutoff = sat_pr;
    
    if (phase(1) == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono1_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM] = init_positioning(time_rx, pr1_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, sbas,  XM,  [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, phase,       2, 0); %#ok<NASGU,ASGLU>
        [ ~, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono1_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR] = init_positioning(time_rx, pr1_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, sbas,  XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, flag_XR, 1); %#ok<NASGU,ASGLU>
        
        err_iono2_M = err_iono1_M .* ionoFactor(sat_pr_M,2);
        err_iono2_R = err_iono1_R .* ionoFactor(sat_pr_R,2);
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono2_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM] = init_positioning(time_rx, pr2_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, sbas,  XM,  [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, phase,       2, 0); %#ok<NASGU,ASGLU>
        [ ~, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono2_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR] = init_positioning(time_rx, pr2_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, sbas,  XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, flag_XR, 1); %#ok<NASGU,ASGLU>
        
        err_iono1_M = err_iono2_M ./ ionoFactor(sat_pr_M,2);
        err_iono1_R = err_iono2_R ./ ionoFactor(sat_pr_R,2);
    end
    
    %keep only satellites that rover and master have in common
    [sat_pr, iR, iM] = intersect(sat_pr_R, sat_pr_M);
    err_tropo_R = err_tropo_R(iR);
    err_tropo_M = err_tropo_M(iM);
    err_iono1_R = err_iono1_R(iR);
    err_iono1_M = err_iono1_M(iM);
    err_iono2_R = err_iono2_R(iR);
    err_iono2_M = err_iono2_M(iM);
    
    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_nocutoff, sat_pr);
    sat(ismember(sat,sat_removed)) = [];
    
    for i = 1:size(sat_pr)
        if (nargin > 22 && ~isempty(dtMdot) && dop1_M(sat_pr(i)) == 0 && sum(sum(Eph)) ~= 0)
            [dop1_M(sat_pr(i)), dop2_M(sat_pr(i))] = doppler_shift_approx(XM, zeros(3,1), XS_tx(i,:)', VS_tx(i,:)', time_tx(i), dtMdot, sat_pr(i), Eph, lambda(sat_pr(i),:));
        end
    end
    
    %----------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
    %----------------------------------------------------------------------------------------
    
    %satellite configuration: code only (-1), both code and phase (+1);
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr) = -1;
    conf_sat(sat) = +1;
    
    %cycle-slip configuration
    conf_cs = zeros(nSatTot,1);
    
    %number of visible satellites (after the cutoffs)
    nsat = size(sat_pr,1);
    n = nsat - 1;
    
    %previous pivot
    if (pivot ~= 0)
        pivot_old = pivot;
    end
    
    %current pivot
    if ~isempty(sat)
        [null_max_elR, pivot_index] = max(elR(sat)); %#ok<ASGLU>
        pivot = sat(pivot_index);
    else
        [null_max_elR, pivot_index] = max(elR(sat_pr)); %#ok<ASGLU>
        pivot = sat_pr(pivot_index);
    end
    
    %if the number of available satellites after the cutoffs is equal or greater than min_nsat
    if (nsat >= min_nsat)
        
        %------------------------------------------------------------------------------------
        % SATELLITE ADDITION/LOSS
        %------------------------------------------------------------------------------------
        
        %search for a lost satellite
        sat_dead = setdiff(sat_old,sat);
        
        if (~isempty(sat_dead))
            
            check_off = 1;

            %for lost satellites it is fundamental to set their N-PIVOT
            % combinations to 0. Furthermore it could be convenient to raise
            %their uncertainty (not necessary - done when a new satellite is
            %added)
            N1 = 0;
            N2 = 0;
            
            if (length(phase) == 2)
                X_t1_t(o3+sat_dead,1) = N1;
                X_t1_t(o3+nSatTot+sat_dead,1) = N2;
            else
                if (phase == 1)
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
            
            %if a new satellite is going to be the pivot, its ambiguity needs
            %to be estimated before applying the pivot change
            if ~isempty(find(sat_born == pivot, 1))
                %if it is not the only satellite with phase
                if (length(sat) > 1)
                    %if the former pivot is still among satellites with phase
                    if ~isempty(find(sat == pivot_old, 1))
                        %set the old pivot as temporary pivot
                        pivot_tmp = pivot_old;
                    else
                        %find the best candidate as temporary pivot
                        sat_tmp = setdiff(sat,pivot);
                        [max_elR, i] = max(elR(sat_tmp)); %#ok<ASGLU>
                        pivot_tmp = sat_tmp(i);
                        %reset the ambiguities of other satellites according to the temporary pivot
                        sat_born = sat;
                        pivot_old = pivot_tmp;
                    end
                    sat_born = setdiff(sat_born,pivot_tmp);
                    sat_slip1 = [];
                    sat_slip2 = [];
                    sat_slip = [];
                    if (length(phase) == 2)
                        [N1_slip, N1_born, sigmaq_N1_slip, sigmaq_N1_born] = ambiguity_init(XR0, XS, pr1_R(sat_pr), pr1_M(sat_pr), ph1_R(sat_pr), ph1_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono1_R, err_iono1_M, pivot_tmp, lambda(sat_pr,1), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R); %#ok<ASGLU>
                        [N2_slip, N2_born, sigmaq_N2_slip, sigmaq_N2_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono2_R, err_iono2_M, pivot_tmp, lambda(sat_pr,2), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R); %#ok<ASGLU>
                        
                        X_t1_t(o3+sat_born,1) = N1_born;
                        X_t1_t(o3+nSatTot+sat_born,1) = N2_born;
                        Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N1_born);
                        Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = diag(sigmaq_N2_born);
                        %Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                        %Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    else
                        if (phase == 1)
                            [N_slip, N_born, sigmaq_N_slip, sigmaq_N_born] = ambiguity_init(XR0, XS, pr1_R(sat_pr), pr1_M(sat_pr), ph1_R(sat_pr), ph1_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono1_R, err_iono1_M, pivot_tmp, lambda(sat_pr,1), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R); %#ok<ASGLU>
                        else
                            [N_slip, N_born, sigmaq_N_slip, sigmaq_N_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono2_R, err_iono2_M, pivot_tmp, lambda(sat_pr,2), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R); %#ok<ASGLU>
                        end
                        
                        X_t1_t(o3+sat_born,1) = N_born;
                        Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N_born);
                        %Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    end
                end
                sat_born = setdiff(sat_born,pivot);
                check_on = 0;
            end
        end
        
        %------------------------------------------------------------------------------------
        % PIVOT CHANGE
        %------------------------------------------------------------------------------------
        
        %search for a possible PIVOT change
        if (pivot ~= pivot_old && pivot_old ~= 0)
            
            check_pivot = 1;
            
            %matrix construction to update the PIVOT change
            %sat: vector with the current visible satellites
            %nsat: current satellites vector dimension
            R = zeros(nSatTot);
            R(sat,sat) = eye(length(sat));
            R(sat,pivot) = -1;
            R(pivot_old,pivot_old) = 0;
            R(pivot,pivot) = 0;
            
            I0 = eye(o3);
            Z_ns_o3 = zeros(nSatTot,o3);
            Z_o3_ns = zeros(o3,nSatTot);
            Z_ns_ns = zeros(nSatTot,nSatTot);
            
            %total matrix construction
            %sat_old, sat
            if (length(phase) == 2)
                A = [I0 Z_o3_ns Z_o3_ns; Z_ns_o3 R Z_ns_ns; Z_ns_o3 Z_ns_ns R];
            else
                A = [I0 Z_o3_ns; Z_ns_o3 R];
            end
            
            %new state estimate
            X_t1_t = A*X_t1_t;
            
            %re-computation of the Cee covariance matrix at the previous epoch
            Cee = A*Cee*A';
        end
        
        %------------------------------------------------------------------------------------
        % CYCLE-SLIP
        %------------------------------------------------------------------------------------
        
        if ~isempty(sat)
            
            %Test presence/absence of a cycle-slip at the current epoch.
            %The state of the system is not changed yet
            if (length(phase) == 2)
                [check_cs1, N_slip1, sat_slip1] = cycle_slip_detection(X_t1_t(o3+1:o3+nSatTot), ph1_R(sat), ph1_M(sat), distR(sat), distM(sat), doppler_pred_range1_R(sat), doppler_pred_range1_M(sat), pivot, sat, sat_born, cs_threshold, lambda(sat,1)); %#ok<ASGLU>
                [check_cs2, N_slip2, sat_slip2] = cycle_slip_detection(X_t1_t(o3+nSatTot+1:o3+nSatTot*2), ph2_R(sat), ph2_M(sat), distR(sat), distM(sat), doppler_pred_range2_R(sat), doppler_pred_range2_M(sat), pivot, sat, sat_born, cs_threshold, lambda(sat,2)); %#ok<ASGLU>
                
                if (check_cs1 || check_cs2)
                    check_cs = 1;
                end
            else
                if (phase == 1)
                    [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o3+1:o3+nSatTot), ph1_R(sat), ph1_M(sat), distR(sat), distM(sat), doppler_pred_range1_R(sat), doppler_pred_range1_M(sat), pivot, sat, sat_born, cs_threshold, lambda(sat,1)); %#ok<ASGLU>
                else
                    [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o3+1:o3+nSatTot), ph2_R(sat), ph2_M(sat), distR(sat), distM(sat), doppler_pred_range2_R(sat), doppler_pred_range2_M(sat), pivot, sat, sat_born, cs_threshold, lambda(sat,2)); %#ok<ASGLU>
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
        
        if (check_on || check_cs)
            if (length(phase) == 2)
                [N1_slip, N1_born, sigmaq_N1_slip, sigmaq_N1_born] = ambiguity_init(XR0, XS, pr1_R(sat_pr), pr1_M(sat_pr), ph1_R(sat_pr), ph1_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono1_R, err_iono1_M, pivot, lambda(sat_pr,1), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R);
                [N2_slip, N2_born, sigmaq_N2_slip, sigmaq_N2_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono2_R, err_iono2_M, pivot, lambda(sat_pr,2), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R);
                
                if (check_on)
                    X_t1_t(o3+sat_born,1) = N1_born;
                    X_t1_t(o3+nSatTot+sat_born,1) = N2_born;
                    Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N1_born);
                    Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = diag(sigmaq_N2_born);
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    %Cvv(o3+nSatTot+sat_born,o3+nSatTot+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                end
                
                if (check_cs1)
                    conf_cs(sat_slip1) = 1;
                    X_t1_t(o3+sat_slip1) = N1_slip;
                    Cvv(o3+sat_slip1,o3+sat_slip1) = diag(sigmaq_N1_slip);
                    %Cvv(o3+sat_slip1,o3+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
                end
                
                if (check_cs2)
                    conf_cs(sat_slip2) = 1;
                    X_t1_t(o3+nSatTot+sat_slip2) = N2_slip;
                    Cvv(o3+nSatTot+sat_slip2,o3+nSatTot+sat_slip2) = diag(sigmaq_N2_slip);
                    %Cvv(o3+nSatTot+sat_slip2,o3+nSatTot+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
                end
            else
                if (phase == 1)
                    [N_slip, N_born, sigmaq_N_slip, sigmaq_N_born] = ambiguity_init(XR0, XS, pr1_R(sat_pr), pr1_M(sat_pr), ph1_R(sat_pr), ph1_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono1_R, err_iono1_M, pivot, lambda(sat_pr,1), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R);
                else
                    [N_slip, N_born, sigmaq_N_slip, sigmaq_N_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono2_R, err_iono2_M, pivot, lambda(sat_pr,2), X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr), sigmaq_pos_R);
                end
                
                if (check_on)
                    X_t1_t(o3+sat_born,1) = N_born;                    
                    Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N_born);
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                end
                
                if (check_cs)
                    conf_cs(sat_slip) = 1;
                    X_t1_t(o3+sat_slip) = N_slip;
                    Cvv(o3+sat_slip,o3+sat_slip) = diag(sigmaq_N_slip);
                    %Cvv(o3+sat_slip,o3+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
                end
            end
        end

        %------------------------------------------------------------------------------------
        % OBSERVATION EQUATIONS
        %------------------------------------------------------------------------------------
        
        %rows in which the phase observation is available
        p = find(ismember(setdiff(sat_pr,pivot),setdiff(sat,pivot))==1);
  
        %compute PCV: phase and code 1 
        [~, index_ph]=intersect(sat_pr,sat);   
        
        if (~isempty(antenna_PCV) && antenna_PCV(2).n_frequency ~= 0) % master
            index_master=2; 
            PCV1_M=PCV_interp(antenna_PCV(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 1);
            pr1_M(sat_pr)=pr1_M(sat_pr)-PCV1_M;
            ph1_M(sat)=ph1_M(sat)-PCV1_M(index_ph)./lambda(sat,1);
            
%             if (length(phase) == 2)
%                 PCV2_M=PCV_interp(antenna_PCV(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 2, XR0, XS);
%                 pr2_M(sat_pr)=pr2_M(sat_pr)-PCV2_M;
%                 ph2_M(sat)=ph2_M(sat)-PCV2_M(index_ph)./lambda(sat,1);
%             end
        end
        
        if (~isempty(antenna_PCV) && antenna_PCV(1).n_frequency ~= 0) % rover
            index_rover=1; 
            PCV1_R=PCV_interp(antenna_PCV(index_rover), 90-elR(sat_pr), azR(sat_pr), sys, 1);
            pr1_R(sat_pr)=pr1_R(sat_pr)-PCV1_R;
            ph1_R(sat)=ph1_R(sat)-PCV1_R(index_ph)./lambda(sat,1);
            
%             if (length(phase) == 2)
%                 PCV2_M=PCV_interp(antenna_PCV(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 2, XR0, XS);
%                 pr2_M(sat_pr)=pr2_M(sat_pr)-PCV2_M;
%                 ph2_M(sat)=ph2_M(sat)-PCV2_M(index_ph)./lambda(sat,1);
%             end
        end        
        
        
        %function that calculates the Kalman filter parameters
        [alpha, probs_pr1, probs_ph1, prapp_pr1, prapp_ph1, probs_pr2, probs_ph2, prapp_pr2, prapp_ph2] = input_kalman(XR0, XS, pr1_R(sat_pr), ph1_R(sat_pr), pr1_M(sat_pr), ph1_M(sat_pr), pr2_R(sat_pr), ph2_R(sat_pr), pr2_M(sat_pr), ph2_M(sat_pr), err_tropo_R, err_iono1_R, err_iono2_R, err_tropo_M, err_iono1_M, err_iono2_M, distR(sat_pr), distM(sat_pr), sat_pr, pivot, lambda(sat_pr,:));

        %zeroes vector useful in matrix definitions
        Z_1_nN = zeros(1,nN);
        Z_n_nN = zeros(n,nN);
        Z_n_om = zeros(n,o1-1);
        Z_1_om = zeros(1,o1-1);
        
        %H matrix computation for the code
        H_cod1 = [alpha(:,1) Z_n_om alpha(:,2) Z_n_om alpha(:,3) Z_n_om Z_n_nN];
        H_cod2 = [alpha(:,1) Z_n_om alpha(:,2) Z_n_om alpha(:,3) Z_n_om Z_n_nN];
        if (length(phase) == 2)
            H_cod = [H_cod1; H_cod2];
        else
            if (phase == 1)
                H_cod = H_cod1;
            else
                H_cod = H_cod2;
            end
        end
        
        %lambda positions computation
        L_pha1 = zeros(n,nSatTot);
        L_pha2 = zeros(n,nSatTot);
        v = 1;
        for u = 1 : n+1 % with the pivot
            if (sat_pr(u) ~= pivot)
                L_pha1(v,sat_pr(u)) = -(lambda(sat_pr(u),1));
                L_pha2(v,sat_pr(u)) = -(lambda(sat_pr(u),2));
                v = v+1;
            end
        end
        
        %H matrix computation for the phase
        if ~isempty(p)
            H_pha1 = [alpha(p,1) Z_n_om(p,:) alpha(p,2) Z_n_om(p,:) alpha(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
            H_pha2 = [alpha(p,1) Z_n_om(p,:) alpha(p,2) Z_n_om(p,:) alpha(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
            if (length(phase) == 2)
                H_pha1(:,o3+1:o3+nSatTot) = L_pha1(p,:);
                H_pha2(:,o3+nSatTot+1:o3+nSatTot*2) = L_pha2(p,:);
                H_pha = [H_pha1; H_pha2];
            else
                if (phase == 1)
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
            H_dtm = [cos(phiR_app)*cos(lamR_app) Z_1_om cos(phiR_app)*sin(lamR_app) Z_1_om sin(phiR_app) Z_1_om Z_1_nN];
        end
        
        %construction of the complete H matrix
        H = [H_cod; H_pha; H_dtm];
        
        %Y0 vector computation for the code
        y0_cod1 = probs_pr1 - prapp_pr1 + alpha(:,1)*X_app + alpha(:,2)*Y_app + alpha(:,3)*Z_app;
        y0_cod2 = probs_pr2 - prapp_pr2 + alpha(:,1)*X_app + alpha(:,2)*Y_app + alpha(:,3)*Z_app;
        
        %Y0 vector computation for the phase
        if ~isempty(p)
            y0_pha1 = probs_ph1(p) - prapp_ph1(p) + alpha(p,1)*X_app + alpha(p,2)*Y_app + alpha(p,3)*Z_app;
            y0_pha2 = probs_ph2(p) - prapp_ph2(p) + alpha(p,1)*X_app + alpha(p,2)*Y_app + alpha(p,3)*Z_app;
        else
            y0_pha1 = [];
            y0_pha2 = [];
        end
        
        %Y0 vector computation for DTM constraint
        y0_dtm = [];
        if (h_dtm ~= tile_header.nodata)
            y0_dtm = h_dtm  - hR_app + cos(phiR_app)*cos(lamR_app)*X_app + cos(phiR_app)*sin(lamR_app)*Y_app + sin(phiR_app)*Z_app;
        end
        
        %construction of the total Y0 vector
        if (length(phase) == 2)
            y0_cod = [y0_cod1; y0_cod2];
            y0_pha = [y0_pha1; y0_pha2];
        else
            if (phase == 1)
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
        Q = cofactor_matrix(elR(sat_pr), elM(sat_pr), snr_R(sat_pr), snr_M(sat_pr), pivot_index);
        
        %zeroes vector useful in matrix definitions
        Z_n_n = zeros(n,n);
        
        %multiplication by the code variance and the phase variance to build the matrix
        if ~isempty(p)
            if (length(phase) == 2)
                Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,:) Z_n_n(:,p) Z_n_n(:,p); Z_n_n(:,:) sigmaq_cod2*Q(:,:) Z_n_n(:,p) Z_n_n(:,p);
                    Z_n_n(p,:) Z_n_n(p,:) sigmaq_ph*Q(p,p) Z_n_n(p,p); Z_n_n(p,:) Z_n_n(p,:) Z_n_n(p,p) sigmaq_ph*Q(p,p)];
            else
                if (phase == 1)
                    Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
                else
                    Cnn = [sigmaq_cod2*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
                end
            end
        else
            if (length(phase) == 2)
                Cnn = [sigmaq_cod1*Q Z_n_n; Z_n_n sigmaq_cod2*Q];
            else
                if (phase == 1)
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
        
        sat_np = sat(sat~=pivot);
        sat_pr_np = sat_pr(sat_pr~=pivot);
        
        y0_residuals=y0;
        H1_residuals=H;
        index_residuals_outlier=[sat_pr_np;nSatTot+sat_np];  %[code;phase]
                
        y0_noamb=y0;
        if (~isempty(sat_np))
            y0_noamb(length(sat_pr_np)+1:end)=y0_noamb(length(sat_pr_np)+1:end)+lambda(sat_np,1).*X_t1_t(o3+sat_np); %add predicted ambiguity to y0
        end
        H1=H(:,[1 o1+1 o2+1]);

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
        doppler_pred_range1_M = zeros(nSatTot,1);
        doppler_pred_range2_M = zeros(nSatTot,1);
        if (dop1_R(sat))
            doppler_pred_range1_R(sat,1) = ph1_R(sat) - dop1_R(sat);
        end
        if (dop2_R(sat))
            doppler_pred_range2_R(sat,1) = ph2_R(sat) - dop2_R(sat);
        end
        if (dop1_M(sat))
            doppler_pred_range1_M(sat,1) = ph1_M(sat) - dop1_M(sat);
        end
        if (dop2_M(sat))
            doppler_pred_range2_M(sat,1) = ph2_M(sat) - dop2_M(sat);
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
    
    %min_ambfloatRMS(t,1) = min(sqrt(diag(Cee(o3+sat_np,o3+sat_np)))); %#ok<NASGU>
    
    Xhat_t_t = (I-G*H)*X_t1_t + G*y0;

    X_t1_t = T*Xhat_t_t;

    Cee = (I-G*H)*K;
    
    %sat_np = sat(sat ~= pivot);
    
else
    %positioning done only by the system dynamics

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';
end

if exist('y0_residuals','var') && exist('sat_np','var')
    X_est = Xhat_t_t([[1 o1+1 o2+1]';o3+sat_np]);
    residuals_float([sat_pr_np;nSatTot+sat_np]) = y0_residuals - H1_residuals(:,[[1 o1+1 o2+1]';o3+sat_np])*X_est;
end

%--------------------------------------------------------------------------------------------
% INTEGER AMBIGUITY SOLVING BY LAMBDA METHOD
%--------------------------------------------------------------------------------------------

if (flag_IAR)
    sat_np = sat(sat ~= pivot);
    if (~isempty(sat_np))
        pos_amb_zero = find(Xhat_t_t(o3+sat_np) == 0);
        if (~isempty(pos_amb_zero))
            sat_np(pos_amb_zero) = [];
        end
        pos_cov_zero = find(sum(abs(Cee(o3+sat_np,o3+sat_np)),1) == 0);
        if (~isempty(pos_cov_zero))
            sat_np(pos_cov_zero) = [];
        end
    end
end

if (flag_IAR && ~isempty(sat_np) && nsat >= min_nsat)
    %try to solve integer ambiguities
    [Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat_np), varNfix, varPosfix] = lambdafix(Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat_np), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), Cee(o3+sat_np,o3+sat_np), Cee([1 o1+1 o2+1],o3+sat_np)); %#ok<ASGLU,NASGU>
        
    %min_ambfixRMS(t,1) = min(sqrt(diag(varNfix)));
else
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    succ_rate = [succ_rate NaN];
    fixed_solution = [fixed_solution 0];
end

if exist('y0_residuals','var') 
    if exist('pos_cov_zero' , 'var') && (~isempty(pos_cov_zero))
            y0_residuals(length(sat_pr_np)+pos_cov_zero) = [];
            H1_residuals(length(sat_pr_np)+pos_cov_zero,:) = [];
            %y0_residuals(pos_cov_zero) = []; % why also on codes?
            %H1_residuals(pos_cov_zero,:) = [];
    end
  
    X_est = Xhat_t_t([[1 o1+1 o2+1]';o3+sat_np]);
    residuals_fixed([sat_pr_np;nSatTot+sat_np]) = y0_residuals - H1_residuals(:,[[1 o1+1 o2+1]';o3+sat_np])*X_est;
end

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
%sigma_rho = sqrt(Cee(1,1) + Cee(o1+1,o1+1) + Cee(o2+1,o2+1));
