function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop_model ...
         (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, dop1_Rsat, ...
          dop1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, dop2_Rsat, dop2_Msat, ...
          snr_R, snr_M, order, phase, dtMdot)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop_model ...
%   (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, dop1_Rsat, ...
%    dop1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, dop2_Rsat, dop2_Msat, ...
%    snr_R, snr_M, order, phase, dtMdot);
%
% INPUT:
%   pos_M = master position (X,Y,Z)
%   time = GPS time
%   Eph = satellite ephemerides
%   iono = ionospheric parameters (not used)
%   pr1_Rsat  = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr1_Msat  = MASTER-SATELLITE code pseudorange (L1 carrier)
%   ph1_Rsat  = ROVER-SATELLITE phase observation (L1 carrier)
%   ph1_Msat  = MASTER-SATELLITE phase observation (L1 carrier)
%   dop1_Rsat = ROVER-SATELLITE Doppler observation (L1 carrier)
%   dop1_Msat = MASTER-SATELLITE Doppler observation (L1 carrier)
%   pr2_Rsat  = ROVER-SATELLITE code pseudorange (L2 carrier)
%   pr2_Msat  = MASTER-SATELLITE code pseudorange (L2 carrier)
%   ph2_Rsat  = ROVER-SATELLITE phase observation (L2 carrier)
%   ph2_Msat  = MASTER-SATELLITE phase observation (L2 carrier)
%   dop2_Rsat = ROVER-SATELLITE Doppler observation (L2 carrier)
%   dop2_Msat = MASTER-SATELLITE Doppler observation (L2 carrier)
%   snr_R = signal-to-noise ratio for ROVER observations
%   snr_M = signal-to-noise ratio for MASTER observations
%   order = dynamical model order (1,2,3)
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%   dtMdot = master receiver clock drift
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
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

global lambda1 lambda2

global sigmaq_vE sigmaq_vN sigmaq_vU sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold o1 o2 o3 nN %order
global tile_header tile_georef dtm_dir
global h_antenna

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global doppler_pred_range1_R doppler_pred_range2_R
global doppler_pred_range1_M doppler_pred_range2_M

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables to point out events (satellite addition, losses, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%azimuth, elevation and ROVER-satellite and MASTER-satellite distances
azR = zeros(32,1);
azM = zeros(32,1);
elR = zeros(32,1);
elM = zeros(32,1);
distR = zeros(32,1);
distM = zeros(32,1);

%--------------------------------------------------------------------------------------------
% DYNAMIC MODEL OF THE KALMAN FILTER
%--------------------------------------------------------------------------------------------

% if the model can change then the largest set of variables is set up
o1 = 3;
o2 = 6;
o3 = 9;

%zeroes vectors useful in the matrices definition
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
%position and velocity equations
T0 = zeros(o1);
T0(1:order,1:order) = eye(order) + diag(ones(order-1,1),1);

%first degree polynomial
% T0 = [1 0 0; 0 0 0; 0 0 0];
%second degree polynomial
% T0 = [1 1 0; 0 1 0; 0 0 0];
%third degree polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%matrix structure of initial comb_N
N0 = eye(nN);

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- for the other two variables Y e Z
%comb_N(t+1) = comb_N(t)

T = [T0      Z_o1_o1 Z_o1_o1 Z_o1_nN;
     Z_o1_o1 T0      Z_o1_o1 Z_o1_nN;
     Z_o1_o1 Z_o1_o1 T0      Z_o1_nN;
     Z_nN_o1 Z_nN_o1 Z_nN_o1 N0];

%construction of an identity matrix of 38 variables (6 for position and
%velocity + 32 or 64 for the satellites number) for the further computations
I = eye(o3+nN);

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
% (if a static model is used, no noise is added)
Cvv = zeros(o3+nN);
if (order > 1)
    Cvv(order,order) = sigmaq_vE;
    Cvv(order+o1,order+o1) = sigmaq_vN;
    Cvv(order+o2,order+o2) = sigmaq_vU;

    %propagate diagonal local cov matrix to global cov matrix
    Cvv(order+[0 o1 o2],order+[0 o1 o2]) = local2globalCov(Cvv(order+[0 o1 o2],order+[0 o1 o2]), X_t1_t([1 o1+1 o2+1]));
end

%------------------------------------------------------------------------------------
% RECEIVER PREDICTED POSITION
%------------------------------------------------------------------------------------
posR_app = X_t1_t([1,o1+1,o2+1]);

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
if (length(phase) == 2)
    sat_pr = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
    sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & (ph1_Rsat ~= 0) & (ph1_Msat ~= 0) & ...
                (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) & (ph2_Rsat ~= 0) & (ph2_Msat ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) );
        sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & ...
                    (ph1_Rsat ~= 0) & (ph1_Msat ~= 0) );
    else
        sat_pr = find( (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
        sat = find( (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) & ...
                    (ph2_Rsat ~= 0) & (ph2_Msat ~= 0) );
    end
end

%only satellites with code and phase
%sat_pr = sat;

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

% if (length(sat_pr) >= 4)
%     %ROVER positioning with code double differences
%     if (phase(1) == 1)
%         [posR_app, cov_pos_R] = code_double_diff(Xhat_t_t([1 o1+1 o2+1]), pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
%     else
%         [posR_app, cov_pos_R] = code_double_diff(Xhat_t_t([1 o1+1 o2+1]), pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
%     end
% 
%     if (phase(1) == 1)
%         [posR_app, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
%     else
%         [posR_app, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
%     end
% 
% else
%     posR_app = X_t1_t([1,o1+1,o2+1]);
% end

% if (sqrt(sum((posR_app - X_t1_t([1,o1+1,o2+1])).^2))) <= 3
%     posR_app = X_t1_t([1,o1+1,o2+1]);
% end

%approximated coordinates X Y Z
X_app = posR_app(1);
Y_app = posR_app(2);
Z_app = posR_app(3);

%----------------------------------------------------------------------------------------
% CONVERSION FROM CARTESIAN TO GEODETIC COORDINATES
%----------------------------------------------------------------------------------------
[phiR_app, lamR_app, hR_app] = cart2geod(X_app, Y_app, Z_app);
[phiM, lamM, hM] = cart2geod(pos_M(1), pos_M(2), pos_M(3));

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
if ( ~isempty(tile_row) & ~isempty(tile_col) )
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
    [tile_height tile_width] = size(tile_buffer);

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
% SATELLITE ELEVATION, PIVOT AND CUT-OFF
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

posS = zeros(3,32);
dtS = zeros(32,1);
pos_S_ttime = zeros(3,32);
vel_S = zeros(3,32);
ttime = zeros(32,1);
prRS_app = zeros(32,1);
prMS_app = zeros(32,1);
err_tropo_RS = zeros(32,1);
err_tropo_MS = zeros(32,1);
err_iono_RS = zeros(32,1);
err_iono_MS = zeros(32,1);

for i = 1:size(sat_pr)

    i_sat = sat_pr(i);

    %satellite position (with clock error and Earth rotation corrections)
    [posS_tmp, dtS_tmp, pos_S_ttime_tmp, vel_S_tmp, ttime_tmp] = sat_corr(Eph, sat_pr(i), time, pr1_Rsat(sat_pr(i)));
    
    if (~isempty(posS_tmp))
        
        posS(:,i_sat) = posS_tmp;
        dtS(i_sat,1) = dtS_tmp;
        pos_S_ttime(:,i_sat) = pos_S_ttime_tmp;
        vel_S(:,i_sat) = vel_S_tmp;
        ttime(i_sat,1) = ttime_tmp;
        
        %computation of the satellite azimuth and elevation
        [azR(i_sat), elR(i_sat), distR(i_sat)] = topocent(posR_app, posS(:,i_sat)');
        [azM(i_sat), elM(i_sat), distM(i_sat)] = topocent(pos_M, posS(:,i_sat)');
        
        %computation of ROVER-SATELLITE approximated pseudorange
        prRS_app(i_sat) = sqrt(sum((posR_app - posS(:,i_sat)).^2));
        prMS_app(i_sat) = sqrt(sum((pos_M - posS(:,i_sat)).^2));
        
        %computation of tropospheric errors
        err_tropo_RS(i_sat) = err_tropo(elR(i_sat), hR_app);
        err_tropo_MS(i_sat) = err_tropo(elM(i_sat), hM);
        
        %computation of ionospheric errors
        err_iono_RS(i_sat) = err_iono(iono, phiR_app*180/pi, lamR_app*180/pi, azR(i_sat), elR(i_sat), time);
        err_iono_MS(i_sat) = err_iono(iono, phiM*180/pi, lamM*180/pi, azM(i_sat), elM(i_sat), time);
    end

    %test ephemerides availability, elevation and signal-to-noise ratio
    if (isempty(posS_tmp) | elR(i_sat) < cutoff | snr_R(i_sat) < snr_threshold)
        bad_sat(j,1) = i_sat;
        j = j + 1;
    end
    
    if (nargin > 20 & ~isempty(dtMdot) & dop1_Msat(i_sat) == 0)
        [dop1_Msat(i_sat), dop2_Msat(i_sat)] = doppler_shift_approx(pos_M, zeros(3,1), pos_S_ttime(:,i_sat), vel_S(:,i_sat), ttime(i_sat,1), dtMdot, i_sat, Eph);
    end
end

%removal of satellites without ephemerides or with elevation or SNR lower than the respective threshold
sat_pr(ismember(sat_pr,bad_sat) == 1) = [];
sat(ismember(sat,bad_sat) == 1) = [];

%previous pivot
if (pivot ~= 0)
    pivot_old = pivot;
end

%current pivot
if ~isempty(sat)
    [max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);
else
    [max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
    pivot = sat_pr(i);
end
%pivot = find(elR == max(elR));

%----------------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%----------------------------------------------------------------------------------------

%previous satellite configuration (with phase measurements)
sat_old = find(conf_sat == 1);

%satellite configuration: code only (-1), both code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slip configuration
conf_cs = zeros(32,1);

%total number of visible satellites
nsat = size(sat_pr,1);
n = nsat - 1;

%if the number of visible satellites is equal or greater than min_nsat
if (nsat >= min_nsat)
    
    %------------------------------------------------------------------------------------
    % SATELLITE ADDITION/LOSS
    %------------------------------------------------------------------------------------

    sat_dead = []; %#ok<NASGU>
    sat_born = [];

    %search for a lost satellite
    if (length(sat) < length(sat_old))

        check_off = 1;

        %save lost satellites
        sat_dead = setdiff(sat_old,sat);

        %for lost satellites it is fundamental to set their N-PIVOT
        % combinations to 0. Furthermore it could be convenient to raise
        %their uncertainty (not necessary - done when a new satellite is
        %added)
        comb_N1 = 0;
        comb_N2 = 0;

        if (length(phase) == 2)
            X_t1_t(o3+sat_dead,1) = comb_N1;
            X_t1_t(o3+32+sat_dead,1) = comb_N2;
        else
            if (phase == 1)
                X_t1_t(o3+sat_dead,1) = comb_N1;
            else
                X_t1_t(o3+sat_dead,1) = comb_N2;
            end
        end
    end

    %search for a new satellite
    if (length(sat) > length(sat_old))
        
        check_on = 1;
        
        %new satellites
        sat_born = setdiff(sat,sat_old);
        
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
                    [N1_slip, N1_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                    [N2_slip, N2_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                    %[N1_slip, N1_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                    %[N2_slip, N2_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>

                    X_t1_t(o3+sat_born,1) = N1_born;
                    X_t1_t(o3+32+sat_born,1) = N2_born;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1_born * eye(size(sat_born,1));
                    %Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq_N2_born * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                else
                    if (phase == 1)
                        [N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                        %[N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                    else
                        [N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                        %[N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                    end
                    
                    X_t1_t(o3+sat_born,1) = N_born;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
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
    if (pivot ~= pivot_old)

        check_pivot = 1;

        %matrix construction to update the PIVOT change
        %sat: vector with the current visible satellites
        %nsat: current satellites vector dimension
        R = zeros(32);
        R(sat,sat) = eye(length(sat));
        R(sat,pivot) = -1;
        R(pivot_old,pivot_old) = 0;
        R(pivot,pivot) = 0;

        I0 = eye(o3);
        Z_32_o3 = zeros(32,o3);
        Z_o3_32 = zeros(o3,32);
        Z_32_32 = zeros(32,32);

        %total matrix construction
        %sat_old, sat
        if (length(phase) == 2)
            A = [I0 Z_o3_32 Z_o3_32; Z_32_o3 R Z_32_32; Z_32_o3 Z_32_32 R];
        else
            A = [I0 Z_o3_32; Z_32_o3 R];
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
            [check_cs1, N_slip1, sat_slip1] = cycle_slip_detection(X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), prRS_app(sat), prMS_app(sat), doppler_pred_range1_R(sat), doppler_pred_range1_M(sat), pivot, sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            [check_cs2, N_slip2, sat_slip2] = cycle_slip_detection(X_t1_t(o3+33:o3+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), prRS_app(sat), prMS_app(sat), doppler_pred_range2_R(sat), doppler_pred_range2_M(sat), pivot, sat, sat_born, cs_threshold, 2); %#ok<ASGLU>

            if (check_cs1 | check_cs2)
                check_cs = 1;
            end
        else
            if (phase == 1)
                [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), prRS_app(sat), prMS_app(sat), doppler_pred_range1_R(sat), doppler_pred_range1_M(sat), pivot, sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            else
                [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o3+1:o3+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), prRS_app(sat), prMS_app(sat), doppler_pred_range2_R(sat), doppler_pred_range2_M(sat), pivot, sat, sat_born, cs_threshold, 2); %#ok<ASGLU>
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

    if (check_on | check_cs)
        if (length(phase) == 2)
            %[N1_slip, N1_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
            %[N2_slip, N2_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
            [N1_slip, N1_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr));
            [N2_slip, N2_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr));
            
            if (check_on)
                X_t1_t(o3+sat_born,1) = N1_born;
                X_t1_t(o3+32+sat_born,1) = N2_born;
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1_born * eye(size(sat_born,1));
                %Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq_N2_born * eye(size(sat_born,1));
                Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
            
            if (check_cs1)
                conf_cs(sat_slip1) = 1;
                X_t1_t(o3+sat_slip1) = N1_slip;
                Cvv(o3+sat_slip1,o3+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
            end
            
            if (check_cs2)
                conf_cs(sat_slip2) = 1;
                X_t1_t(o3+32+sat_slip2) = N2_slip;
                Cvv(o3+32+sat_slip2,o3+32+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
            end
        else
            if (phase == 1)
                [N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
                %[N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr));
            else
                [N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
                %[N_slip, N_born] = ambiguity_init(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o3+sat_pr));
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
    end
    
    %------------------------------------------------------------------------------------
    % OBSERVATION EQUATIONS
    %------------------------------------------------------------------------------------

    %rows in which the phase observation is available
    p = find(ismember(setdiff(sat_pr,pivot),setdiff(sat,pivot))==1);

    %function that calculates the Kalman filter parameters
    [alfa1, prstim_pr1, prstim_ph1, ddc1, ddp1] = input_kalman(posR_app, posS(:,sat_pr), prRS_app(sat_pr), prMS_app(sat_pr), pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Msat(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_MS(sat_pr), sat_pr, pivot, 1);
    [alfa2, prstim_pr2, prstim_ph2, ddc2, ddp2] = input_kalman(posR_app, posS(:,sat_pr), prRS_app(sat_pr), prMS_app(sat_pr), pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Msat(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_MS(sat_pr), sat_pr, pivot, 2);

    %zeroes vector useful in matrix definitions
    Z_1_nN = zeros(1,nN);
    Z_n_nN = zeros(n,nN);
    Z_n_om = zeros(n,o1-1);
    Z_1_om = zeros(1,o1-1);

    %H matrix computation for the code
    H_cod1 = [alfa1(:,1) Z_n_om alfa1(:,2) Z_n_om alfa1(:,3) Z_n_om Z_n_nN];
    H_cod2 = [alfa2(:,1) Z_n_om alfa2(:,2) Z_n_om alfa2(:,3) Z_n_om Z_n_nN];
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
    L_fas1 = zeros(n,32);
    L_fas2 = zeros(n,32);
    v = 1;
    for u = 1 : n+1 % with the pivot
        if (sat_pr(u) ~= pivot)
            L_fas1(v,sat_pr(u)) = -(lambda1);
            L_fas2(v,sat_pr(u)) = -(lambda2);
            v = v+1;
        end
    end

    %H matrix computation for the phase
    if ~isempty(p)
        H_fas1 = [alfa1(p,1) Z_n_om(p,:) alfa1(p,2) Z_n_om(p,:) alfa1(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
        H_fas2 = [alfa2(p,1) Z_n_om(p,:) alfa2(p,2) Z_n_om(p,:) alfa2(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
        if (length(phase) == 2)
            H_fas1(:,o3+1:o3+32) = L_fas1(p,:);
            H_fas2(:,o3+33:o3+64) = L_fas2(p,:);
            H_fas = [H_fas1; H_fas2];
        else
            if (phase == 1)
                H_fas1(:,o3+1:o3+32) = L_fas1(p,:);
                H_fas = H_fas1;
            else
                H_fas2(:,o3+1:o3+32) = L_fas2(p,:);
                H_fas = H_fas2;
            end
        end
    else
        H_fas = [];
    end

    %H matrix computation for the DTM pseudo-observation
    H_dtm = [];
    if (h_dtm ~= tile_header.nodata)
        H_dtm = [cos(phiR_app)*cos(lamR_app) Z_1_om cos(phiR_app)*sin(lamR_app) Z_1_om sin(phiR_app) Z_1_om Z_1_nN];
    end

    %construction of the complete H matrix
    H = [H_cod; H_fas; H_dtm];

    %Y0 vector computation for the code
    y0_cod1 = ddc1 - prstim_pr1 + alfa1(:,1)*X_app + alfa1(:,2)*Y_app + alfa1(:,3)*Z_app;
    y0_cod2 = ddc2 - prstim_pr2 + alfa2(:,1)*X_app + alfa2(:,2)*Y_app + alfa2(:,3)*Z_app;

    %Y0 vector computation for the phase
    if ~isempty(p)
        y0_fas1 = ddp1(p) - prstim_ph1(p) + alfa1(p,1)*X_app + alfa1(p,2)*Y_app + alfa1(p,3)*Z_app;
        y0_fas2 = ddp2(p) - prstim_ph2(p) + alfa2(p,1)*X_app + alfa2(p,2)*Y_app + alfa2(p,3)*Z_app;
    else
        y0_fas1 = [];
        y0_fas2 = [];
    end

    %Y0 vector computation for DTM constrain
    y0_dtm = [];
    if (h_dtm ~= tile_header.nodata)
        y0_dtm = h_dtm  - hR_app + cos(phiR_app)*cos(lamR_app)*X_app + cos(phiR_app)*sin(lamR_app)*Y_app + sin(phiR_app)*Z_app;
    end

    %construction of the total Y0 vector
    if (length(phase) == 2)
        y0_cod = [y0_cod1; y0_cod2];
        y0_fas = [y0_fas1; y0_fas2];
    else
        if (phase == 1)
            y0_cod = y0_cod1;
            y0_fas = y0_fas1;
        else
            y0_cod = y0_cod2;
            y0_fas = y0_fas2;
        end
    end
    y0 = [y0_cod; y0_fas; y0_dtm];

    %------------------------------------------------------------------------------------
    % OBSERVATION COVARIANCE MATRIX
    %------------------------------------------------------------------------------------

    %construction of the cofactor matrix
    Q = cofactor_matrix(elR(sat_pr), elM(sat_pr), snr_R(sat_pr), snr_M(sat_pr), sat_pr, pivot);

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
    % DILUTION OF PRECISION
    %------------------------------------------------------------------------------------
    
    cov_XYZ = (alfa1'*alfa1)^-1;
    cov_ENU = global2localCov(cov_XYZ, posR_app);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));

    %--------------------------------------------------------------------------------------------
    % DOPPLER-BASED PREDICTION OF PHASE RANGES
    %--------------------------------------------------------------------------------------------
    doppler_pred_range1_R = zeros(32,1);
    doppler_pred_range2_R = zeros(32,1);
    doppler_pred_range1_M = zeros(32,1);
    doppler_pred_range2_M = zeros(32,1);
    if (dop1_Rsat(sat))
        doppler_pred_range1_R(sat,1) = ph1_Rsat(sat) - dop1_Rsat(sat);
    end
    if (dop2_Rsat(sat))
        doppler_pred_range2_R(sat,1) = ph2_Rsat(sat) - dop2_Rsat(sat);
    end
    if (dop1_Msat(sat))
        doppler_pred_range1_M(sat,1) = ph1_Msat(sat) - dop1_Msat(sat);
    end
    if (dop2_Msat(sat))
        doppler_pred_range2_M(sat,1) = ph2_Msat(sat) - dop2_Msat(sat);
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
    %positioning done only by the system dynamics

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';
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
