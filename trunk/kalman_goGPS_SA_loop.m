function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_SA_loop(time, Eph, iono, pr1_Rsat, ph1_Rsat, pr2_Rsat, ph2_Rsat, snr_R, phase)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_SA_loop(time, Eph, iono, pr1_Rsat, ph1_Rsat, pr2_Rsat, ph2_Rsat, snr_R, phase);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemerides
%   iono = ionospheric parameters
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (L1 carrier)
%   ph1_Rsat = ROVER-SATELLITE phase observation (carrier L1)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (L2 carrier)
%   ph2_Rsat = ROVER-SATELLITE phase observation (carrier L2)
%   snr_R = signal-to-noise ratio for ROVER observations
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
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
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

global lambda1 lambda2

global sigmaq_velx sigmaq_vely sigmaq_velz sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold o1 o2 o3 nN
global tile_header tile_georef dtm_dir
global h_antenna

global Xhat_t_t X_t1_t T I Cee nsat conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP

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

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
Cvv = zeros(o3+nN);
Cvv(o1,o1) = sigmaq_velx;
Cvv(o2,o2) = sigmaq_vely;
Cvv(o3,o3) = sigmaq_velz;

if o1 == 1
    Cvv = zeros(o3+nN);
end
%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
if (length(phase) == 2)
    sat_pr = find( (pr1_Rsat(:,1) ~= 0) & (pr2_Rsat(:,1) ~= 0) );
    sat = find( (pr1_Rsat(:,1) ~= 0) & (ph1_Rsat(:,1) ~= 0) & ...
                (pr2_Rsat(:,1) ~= 0) & (ph2_Rsat(:,1) ~= 0) );
else
    if (phase == 1)
        sat_pr = find(pr1_Rsat(:,1) ~= 0);
        sat = find( (pr1_Rsat(:,1) ~= 0) & (ph1_Rsat(:,1) ~= 0) );
    else
        sat_pr = find(pr2_Rsat(:,1) ~= 0);
        sat = find( (pr2_Rsat(:,1) ~= 0) & (ph2_Rsat(:,1) ~= 0) );
    end
end

%only satellites with code and phase
%sat_pr = sat;

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

% if (length(sat_pr) >= 4)
% 
%     %ROVER positioning with code double differences
%     if (phase(1) == 1)
%         posR_app = code_SA(Xhat_t_t([1,o1+1,o2+1]), pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     else
%         posR_app = code_SA(Xhat_t_t([1,o1+1,o2+1]), pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     end
% 
%     if (phase(1) == 1)
%         posR_app = code_SA(posR_app, pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     else
%         posR_app = code_SA(posR_app, pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     end
% 
% else
%     posR_app = X_t1_t([1,o1+1,o2+1]);
% end

% if (sqrt(sum((posR_app - X_t1_t([1,o1+1,o2+1])).^2))) <= 3
    posR_app = X_t1_t([1,o1+1,o2+1]);
% end    

%approximated coordinates X Y Z
X_app = posR_app(1);
Y_app = posR_app(2);
Z_app = posR_app(3);

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

posS = zeros(32,3);
dtS = zeros(32,1);
prRS_app = zeros(32,1);
err_tropo_RS = zeros(32,1);
err_iono_RS = zeros(32,1);

for i = 1:size(sat_pr)
    
    i_sat = sat_pr(i);

    %satellite position (with clock error and Earth rotation corrections)
    [posS(i_sat,:) dtS(i_sat)] = sat_corr(Eph, i_sat, time, pr1_Rsat(i_sat), posR_app);
    
    if (~isempty(posS))

        %computation of the satellite azimuth and elevation
        [azR(i_sat), elR(i_sat), distR(i_sat)] = topocent(posR_app, posS(i_sat,:));
        
        %computation of ROVER-SATELLITE approximated pseudorange
        prRS_app(i_sat) = sqrt(sum((posR_app - posS(i_sat,:)').^2));
        
        %computation of tropospheric errors
        err_tropo_RS(i_sat) = err_tropo(elR(i_sat), h_app);
        
        %computation of ionospheric errors
        err_iono_RS(i_sat) = err_iono(iono, phi_app*180/pi, lam_app*180/pi, azR(i_sat), elR(i_sat), time);
    end
    
    %test ephemerides availability, elevation and signal-to-noise ratio
    if (isempty(posS) | elR(i_sat) < cutoff | snr_R(i_sat) < snr_threshold)
        bad_sat(j,1) = i_sat;
        j = j + 1;
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

%satellite configuration
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slip configuration
conf_cs = zeros(32,1);

%number of visible satellites
nsat = size(sat_pr,1);
n = nsat;

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

        %for lost satellites it is fundamental to set them to 0.
        %Furthermore it could be convenient to raise
        %their uncertainty (not necessary - done when a new satellite is
        %added)
        N1 = 0;
        N2 = 0;

        if (length(phase) == 2)
            X_t1_t(o3+sat_dead,1) = N1;
            X_t1_t(o3+32+sat_dead,1) = N2;
        else
            if (phase == 1)
                X_t1_t(o3+sat_dead,1) = N1;
            else
                X_t1_t(o3+sat_dead,1) = N2;
            end
        end
    end

    %search for a new satellite
    if (length(sat) > length(sat_old))

        check_on = 1;

        %new satellites
        sat_born = setdiff(sat,sat_old);
    end

    %------------------------------------------------------------------------------------
    % PIVOT CHANGE (NOT REQUIRED FOR STANDALONE)
    %------------------------------------------------------------------------------------

    %search for a possible PIVOT change
    if (pivot ~= pivot_old)

        check_pivot = 1;

    end

    %------------------------------------------------------------------------------------
    % CYCLE-SLIP
    %------------------------------------------------------------------------------------
    
    if ~isempty(sat)
        
        %Test presence/absence of a cycle-slip at the current epoch.
        %The state of the system is changed only for phase ambiguities
        if (length(phase) == 2)
            [check_cs1, N_slip1, sat_slip1] = cycle_slip_kalman_SA(X_t1_t(o3+1:o3+32), pr1_Rsat(sat), ph1_Rsat(sat), err_iono_RS, sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            [check_cs2, N_slip2, sat_slip2] = cycle_slip_kalman_SA(X_t1_t(o3+33:o3+64), pr2_Rsat(sat), ph2_Rsat(sat), (lambda2/lambda1)^2 * err_iono_RS, sat, sat_born, cs_threshold, 2); %#ok<ASGLU>
            
            if (check_cs1 | check_cs2)
                check_cs = 1;
            end
        else
            if (phase == 1)
                [check_cs, N_slip, sat_slip] = cycle_slip_kalman_SA(X_t1_t(o3+1:o3+32), pr1_Rsat(sat), ph1_Rsat(sat), err_iono_RS(sat), sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            else
                [check_cs, N_slip, sat_slip] = cycle_slip_kalman_SA(X_t1_t(o3+33:o3+64), pr2_Rsat(sat), ph2_Rsat(sat), (lambda2/lambda1)^2 * err_iono_RS(sat), sat, sat_born, cs_threshold, 2); %#ok<ASGLU>
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
    % PHASE AMBIGUITY AND RECEIVER CLOCK ERROR ESTIMATION
    %------------------------------------------------------------------------------------
    
    if (length(phase) == 2)
        [N1_slip, N1_born, dtR1] = amb_estimate_LS_SA(posR_app, posS(sat_pr,:), dtS(sat_pr), pr1_Rat(sat_pr), ph1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, sat_slip1, sat_born, prRS_app(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
        [N2_slip, N2_born, dtR2] = amb_estimate_LS_SA(posR_app, posS(sat_pr,:), dtS(sat_pr), pr2_Rat(sat_pr), ph2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, sat_slip2, sat_born, prRS_app(sat_pr), err_tropo_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<NASGU>
        
        %choose one of the two estimates
        dtR = dtR1;
        
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
            [N_slip, N_born, dtR] = amb_estimate_LS_SA(posR_app, posS(sat_pr,:), dtS(sat_pr), pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), snr_R(sat_pr), elR(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
        else
            [N_slip, N_born, dtR] = amb_estimate_LS_SA(posR_app, posS(sat_pr,:), dtS(sat_pr), pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), snr_R(sat_pr), elR(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), err_tropo_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr));
        end
        
        if (check_on)
            X_t1_t(o3+sat_born,1) = N_born;
            %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
            Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
        end
        
        if (check_cs == 1)
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

    %function that calculates the Kalman filter parameters
    [alfa, prstim_pr1, prstim_ph1, prstim_pr2, prstim_ph2] = input_kalman_SA(posR_app, posS(sat_pr,:), prRS_app(sat_pr), dtR, dtS(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr));

    %zeroes vector useful in matrix definitions
    Z_1_nN = zeros(1,nN);
    Z_n_nN = zeros(n,nN);
    Z_n_om = zeros(n,o1-1);
    Z_1_om = zeros(1,o1-1);

    %H matrix computation for the code
    H_cod1 = [alfa(:,1) Z_n_om alfa(:,2) Z_n_om alfa(:,3) Z_n_om Z_n_nN];
    H_cod2 = [alfa(:,1) Z_n_om alfa(:,2) Z_n_om alfa(:,3) Z_n_om Z_n_nN];
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
    for u = 1 : n
        L_fas1(u,sat_pr(u)) = -(lambda1);
        L_fas2(u,sat_pr(u)) = -(lambda2);
    end

    %H matrix computation for the phase
    if ~isempty(p)
        H_fas1 = [alfa(p,1) Z_n_om(p,:) alfa(p,2) Z_n_om(p,:) alfa(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
        H_fas2 = [alfa(p,1) Z_n_om(p,:) alfa(p,2) Z_n_om(p,:) alfa(p,3) Z_n_om(p,:) Z_n_nN(p,:)];
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
        H_dtm = [cos(phi_app)*cos(lam_app) Z_1_om cos(phi_app)*sin(lam_app) Z_1_om sin(phi_app) Z_1_om Z_1_nN];
    end

    %construction of the complete H matrix
    H = [H_cod; H_fas; H_dtm];

    %Y0 vector computation for the code
    y0_cod1 = pr1_Rsat(sat_pr) - prstim_pr1 + alfa(:,1)*X_app + alfa(:,2)*Y_app + alfa(:,3)*Z_app;
    y0_cod2 = pr2_Rsat(sat_pr) - prstim_pr2 + alfa(:,1)*X_app + alfa(:,2)*Y_app + alfa(:,3)*Z_app;

    %Y0 vector computation for the phase
    if ~isempty(p)
        y0_fas1 = ph1_Rsat(sat)*lambda1 - prstim_ph1(p) + alfa(p,1)*X_app + alfa(p,2)*Y_app + alfa(p,3)*Z_app;
        y0_fas2 = ph2_Rsat(sat)*lambda2 - prstim_ph2(p) + alfa(p,1)*X_app + alfa(p,2)*Y_app + alfa(p,3)*Z_app;
    else
        y0_fas1 = [];
        y0_fas2 = [];
    end

    %Y0 vector computation for DTM constrain
    y0_dtm = [];
    if (h_dtm ~= tile_header.nodata)
        y0_dtm = h_dtm + cos(phi_app)*cos(lam_app)*X_app + cos(phi_app)*sin(lam_app)*Y_app + sin(phi_app)*Z_app - h_app;
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
    Q = cofactor_matrix_SA(elR(sat_pr), snr_R(sat_pr), sat_pr);

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

    cov_XYZ = (alfa'*alfa)^-1;
    cov_ENU = global2localCov(cov_XYZ, posR_app);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));

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

%   vvX = Xhat_t_t(2,end);
%   vvY = Xhat_t_t(o1+2,end);
%   vvZ = Xhat_t_t(o2+2,end);
%   vvv = sqrt(vvX(end)^2 + vvY(end)^2 + vvZ(end)^2);

%positioning error
%sigma_rho = sqrt(Cee(1,1,end) + Cee(o1+1,o1+1,end) + Cee(o2+1,o2+1,end));
