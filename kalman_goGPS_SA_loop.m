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
%                           goGPS v0.1 beta
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
global X_comb
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global flag_LS_N_estim
global sm_weight

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
    sat_pr = find( (pr1_Rsat(:,2) ~= 0) & (pr2_Rsat(:,2) ~= 0) );
    sat = find( (pr1_Rsat(:,2) ~= 0) & (ph1_Rsat(:,2) ~= 0) & ...
                (pr2_Rsat(:,2) ~= 0) & (ph2_Rsat(:,2) ~= 0) );
else
    if (phase == 1)
        sat_pr = find(pr1_Rsat(:,2) ~= 0);
        sat = find( (pr1_Rsat(:,2) ~= 0) & (ph1_Rsat(:,2) ~= 0) );
    else
        sat_pr = find(pr2_Rsat(:,2) ~= 0);
        sat = find( (pr2_Rsat(:,2) ~= 0) & (ph2_Rsat(:,2) ~= 0) );
    end
end

%only satellites with code and phase
%sat_pr = sat;

%satellites that can be phase-smoothed
if (length(phase) == 2)
    sat_sm = find( (pr1_Rsat(:,1) ~= 0) & (pr2_Rsat(:,1) ~= 0) & (ph1_Rsat(:,1) ~= 0) & (ph2_Rsat(:,1) ~= 0) & ...
                   (pr1_Rsat(:,2) ~= 0) & (pr2_Rsat(:,2) ~= 0) & (ph1_Rsat(:,2) ~= 0) & (ph2_Rsat(:,2) ~= 0));
else
    if (phase == 1)
        sat_sm = find( (pr1_Rsat(:,1) ~= 0) & (ph1_Rsat(:,1) ~= 0) & ...
                       (pr1_Rsat(:,2) ~= 0) & (ph1_Rsat(:,2) ~= 0));
    else
        sat_sm = find( (pr2_Rsat(:,1) ~= 0) & (ph2_Rsat(:,1) ~= 0) & ...
                       (pr2_Rsat(:,2) ~= 0) & (ph2_Rsat(:,2) ~= 0));
    end
end

if (sm_weight ~= 0)
    sm_weight = sm_weight - 0.01;
end

%computation of phase-smoothed code
for i = 1:size(sat_sm)
    if (phase(1) == 1)
        R_ex = pr1_Rsat(sat_sm(i),1)/lambda1 + (ph1_Rsat(sat_sm(i),2) - ph1_Rsat(sat_sm(i),1));
        pr1_Rsat(sat_sm(i),2) = (sm_weight * (pr1_Rsat(sat_sm(i),2)/lambda1) + (1-sm_weight) * R_ex) * lambda1;
    else
        R_ex = pr2_Rsat(sat_sm(i),1)/lambda2 + (ph2_Rsat(sat_sm(i),2) - ph2_Rsat(sat_sm(i),1));
        pr2_Rsat(sat_sm(i),2) = (sm_weight * (pr2_Rsat(sat_sm(i),2)/lambda2) + (1-sm_weight) * R_ex) * lambda2;
    end
end

%drop the measurements of the previous epoch
pr1_Rsat(:,1) = [];
ph1_Rsat(:,1) = [];
pr2_Rsat(:,1) = [];
ph2_Rsat(:,1) = [];

%------------------------------------------------------------------------------------
% SATELLITE ELEVATION, PIVOT AND CUT-OFF
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

for i = 1:size(sat_pr)
    
    %satellite position correction (clock and Earth rotation)
    Rot_X = sat_corr(Eph, sat_pr(i), time, pr1_Rsat(i), X_t1_t([1,o1+1,o2+1])');
    
    if (~isempty(Rot_X))
        %azimuth, elevation, ROVER-SATELLITE distance computation
        [azR(sat_pr(i)), elR(sat_pr(i)), distR(sat_pr(i))] = topocent(X_t1_t([1,o1+1,o2+1]), Rot_X');
    end
    
    %test ephemerides availability, elevation and signal-to-noise ratio
    if (isempty(Rot_X) | elR(sat_pr(i)) < cutoff | snr_R(sat_pr(i)) < snr_threshold)
        bad_sat(j,1) = sat_pr(i);
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

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

% if (length(sat_pr) >= 4)
% 
%     %ROVER positioning with code double differences
%     if (phase(1) == 1)
%         [pos_R, cov_pos_R] = code_SA(Xhat_t_t([1,o1+1,o2+1]), pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     else
%         [pos_R, cov_pos_R] = code_SA(Xhat_t_t([1,o1+1,o2+1]), pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     end
% 
%     if (phase(1) == 1)
%         [pos_R, cov_pos_R] = code_SA(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     else
%         [pos_R, cov_pos_R] = code_SA(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph); %#ok<NASGU>
%     end
% 
% else
%     pos_R = X_t1_t([1,o1+1,o2+1]);
% end

% if (sqrt(sum((pos_R - X_t1_t([1,o1+1,o2+1])).^2))) <= 3
    pos_R = X_t1_t([1,o1+1,o2+1]);
% end    

%approximated coordinates X Y Z
X_app = pos_R(1);
Y_app = pos_R(2);
Z_app = pos_R(3);

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
% OBSERVATION EQUATIONS
%------------------------------------------------------------------------------------

%if the number of visible satellites is equal or greater than min_nsat
if (nsat >= min_nsat)

    %rows in which the phase observation is available
    p = find(ismember(sat_pr,sat)==1);

    %function that calculates the Kalman filter parameters
    [alfa1, prstim1, err_iono] = input_kalman_SA(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), time, sat_pr, Eph, iono);
    [alfa2, prstim2          ] = input_kalman_SA(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), time, sat_pr, Eph, iono);

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
    for u = 1 : n
        L_fas1(u,sat_pr(u)) = -(lambda1);
        L_fas2(u,sat_pr(u)) = -(lambda2);
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
        H_dtm = [cos(phi_app)*cos(lam_app) Z_1_om cos(phi_app)*sin(lam_app) Z_1_om sin(phi_app) Z_1_om Z_1_nN];
    end

    %construction of the complete H matrix
    H = [H_cod; H_fas; H_dtm];

    %Y0 vector computation for the code
    y0_cod1 = pr1_Rsat(sat_pr) - prstim1 + alfa1(:,1)*X_app + alfa1(:,2)*Y_app + alfa1(:,3)*Z_app;
    y0_cod2 = pr2_Rsat(sat_pr) - prstim2 + alfa2(:,1)*X_app + alfa2(:,2)*Y_app + alfa2(:,3)*Z_app;

    %Y0 vector computation for the phase
    if ~isempty(p)
        y0_fas1 = ph1_Rsat(sat)*lambda1 - (prstim1(p) - 2*err_iono(p)') + alfa1(p,1)*X_app + alfa1(p,2)*Y_app + alfa1(p,3)*Z_app;
        y0_fas2 = ph2_Rsat(sat)*lambda2 - (prstim2(p) - 2*err_iono(p)') + alfa2(p,1)*X_app + alfa2(p,2)*Y_app + alfa2(p,3)*Z_app;
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

    cov_XYZ = (alfa1'*alfa1)^-1;
    cov_ENU = global2localCov(cov_XYZ, pos_R);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));

    %------------------------------------------------------------------------------------
    % SATELLITE ADDITION/LOSS
    %------------------------------------------------------------------------------------

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

        %do not use least squares ambiguity estimation
        % NOTE: LS amb. estimation is automatically switched off if the number of
        % satellites with phase available is not sufficient
        if(~flag_LS_N_estim) | (size(sat) < 4)
            check_on = 1;
            
            %new satellites
            sat_born = setdiff(sat,sat_old);
            
            %sigmaq_pos_R = diag(T*Cee*T');
            %sigmaq_pos_R = sigmaq_pos_R([1,o1+1,o2+1]);
            
            %N combination estimation
            [N1, sigmaq_N1] = amb_estimate_observ_SA(pr1_Rsat(sat), ph1_Rsat(sat), 1); %#ok<NASGU>
            [N2, sigmaq_N2] = amb_estimate_observ_SA(pr2_Rsat(sat), ph2_Rsat(sat), 2); %#ok<NASGU>
            %[N1, sigmaq_N1] = amb_estimate_approx_SA(pos_R, sigmaq_pos_R, pr1_Rsat(sat), ph1_Rsat(sat), Eph, time, sat, 1); %#ok<NASGU>
            %[N2, sigmaq_N2] = amb_estimate_approx_SA(pos_R, sigmaq_pos_R, pr2_Rsat(sat), ph2_Rsat(sat), Eph, time, sat, 2); %#ok<NASGU>

            index = find(ismember(sat,sat_born) == 0);
            N1(index) = [];
            N2(index) = [];
            
            %estimated parameters
            if (length(phase) == 2)
                X_t1_t(o3+sat_born,1) = N1;
                X_t1_t(o3+32+sat_born,1) = N2;
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
                %Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
                Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            else
                if (phase == 1)
                    X_t1_t(o3+sat_born,1) = N1;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                    %Cvv(o3+(1:32),o3+(1:32)) = sigmaq0_N * eye(32);
                else
                    X_t1_t(o3+sat_born,1) = N2;
                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
                    Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                end
            end

        %use least squares ambiguity estimation
        else
            
            check_on = 1;
            
            %new satellites
            sat_born = setdiff(sat,sat_old);
            
            %N combination estimation (least squares)
            [null_pos_R, null_cov_pos_R, N1_stim, cov_N1_stim] = code_phase_SA(pos_R, pr1_Rsat(sat), ph1_Rsat(sat), snr_R(sat), sat, time, Eph, 1, iono); %#ok<ASGLU>
            [null_pos_R, null_cov_pos_R, N2_stim, cov_N2_stim] = code_phase_SA(pos_R, pr2_Rsat(sat), ph2_Rsat(sat), snr_R(sat), sat, time, Eph, 2, iono); %#ok<ASGLU>
            
            if isempty(cov_N1_stim) %if it was not possible to compute the covariance matrix
                cov_N1_stim = sigmaq0_N * eye(length(sat));
            end
            if isempty(cov_N2_stim) %if it was not possible to compute the covariance matrix
                cov_N2_stim = sigmaq0_N * eye(length(sat));
            end
            sigmaq_N1 = diag(cov_N1_stim);
            sigmaq_N2 = diag(cov_N2_stim);
            
            index = find(ismember(sat,sat_born) == 0);
            N1_stim(index)   = [];
            N2_stim(index)   = [];
            sigmaq_N1(index) = [];
            sigmaq_N2(index) = [];
            
            %save estimated parameters
            if (length(phase) == 2)
                X_t1_t(o3+sat_born,1) = N1_stim;
                X_t1_t(o3+32+sat_born,1) = N2_stim;
                Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N1);
                Cvv(o3+32+sat_born,o3+32+sat_born) = diag(sigmaq_N2);
            else
                if (phase == 1)
                    X_t1_t(o3+sat_born,1) = N1_stim;
                    Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N1);
                else
                    X_t1_t(o3+sat_born,1) = N2_stim;
                    Cvv(o3+sat_born,o3+sat_born) = diag(sigmaq_N2);
                end
            end
        end
    end

    %------------------------------------------------------------------------------------
    % PIVOT CHANGE
    %------------------------------------------------------------------------------------

    %search for a possible PIVOT change
    if (pivot ~= pivot_old)

        check_pivot = 1;

    end

    %------------------------------------------------------------------------------------
    % CYCLE-SLIP
    %------------------------------------------------------------------------------------
    
    if ~isempty(sat)
        
        %do not use least squares ambiguity estimation
        % NOTE: LS amb. estimation is automatically switched off if the number of
        % satellites with phase available is not sufficient
        if(~flag_LS_N_estim) | (size(sat) < 4)

            %Test presence/absence of a cycle-slip at the current epoch.
            %The state of the system is changed only for phase ambiguities
            if (length(phase) == 2)
                [cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_kalman_SA(pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), pr1_Rsat(sat), Eph, time, sat, 3, 1);
                [cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_kalman_SA(pos_R, X_t1_t(o3+33:o3+64), ph2_Rsat(sat), pr2_Rsat(sat), Eph, time, sat, 3, 2);
                
                if (cycle_slip_found1 == 1)
                    check_cs = 1;
                    conf_cs(sat_slip1) = 1;
                    X_t1_t(o3+sat_slip1) = N_slip1;
                    Cvv(o3+sat_slip1,o3+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
                end
                if (cycle_slip_found2 == 1)
                    check_cs = 1;
                    conf_cs(sat_slip2) = 1;
                    X_t1_t(o3+32+sat_slip2) = N_slip2;
                    Cvv(o3+32+sat_slip2,o3+32+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
                end
            else
                if (phase == 1)
                    [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman_SA(pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), pr1_Rsat(sat), Eph, time, sat, 3, 1);
                else
                    [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman_SA(pos_R, X_t1_t(o3+1:o3+32), ph2_Rsat(sat), pr2_Rsat(sat), Eph, time, sat, 3, 2);
                end
                if (cycle_slip_found == 1)
                    check_cs = 1;
                    conf_cs(sat_slip) = 1;
                    X_t1_t(o3+sat_slip) = N_slip;
                    Cvv(o3+sat_slip,o3+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
                end
            end
            
        %use least squares ambiguity estimation
        else
            %Test presence/absence of a cycle-slip at the current epoch.
            %The state of the system is changed only for phase ambiguities
            if (length(phase) == 2)
                
                [cycle_slip_found1, N_slip1, sat_slip1, sigmaq_N_slip1] = cycle_slip_LS_N_SA(pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), pr1_Rsat(sat), snr_R(sat), Eph, time, sat, iono, cs_threshold, 1);
                [cycle_slip_found2, N_slip2, sat_slip2, sigmaq_N_slip2] = cycle_slip_LS_N_SA(pos_R, X_t1_t(o3+33:o3+64), ph2_Rsat(sat), pr2_Rsat(sat), snr_R(sat), Eph, time, sat, iono, cs_threshold, 2);
                
                if (cycle_slip_found1 == 1)
                    check_cs = 1;
                    conf_cs(sat_slip1) = 1;
                    X_t1_t(o3+sat_slip1) = N_slip1;
                    Cvv(o3+sat_slip1,o3+sat_slip1) = diag(sigmaq_N_slip1);
                end
                if (cycle_slip_found2 == 1)
                    check_cs = 1;
                    conf_cs(sat_slip2) = 1;
                    X_t1_t(o3+32+sat_slip2) = N_slip2;
                    Cvv(o3+32+sat_slip2,o3+32+sat_slip2) = diag(sigmaq_N_slip2);
                end
            else
                if (phase == 1)
                    [cycle_slip_found, N_slip, sat_slip, sigmaq_N_slip] = cycle_slip_LS_N_SA(pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), pr1_Rsat(sat), snr_R(sat), Eph, time, sat, iono, cs_threshold, 1);
                else
                    [cycle_slip_found, N_slip, sat_slip, sigmaq_N_slip] = cycle_slip_LS_N_SA(pos_R, X_t1_t(o3+1:o3+32), ph2_Rsat(sat), pr2_Rsat(sat), snr_R(sat), Eph, time, sat, iono, cs_threshold, 2);
                end
                if (cycle_slip_found == 1)
                    check_cs = 1;
                    conf_cs(sat_slip) = 1;
                    X_t1_t(o3+sat_slip) = N_slip;
                    Cvv(o3+sat_slip,o3+sat_slip) = diag(sigmaq_N_slip);
                end
            end
        end
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

%--------------------------------------------------------------------------------------------
% STATIC POSITIONING
%--------------------------------------------------------------------------------------------

if o1 == 1
    X_comb = Xhat_t_t(1:3,1);
end