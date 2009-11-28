function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop ...
         (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, ...
         pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R, snr_M, phase) %#ok<INUSL>

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop ...
%   (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, ...
%   pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R, snr_M, phase);
%
% INPUT:
%   pos_M = master position (X,Y,Z)
%   time = GPS time 
%   Eph = satellite ephemerides
%   iono = ionospheric parameters (not used)
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr1_Msat = MASTER-SATELLITE code pseudorange (L1 carrier)
%   ph1_Rsat = ROVER-SATELLITE phase observation (L1 carrier)
%   ph1_Msat = MASTER-SATELLITE phase observation (L1 carrier)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (L2 carrier)
%   pr2_Msat = MASTER-SATELLITE code pseudorange (L2 carrier)
%   ph2_Rsat = ROVER-SATELLITE phase observation (L2 carrier)
%   ph2_Msat = MASTER-SATELLITE phase observation (L2 carrier)
%   snr_R = signal-to-noise ratio for ROVER observations
%   snr_M = signal-to-noise ratio for MASTER observations
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
%   Addition and loss of satellites, cycle slips e pivot changes are considered.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
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

%--------------------------------------------------------------------------------------------
% KALMAN FILTER PARAMETERS
%--------------------------------------------------------------------------------------------

global lambda1 lambda2 a f

global sigmaq_velx sigmaq_vely sigmaq_velz sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq_dtm weights
global min_nsat cutoff snr_threshold o1 o2 o3 nN
global tile_header tile_georef dtm_dir
global h_antenna

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM 

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
% SATELLITE ELEVATION, PIVOT AND CUT-OFF
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

for i = 1:size(sat_pr)

    %satellite position correction (clock and Earth rotation)
    Rot_X = sat_corr(Eph, sat_pr(i), time, pr1_Rsat(sat_pr(i)), X_t1_t([1,o1+1,o2+1])');

    %azimuth, elevation, ROVER-SATELLITE distance computation
    [azR(sat_pr(i)), elR(sat_pr(i)), distR(sat_pr(i))] = topocent(X_t1_t([1,o1+1,o2+1]), Rot_X', a, f);

    %azimuth, elevation, MASTER-SATELLITE distance computation
    [azM(sat_pr(i)), elM(sat_pr(i)), distM(sat_pr(i))] = topocent(pos_M, Rot_X', a, f);

    %test on elevation and on signal-to-noise ratio
    if (elR(sat_pr(i)) < cutoff) | (snr_R(sat_pr(i)) < snr_threshold)
        bad_sat(j,1) = sat_pr(i);
        j = j + 1;
    end

end

%removal of satellites with elevation or SNR lower than the respective threshold
sat_pr(ismember(sat_pr,bad_sat) == 1) = [];
sat(ismember(sat,bad_sat) == 1) = [];

%previous pivot 
if (pivot ~= 0)
    pivot_old = pivot;
end

%current pivot
if ~isempty(sat)
    [~, i] = max(elR(sat));
    pivot = sat(i);
else
    [~, i] = max(elR(sat_pr));
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

%------------------------------------------------------------------------------------
% LINEARIZATION POINT (APPROXIMATE COORDINATES)
%------------------------------------------------------------------------------------

if (length(sat_pr) >= 4)

    pos_R = X_t1_t([1,o1+1,o2+1]);

%     %ROVER positioning with code double differences
%     if (phase(1) == 1)
%         [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph); %#ok<NASGU>
%     else
%         [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph); %#ok<NASGU>
%     end
% 
%     if (phase(1) == 1)
%         [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph); %#ok<NASGU>
%     else
%         [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph); %#ok<NASGU>
%     end

else
    pos_R = X_t1_t([1,o1+1,o2+1]);
end

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
    [h_dtm] = dtm_bilin_interp(E_app, N_app, tile_buffer, tile_header.ncols*3, tile_header.nrows*3, tile_header.cellsize, Ell, Nll, tile_header.nodata);
    
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
    p = find(ismember(setdiff(sat_pr,pivot),setdiff(sat,pivot))==1);

    %function that calculates the Kalman filter parameters
    [alfa1, prstim1, ddc1, ddp1] = input_kalman (pos_R, pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), ph1_Msat(sat_pr), time, sat_pr, pivot, Eph, 1);
    [alfa2, prstim2, ddc2, ddp2] = input_kalman (pos_R, pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), ph2_Msat(sat_pr), time, sat_pr, pivot, Eph, 2);
    %[alfa1, prstim1, ddc1, ddp1] = input_kalman (X_t1_t([1,o1+1,o2+1]), pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), ph1_Msat(sat_pr), time, sat_pr, pivot, Eph, 1);
    %[alfa2, prstim2, ddc2, ddp2] = input_kalman (X_t1_t([1,o1+1,o2+1]), pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), ph2_Msat(sat_pr), time, sat_pr, pivot, Eph, 2);

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
        H_dtm = [cos(phi_app)*cos(lam_app) Z_1_om cos(phi_app)*sin(lam_app) Z_1_om sin(phi_app) Z_1_om Z_1_nN];
    end

    %construction of the complete H matrix
    H = [H_cod; H_fas; H_dtm];

    %Y0 vector computation for the code
    y0_cod1 = ddc1 - prstim1 + alfa1(:,1)*X_app + alfa1(:,2)*Y_app + alfa1(:,3)*Z_app;
    y0_cod2 = ddc2 - prstim2 + alfa2(:,1)*X_app + alfa2(:,2)*Y_app + alfa2(:,3)*Z_app;

    %Y0 vector computation for the phase
    if ~isempty(p)
        y0_fas1 = ddp1(p) - prstim1(p) + alfa1(p,1)*X_app + alfa1(p,2)*Y_app + alfa1(p,3)*Z_app;
        y0_fas2 = ddp2(p) - prstim2(p) + alfa2(p,1)*X_app + alfa2(p,2)*Y_app + alfa2(p,3)*Z_app;
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

    %zeroes vector useful in matrix definitions
    Z_n_n = zeros(n,n);
    
    %weight function parameters
    snr_a = 30;
    snr_0 = 10;
    snr_1 = 50;
    snr_A = 30;

    if (weights == 0)

        %measurament noise covariance matrix
        %code-code or phase-phase co-factor matrix Q construction
        Q = 2*ones(n) + 2*eye(n);

    else
        if (weights == 1)

            %weight vectors (elevation)
            q_R = 1 ./ (sin(elR * pi/180).^2);
            q_M = 1 ./ (sin(elM * pi/180).^2);

        elseif (weights == 2)

            %weight vectors (signal-to-noise ratio)
            q_R = 10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1);
            q_R(snr_R >= snr_1) = 1;
            q_M = 10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1);
            q_M(snr_M >= snr_1) = 1;

        elseif (weights == 3)

            %weight vectors (elevation and signal-to-noise ratio)
            q_R = 1 ./ (sin(elR * pi/180).^2) .* (10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1));
            q_R(snr_R >= snr_1) = 1;
            q_M = 1 ./ (sin(elM * pi/180).^2) .* (10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1));
            q_M(snr_M >= snr_1) = 1;

        end

        q_RP = q_R(pivot,1);                  % ROVER-PIVOT
        q_MP = q_M(pivot,1);                  % MASTER-PIVOT
        q_RS = q_R(sat_pr);                   % ROVER-generic satellite 
        q_MS = q_M(sat_pr);                   % MASTER-generic satellite 
        q_RS(sat_pr==pivot) = [];       % ROVER-generic satellite (without pivot)
        q_MS(sat_pr==pivot) = [];       % MASTER-generic satellite (without pivot)

        %measurement noise covariance matrix
        %code-code or phase-phase co-factor matrix Q construction
        Q = (q_RP + q_MP) * ones(n) + diag(q_RS + q_MS);
    end

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
    % SATELLITE ADDITION/LOSS
    %------------------------------------------------------------------------------------

    %search for a lost satellite
    if (length(sat) < length(sat_old))

        check_off = 1;

        %save lost satellites
        sat_dead = setdiff(sat_old,sat);

        %for a lost satellite is necessary to set its N-Pivot combination to 0
        %It is useful to increase its uncertainty (not necessary - done when a new satellite is added)
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

        %display lost satellites
        ['Lost satellites at time ' num2str(time) ': ' num2str(sat_dead')]; %#ok<VUNUS>
    end

    %search for a new satellite
    if (length(sat) > length(sat_old))

        check_on = 1;

        %save new satellites
        sat_born = setdiff(sat,sat_old);

        %for a new satellites it is necessary to estimate its N-pivot combination
        %It's necessary to increase its uncertainty because only after 10 epochs 
		%it's possible to have a correct estimation
        sigmaq_pos_R = diag(T*Cee*T');
        sigmaq_pos_R = sigmaq_pos_R([1,o1+1,o2+1]);

        %N combination estimation
        [comb_N1, sigmaq_N1] = amb_estimate_approx(pos_R, pos_M, sigmaq_pos_R, pr1_Rsat(sat), pr1_Msat(sat), ph1_Rsat(sat), ph1_Msat(sat), Eph, time, pivot, sat, 1); %#ok<NASGU>
        [comb_N2, sigmaq_N2] = amb_estimate_approx(pos_R, pos_M, sigmaq_pos_R, pr2_Rsat(sat), pr2_Msat(sat), ph2_Rsat(sat), ph2_Msat(sat), Eph, time, pivot, sat, 2); %#ok<NASGU>
        %[comb_N1, sigmaq_N1] = amb_estimate_approx(X_t1_t([1,o1+1,o2+1]), pos_M, sigmaq_pos_R, pr1_Rsat(sat), pr1_Msat(sat), ph1_Rsat(sat), ph1_Msat(sat), Eph, time, pivot, sat, 1);
        %[comb_N2, sigmaq_N2] = amb_estimate_approx(X_t1_t([1,o1+1,o2+1]), pos_M, sigmaq_pos_R, pr2_Rsat(sat), pr2_Msat(sat), ph2_Rsat(sat), ph2_Msat(sat), Eph, time, pivot, sat, 2);

        index = find(ismember(sat,sat_born) == 0);
        comb_N1(index) = [];
        comb_N2(index) = [];

        %save estimated parameters
        if (length(phase) == 2)
            X_t1_t(o3+sat_born,1) = comb_N1;
            X_t1_t(o3+32+sat_born,1) = comb_N2;
            %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
            %Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
            Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
        else
            if (phase == 1)
                X_t1_t(o3+sat_born,1) = comb_N1;
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                Cvv(o3+(1:32),o3+(1:32)) = sigmaq0_N * eye(32);
            else
                X_t1_t(o3+sat_born,1) = comb_N2;
                %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
                Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
        end

        %display new satellites
        ['New satellites at time ' num2str(time) ': ' num2str(sat_born')]; %#ok<VUNUS>
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

        %display the PIVOT change
        ['PIVOT change at time ' num2str(time) ' from ' num2str(pivot_old) ' to ' num2str(pivot)]; %#ok<VUNUS>
    end

    %------------------------------------------------------------------------------------
    % CYCLE-SLIP
    %------------------------------------------------------------------------------------

    if ~isempty(sat)

        %Test presence/absence of a cycle-slip at the current epoch.
        %The state of the system is changed only for phase ambiguities
        if (length(phase) == 2)
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_kalman(pos_M, pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_kalman(pos_M, pos_R, X_t1_t(o3+33:o3+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
            [cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_kalman(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
            [cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_kalman(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+33:o3+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_observv(X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_observv(X_t1_t(o3+33:o3+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_sigma(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), diag(Cee(o3+1:o3+32,o3+1:o3+32)), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_sigma(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+33:o3+64), diag(Cee(o3+33:o3+64,o3+33:o3+64)), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time, 2);
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time-1, time, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+33:o3+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time-1, time, 2);

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
                if sqrt(sum((pos_R - X_t1_t([1,o1+1,o2+1])).^2)) > 3
                [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, pos_R, X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 10, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_observv(X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_sigma(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), diag(Cee(o3+1:o3+32,o3+1:o3+32)), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time-1, time, 1);
                else
                    [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
                end
            else
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, pos_R, X_t1_t(o3+1:o3+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
                [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_observv(X_t1_t(o3+1:o3+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_sigma(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), diag(Cee(o3+1:o3+32,o3+1:o3+32)), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time, 2);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip(pos_M, X_t1_t([1,o1+1,o2+1]), X_t1_t(o3+1:o3+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time-1, time, 2);
            end
            if (cycle_slip_found == 1)
                check_cs = 1;
                conf_cs(sat_slip) = 1;
                X_t1_t(o3+sat_slip) = N_slip;
                Cvv(o3+sat_slip,o3+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
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
    %positioning done only by the system dynamics

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';
end

%   vvX = Xhat_t_t(2);
%   vvY = Xhat_t_t(o1+2);
%   vvZ = Xhat_t_t(o2+2);
%   vvv = sqrt(vvX^2 + vvY^2 + vvZ^2);

%positioning error
%sigma_rho = sqrt(Cee(1,1) + Cee(o1+1,o1+1) + Cee(o2+1,o2+1));