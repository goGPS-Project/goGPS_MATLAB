%  CLASS LS_Manipulator_new
% =========================================================================
%
% DESCRIPTION
%   Manipulate least squares system in a new and more rational way w.r.t
%   the old class
%
% EXAMPLE
%   LSM = LS_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     Andrea Gatti
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef LS_Manipulator_new < handle
    
    % ==================================================================================================================================================
    %% Parameter columns id (order)
    % ==================================================================================================================================================
    
    properties (Constant)
        PAR_REC_X = 1;
        PAR_REC_Y = 2;
        PAR_REC_Z = 3;
        PAR_REC_EB = 4;
        PAR_AMB = 5;
        PAR_REC_CLK = 6;
        PAR_TROPO = 7;
        PAR_TROPO_N = 8;
        PAR_TROPO_E = 9;
        PAR_TROPO_V = 10;
        PAR_TROPO_Z = 24;
        PAR_SAT_CLK = 11;
        PAR_ANT_MP = 12;
        PAR_IONO = 13;
        PAR_TROPO_S = 14;
        PAR_SAT_X = 15;
        PAR_SAT_Y = 16;
        PAR_SAT_Z = 17;
        PAR_SAT_EB = 18;
        PAR_REC_EB_LIN = 19;
        PAR_REC_PPB = 20; % phase psudorange bias (to make ambiguity integer with common phase pseudorange obervables)
        PAR_SAT_PPB = 21; % phase psudorange bias (to make ambiguity integer with common phase pseudorange obervables)
        PAR_REC_EBFR = 22; % electrinic bias frequency depandant
        PAR_SAT_EBFR = 23; % electrinic bias frequency depandant
        PAR_REC_CLK_PH = 25;
        PAR_REC_CLK_PR = 26;
        PAR_SAT_CLK_PH = 27;
        PAR_SAT_CLK_PR = 28;
        
        
        CLASS_NAME = {'REC_X', 'REC_Y', 'REC_Z', 'REC_EB', 'AMB', 'REC_CLK', 'TROPO', 'TROPO_N', 'TROPO_E', 'TROPO_V', 'SAT_CLK', 'ANT_MP', 'IONO', 'TROPO_S', 'SAT_X', 'SAT_Y', 'SAT_Z', 'SAT_EB', 'REC_EB_LIN', 'REC_PPB', 'SAT_PPB', 'REC_EBFR', 'SAT_EBFR'};
    end
    
    properties
        %%% sparse A , this in more efficient than a sparse(slow to index in matlab) because
        %%% normally all the observation euation have the same numebr of
        %%% entry
        A % observations
        A_idx
        A_full 
        obs
        res
        param_class % class id of the column of A
        param_class_o % class id of the originallly requested paramaters
        
        time_obs   % epoch of the observation (GPS_Time)
        ref_time_obs % epoch of the observation since time mimimum (seconds)
        rate_obs     % rate of observations
        satellite_obs % goid satellite of the observations
        receiver_obs % receiver of the observations
        azimuth_obs % azimuth of the observations
        elevation_obs % elevation of the observations
        variance_obs % variance of the observations
        reweight_obs % rewight factor of observations
        obs_codes_id_obs % id of the signal used
        phase_obs % logical to tell if obs are phase or code
        wl_id_obs % id of the wavelength
        outlier_obs % is obs an outlier?
        
        A_pseudo
        A_idx_pseudo
        time_pseudo   % epoch of the pseudo-observation (GPS_Time)
        satellite_pseudo % goid satellite of the pseudo-observation
        receiver_pseudo % receiver of the pseudo-observation
        variance_pseudo % varaince of the pseudo-observation
        
        unique_obs_codes % uniques ids (cell) of the signals / since lot of combinations are possible they will be generated dynamically)
        unique_obs_codes_band  % band of the unique ob code
        unique_obs_codes_sys_c % system of the unique ob code
        unique_wl % set of unique wavelength
        rec_xyz % receiver coordinates to be used in
        unique_rec_name % names of the receivers
        unique_sat_goid % unique satellite goids
        cycle_slips = cell(1) % epoch of the cycle slip
        unique_time % unique epoch of the system
        
        time_par   % time of the parameter!!! very important the parameters (within the same class e.g. time for satellite s ) MUST be ordered in chronological order
        param_par  % parametrization of the parameter
        time_min   % ref_time_of the parameter
        obs_rate
        rec_par    % receiver of the parameter
        sat_par    % satellite of the parameter
        class_par  % class of the parameter
        obs_codes_id_par  % obs code id fo the parameter
        wl_id_par % wl id of the parameter
        out_par % paramters that are observed only by outlier observation
        phase_par % is pahse coode or both
        
        rec_set % set of receivers
        sat_set % set of satellites
        ch_set  % set of observation codes
        
        param_ch_set % ch set paramterization
        
        rec_amb_jmp % set of recievr ambiguity jmps
        sat_amb_jmp % set of satellite ambiguity jmps
        
        N
        idx_rd % idx parameter removed to silve the rank deficency
        ls_parametrization;
        
        free_tropo = false;% free network tropo paramters
        
        x
        
        coo_vcv; % variance covariance matrix of the coordinates
        
        log
    end
    
    methods
        function this = LS_Manipulator_new()
            this.log = Core.getLogger;
        end
    end
    
    methods
        function addObsEq(this, rec, obs_set, param_selection)
            % add observation equations to the matrices
            %
            % SYNTAX:
            %    this.addObsEq(rec, obs_set)
            n_par = length(param_selection);
            
            if isempty(this.param_class)
                this.param_class = param_selection;
            end
            if nargin < 4
                param_selection = this.param_class;
            end
            
            % --- check all parameters presence and their order -----
            par_rec_x_lid = param_selection == this.PAR_REC_X;
            par_rec_x = sum(par_rec_x_lid) > 0;
            par_rec_y_lid = param_selection == this.PAR_REC_Y;
            par_rec_y = sum(par_rec_y_lid) > 0;
            par_rec_z_lid = param_selection == this.PAR_REC_Z;
            par_rec_z = sum(par_rec_z_lid) > 0;
            
            par_sat_x_lid = param_selection == this.PAR_SAT_X;
            par_sat_x = sum(par_sat_x_lid) > 0;
            par_sat_y_lid = param_selection == this.PAR_SAT_Y;
            par_sat_y = sum(par_sat_y_lid) > 0;
            par_sat_z_lid = param_selection == this.PAR_SAT_Z;
            par_sat_z = sum(par_sat_z_lid) > 0;
            
            par_rec_eb_lid = param_selection == this.PAR_REC_EB;
            par_rec_eb = sum(par_rec_eb_lid) > 0;
            par_rec_ppb_lid = param_selection == this.PAR_REC_PPB;
            par_rec_ppb = sum(par_rec_ppb_lid) > 0;
            par_rec_ebfr_lid = param_selection == this.PAR_REC_EBFR;
            par_rec_ebfr = sum(par_rec_ebfr_lid) > 0;
            par_rec_eb_lin_lid = param_selection == this.PAR_REC_EB_LIN;
            par_rec_eb_lin = sum(par_rec_eb_lin_lid) > 0;
            par_sat_eb_lid = param_selection == this.PAR_SAT_EB;
            par_sat_eb = sum(par_sat_eb_lid) > 0;
            par_sat_ppb_lid = param_selection == this.PAR_SAT_PPB;
            par_sat_ppb = sum(par_sat_ppb_lid) > 0;
            par_sat_ebfr_lid = param_selection == this.PAR_SAT_EBFR;
            par_sat_ebfr = sum(par_sat_ebfr_lid) > 0;
            
            par_amb_lid = param_selection == this.PAR_AMB;
            par_amb = sum(par_amb_lid) > 0;
            
            par_rec_clk_lid = param_selection == this.PAR_REC_CLK;
            par_rec_clk = sum(par_rec_clk_lid) > 0;
            
            par_sat_clk_lid = param_selection == this.PAR_SAT_CLK;
            par_sat_clk = sum(par_sat_clk_lid) > 0;
            
            par_tropo_lid = param_selection == this.PAR_TROPO;
            par_tropo = sum(par_tropo_lid) > 0;
            
            par_tropo_e_lid = param_selection == this.PAR_TROPO_E;
            par_tropo_e = sum(par_tropo_e_lid) > 0;
            
            par_tropo_n_lid = param_selection == this.PAR_TROPO_N;
            par_tropo_n = sum(par_tropo_n_lid) > 0;
            
            par_tropo_v_lid = param_selection == this.PAR_TROPO_V;
            par_tropo_v = sum(par_tropo_v_lid) > 0;
            
            par_tropo_s_lid = param_selection == this.PAR_TROPO_S;
            par_tropo_s = sum(par_tropo_s_lid) > 0;
            
            par_tropo_z_lid = param_selection == this.PAR_TROPO_Z;
            par_tropo_z = sum(par_tropo_z_lid) > 0;
            
            par_iono_lid = param_selection == this.PAR_IONO;
            par_iono = sum(par_iono_lid) > 0;
            
            par_rec_clk_pr_lid = param_selection == this.PAR_REC_CLK_PR;
            par_rec_clk_pr = sum(par_rec_clk_pr_lid) > 0;
            
            par_rec_clk_ph_lid = param_selection == this.PAR_REC_CLK_PH;
            par_rec_clk_ph = sum(par_rec_clk_ph_lid) > 0;
            
            par_sat_clk_pr_lid = param_selection == this.PAR_SAT_CLK_PR;
            par_sat_clk_pr = sum(par_sat_clk_pr_lid) > 0;
            
            par_sat_clk_ph_lid = param_selection == this.PAR_SAT_CLK_PH;
            par_sat_clk_ph = sum(par_sat_clk_ph_lid) > 0;
            
            % ---- add the recever to the recibers
            if Core_Utils.findAinB(rec.parent.getMarkerName,this.unique_rec_name) == 0
                this.unique_rec_name{end+1} = rec.parent.getMarkerName;
                this.rec_xyz = [this.rec_xyz; rec.getMedianPosXYZ];
                r = size(this.rec_xyz,1);
            else
                r = Core_Utils.findAinB(rec.parent.getMarkerName,this.unique_rec_name);
            end
            
            % --- remove syntetic measurement form the observations -----
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            xs_loc = zero2nan(xs_loc);
            diff_obs = nan2zero(zero2nan(obs_set.obs) - zero2nan(synt_obs));
            n_obs = sum(sum(diff_obs ~= 0));
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            n_stream = size(diff_obs, 2); % number of satellites present in the observation set
            
            
            
            % add signal to the unique
            u_obs_code = unique(cellstr(obs_set.obs_code));
            for i = 1 : length(u_obs_code)
                if Core_Utils.findAinB(u_obs_code(i),this.unique_obs_codes) == 0
                    this.unique_obs_codes{end+1} = u_obs_code{i};
                end
            end
            obs_code_id = Core_Utils.findAinB(cellstr(obs_set.obs_code),this.unique_obs_codes);
            
            
            % add wavelength to the unique
            u_wl = unique(round(obs_set.wl*1e15))/1e15;
            this.unique_wl = unique([this.unique_wl; u_wl]*1e15)/1e15;
            [~, obs_wl_id] = ismember(round(obs_set.wl*1e7), round(this.unique_wl*1e7));
            
            
            mfw_on =  sum(param_selection == this.PAR_TROPO | param_selection == this.PAR_TROPO_E | param_selection == this.PAR_TROPO_N | param_selection == this.PAR_TROPO_V | param_selection == this.PAR_TROPO_Z) > 0;
            % get the mapping function for tropo
            if mfw_on
                id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
                [~, mfw] = rec.getSlantMF(id_sync_out);
               % mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
            % check whivh observations are phase ones
            phase_s = obs_set.obs_code(:,2) == 'L'; % first char is the system second MUST be the bervation type (Phase pseudorance doppler snr)
            
            % initliaze the matrices
            A = zeros(n_obs, n_par);
            [obs,satellite_obs, azimuth_obs, elevation_obs, variance_obs, wl_obs] = deal(zeros(n_obs, 1));
            [ obs_codes_id_obs,  wl_obs] = deal(zeros(n_obs, 1,'uint8'));
            [ satellite_obs] = deal(zeros(n_obs, 1,'uint16'));
            
            phase_obs = false(n_obs,1);
            time_obs = GPS_Time();
            obs_count = 1;
            if length(this.cycle_slips) < r || isempty(this.cycle_slips{r})
                this.cycle_slips{r} = {};
                this.cycle_slips{r} = cell(max([this.unique_sat_goid obs_set.go_id']),1);
                for s = 1 : max([this.unique_sat_goid obs_set.go_id'])
                    this.cycle_slips{r}{s} = cell(length(this.unique_obs_codes),1);
                end
            end
            iono_const = 40.3*10^16;%GPS_SS.L_VEC(1)^2;
            % fill the A matrix per satellite
            for s = 1 : n_stream
                id_ok_stream = diff_obs(:, s) ~= 0; % check observation existence -> logical array for a "s" stream
                if any(id_ok_stream)
                    %---- some quantities related to the stream
                    obs_stream = diff_obs(id_ok_stream, s);
                    n_obs_stream = length(obs_stream);
                    lines_stream = obs_count + (0:(n_obs_stream - 1));
                    
                    el_stream = obs_set.el(id_ok_stream, s) / 180 * pi;
                    az_stream = obs_set.az(id_ok_stream, s) / 180 * pi;
                    s_go_id   = obs_set.go_id(s);
                    s_s_id    = obs_code_id(s);
                    wl_id     = obs_wl_id(s);
                    if mfw_on
                        mfw_stream = mfw(id_ok_stream, s_go_id);
                    end
                    this.unique_sat_goid = unique([this.unique_sat_goid  s_go_id]);
                    
                    
                    xs_loc_stream = permute(xs_loc(id_ok_stream, s, :), [1, 3, 2]);
                    los_stream = rowNormalize(xs_loc_stream);
                    
                    %--- Fill Observation related vectors------------
                    obs(lines_stream) = obs_stream;
                    satellite_obs(lines_stream) = s_go_id;
                    variance_obs(lines_stream) =  obs_set.sigma(s)^2;
                    azimuth_obs(lines_stream) =  obs_set.az(id_ok_stream,s);
                    elevation_obs(lines_stream) =  obs_set.el(id_ok_stream,s);
                    obs_codes_id_obs(lines_stream) = s_s_id;
                    phase_obs(lines_stream) = phase_s(s);
                    wl_obs(lines_stream) = wl_id;
                    time_obs.addEpoch(obs_set.time.getEpoch(id_ok_stream).getMatlabTime);
                    
                    % ----------- save the cycle slips
                    if phase_s(s)
                        %obs_set.cycle_slip(,s) = true;
                        if any(obs_set.obs(:,s))
                            obs_set.cycle_slip(find(obs_set.obs(:,s) ~= 0, 1, 'first'), s) = true;
                        end
                        this.cycle_slips{r}{s_go_id}{s_s_id} = obs_set.time.getEpoch(find(obs_set.cycle_slip(:, s) & obs_set.obs(:, s) ~= 0));
                    end
                    
                    % ----------- FILL IMAGE MATRIX ------------
                    % ----------- coordinates ------------------
                    if par_rec_x
                        A(lines_stream, par_rec_x_lid) = - los_stream(:,1);
                    end
                    
                    if par_rec_y
                        A(lines_stream, par_rec_y_lid) = - los_stream(:,2);
                    end
                    
                    if par_rec_z
                        A(lines_stream, par_rec_z_lid) = - los_stream(:,3);
                    end
                    
                    if par_sat_x
                        A(lines_stream, par_sat_x_lid) =  los_stream(:,1);
                    end
                    
                    if par_sat_y
                        A(lines_stream, par_sat_y_lid) =  los_stream(:,2);
                    end
                    
                    if par_sat_z
                        A(lines_stream, par_sat_z_lid) =  los_stream(:,3);
                    end
                    % ----------- electronic bias ------------------
                    if par_rec_eb
                        A(lines_stream, par_rec_eb_lid) = 1;
                    end
                    if par_rec_ebfr
                        A(lines_stream, par_rec_ebfr_lid) = 1;
                    end
                    if par_rec_ppb && phase_s(s)
                        A(lines_stream, par_rec_ppb_lid) = 1;
                    end
                    if par_rec_eb_lin
                        A(lines_stream, par_rec_eb_lin_lid) = 1/obs_set.wl(s);
                    end
                    if par_sat_eb
                        A(lines_stream, par_sat_eb_lid) = 1;
                    end
                    if par_sat_ebfr
                        A(lines_stream, par_sat_ebfr_lid) = 1;
                    end
                    if par_sat_ppb && phase_s(s)
                        A(lines_stream, par_sat_ppb_lid) = 1;
                    end
                    % ----------- Abiguity ------------------
                    if par_amb && phase_s(s)
                        A(lines_stream, par_amb_lid) = obs_set.wl(s);
                    end
                    % ----------- Clock ------------------
                    if par_rec_clk
                        A(lines_stream, par_rec_clk_lid) = 1;
                    end
                    if par_sat_clk
                        A(lines_stream, par_sat_clk_lid) = 1;
                    end
                    if par_rec_clk_pr  && ~phase_s(s)
                        A(lines_stream, par_rec_clk_pr_lid) = 1;
                    end
                    if par_sat_clk_pr  && ~phase_s(s)
                        A(lines_stream, par_sat_clk_pr_lid) = 1;
                    end
                    if par_rec_clk_ph  && phase_s(s)
                        A(lines_stream, par_rec_clk_ph_lid) = 1;
                    end
                    if par_sat_clk_ph  && phase_s(s)
                        A(lines_stream, par_sat_clk_ph_lid) = 1;
                    end
                    % ----------- ZTD ------------------
                    if par_tropo
                        A(lines_stream, par_tropo_lid) = mfw_stream;
                    end
                    % ----------- ZTD gradients ------------------
                    if par_tropo_n || par_tropo_e
                        cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                        if par_tropo_e
                            A(lines_stream, par_tropo_e_lid) = sin(az_stream) .* cotan_term /10; % east gradient
                        end
                        if par_tropo_n
                            A(lines_stream, par_tropo_n_lid) = cos(az_stream) .* cotan_term /10; % noth gradient
                        end
                    end
                    if par_tropo_v
                        A(lines_stream, par_tropo_v_lid) = mfw_stream*rec.h_ellips;
                    end
                    if par_tropo_z
                        n_pol = sum(par_tropo_z_lid);
                        degree = ceil(-3/2 + sqrt(9/4 + 2*(n_pol -1)));
                        A(lines_stream, par_tropo_z_lid) = repmat(mfw_stream,1,n_pol).*Core_Utils.getAllZernike(degree, az_stream, el_stream);
                    end
                    % ----------- Ionosphere delay --------------------
                    if par_iono
                        if phase_s(s)
                            A(lines_stream, par_iono_lid) = - iono_const*(obs_set.wl(s)/Core_Utils.V_LIGHT).^2; %obs_set.wl(s)^2/iono_const;
                        else
                            A(lines_stream, par_iono_lid) =  iono_const*(obs_set.wl(s)/Core_Utils.V_LIGHT).^2; obs_set.wl(s)^2/iono_const;
                        end
                    end
                    obs_count = obs_count + n_obs_stream;
                end
            end
            this.A = [this.A; A];
            this.obs = [this.obs; obs];
            this.satellite_obs = [this.satellite_obs; satellite_obs];
            this.azimuth_obs = [this.azimuth_obs; azimuth_obs];
            this.elevation_obs = [this.elevation_obs; elevation_obs];
            this.obs_codes_id_obs = [this.obs_codes_id_obs; obs_codes_id_obs];
            this.variance_obs = [this.variance_obs; variance_obs];
            this.phase_obs = [this.phase_obs; phase_obs];
            this.wl_id_obs = [this.wl_id_obs; wl_obs];
            
            if isempty(this.time_obs)
                this.time_obs = time_obs;
            else
                this.time_obs.addEpoch(time_obs.getMatlabTime);
            end
            this.receiver_obs = [this.receiver_obs; r*ones(size(phase_obs))];
        end
        
        function bondParamsGenerateIdx(this, ls_parametrization)
            % bond parameters (splines or other models) and generate idx
            %
            % SYNTAX
            %    this.bondParamGenerateIdx(parametrization)
            
            % generate jmps idx in for both receiver and satellites
            this.log = Core.getLogger();
            this.ls_parametrization = ls_parametrization;
            this.computeRefTimeObs();
            this.computeAmbJmps();
            
            
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            n_obs = size(this.A,1);
            this.A_idx = zeros(size(this.A),'uint32');
            
          
            time_obs = round(this.ref_time_obs/this.obs_rate);
            obs_rate = this.obs_rate;
            time_min = this.time_min.getMatlabTime;
            
            % generate system and band of the uniqes obs codes
            this.unique_obs_codes_band = char(zeros(size(this.unique_obs_codes)));
            this.unique_obs_codes_sys_c = char(zeros(size(this.unique_obs_codes)));
            for o = 1 : length(this.unique_obs_codes)
                this.unique_obs_codes_sys_c(o) = this.unique_obs_codes{o}(1);
                this.unique_obs_codes_band(o) = this.unique_obs_codes{o}(3);
                
            end
            
            cumulative_idx = 0;
            i_col = 1;
            
            % is the observatio  pahse or code
            phpr_unique_obs_codes = zeros(size(this.unique_obs_codes));
            for o = 1 : length(this.unique_obs_codes)
                if this.unique_obs_codes{o}(2) == 'C'
                    phpr_unique_obs_codes(o) = 1;
                elseif this.unique_obs_codes{o}(2) == 'L'
                    phpr_unique_obs_codes(o) = 2;
                end
            end
            % ---- preallocation to speed up
            ch_lid = false(size(this.A,1),1);
            obs_lid = false(size(this.A,1),1);
            
            
            
            this.time_par = zeros(n_obs,2,'uint32');
            this.param_par = zeros(n_obs,4,'uint8'); % time of the parameter
            this.rec_par = zeros(n_obs,1,'int16'); % receiver of the parameters
            this.sat_par = zeros(n_obs,1,'int8');  % receiver of the parameters
            this.class_par = zeros(n_obs,1,'uint8');  % class of the parameter
            this.obs_codes_id_par= zeros(n_obs,1,'int8');  % obs_code id parameters
            this.wl_id_par= zeros(n_obs,1,'uint8');  % obs_code id parameters
            this.phase_par= zeros(n_obs,1,'uint8');  % phse code or both parameter
            
            this.outlier_obs = false(size(this.obs));
            this.param_class_o = this.param_class;
            ch_set_old = []; % avoid repating expensive task
            i_p = 1;
            i_p_o = 1;
            while i_p <= length(this.param_class)
                has_not_expanded = true;
                col_incr = 0;
                p = this.param_class(i_p);
                [parametriz, opt] = ls_parametrization.getParametrization(p);
                % ----------------- defining receiver sets ---------
                if parametriz(2) == ls_parametrization.SING_REC
                    n_rec_set = size(this.rec_xyz,1);
                    rec_set = num2cell(1:n_rec_set);
                elseif parametriz(2) == ls_parametrization.ALL_REC
                    n_rec_set = 1;
                    rec_set = {1:size(this.rec_xyz,1)};
                elseif parametriz(2) == ls_parametrization.MULTI_REC
                    n_rec_set = length(opt.rec_sets);
                    rec_set = opt.rec_sets;
                end
                % ----------------- defining satellites sets ---------
                if parametriz(3) == ls_parametrization.SING_SAT
                    n_sat_set = length(this.unique_sat_goid);
                    sat_set = num2cell(this.unique_sat_goid);
                elseif parametriz(3) == ls_parametrization.ALL_SAT
                    n_sat_set = 1;
                    sat_set = {this.unique_sat_goid};
                elseif parametriz(3) == ls_parametrization.MULTI_SAT
                    n_sat_set = length(opt.sat_sets);
                    sat_set = opt.sat_sets;
                end
                
                % ----------------- defining channel sets ---------
                % construct channel to frequency mapping
                sig2wl = zeros(size(this.unique_obs_codes));
                sig_p_id = 1:numel(this.unique_obs_codes);
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        sig2wl(c) = this.wl_id_obs(idx1);
                    end
                end
                % construct channel to phase mapping
                sig2phase = zeros(size(this.unique_obs_codes));
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        sig2phase(c) = this.phase_obs(idx1);
                    end
                end
                % construct channel to constellation mapping
                sig2const = char(zeros(size(this.unique_obs_codes)));
                sig2const_num = zeros(size(this.unique_obs_codes));
                
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        ant_id = Core.getConstellationCollector.getAntennaId(this.satellite_obs(idx1));
                        sig2const(c) = ant_id(1);
                        sig2const_num(c) = find(Core.getConstellationCollector.SYS_C == ant_id(1));
                    end
                end
                if parametriz(4) == ls_parametrization.SING_TRACK
                    n_ch_set = numel(this.unique_obs_codes);
                    ch_set = num2cell(uint8(1:n_ch_set));
                elseif parametriz(4) == ls_parametrization.ALL_FREQ
                    n_ch_set = 1;
                    ch_set = {uint8(1:numel(this.unique_obs_codes))};
                elseif parametriz(4) == ls_parametrization.SING_FREQ
                    n_ch_set = length(this.unique_wl);
                    ch_set = {};
                    for c = 1 : n_ch_set
                        ch_set{c} = uint8(sig_p_id(sig2wl == c));
                    end
                elseif parametriz(4) == ls_parametrization.SING_BAND
                    u_bnd = unique(Core_Utils.code2Char2Num([this.unique_obs_codes_sys_c' this.unique_obs_codes_band']));
                    n_ch_set = length(u_bnd);
                    ch_set = {};
                    for c = 1 : n_ch_set
                        ch_set{c} = uint8(sig_p_id(Core_Utils.code2Char2Num([this.unique_obs_codes_sys_c' this.unique_obs_codes_band']) == u_bnd(c)));
                    end
                elseif parametriz(4) == ls_parametrization.FREQ_CONST
                    u_wl_const = unique(sig2wl*100+sig2const_num)
                    n_ch_set = length(u_wl_const);
                    ch_set = {};
                    for c = 1 : n_ch_set
                        const = rem(u_wl_const(c),100);
                        wl = floor(u_wl_const(c)/100);
                        ch_set{c} = uint8(sig_p_id(sig2const_num == const & sig2wl == wl));
                    end
                elseif parametriz(4) == ls_parametrization.SING_FREQ_BIN
                    n_ch_set = length(this.unique_wl);
                    ch_set = {};
                    % waveleght are binned into theese bands
                    fr_bin = [1176.450*1e6+[-1e5  +1e5];  %  L5a   -> G5 , I5 , J5, E5 , S5
                        1191.795*1e6+[-1e5  +1e5];  %  L5ab  -> E8
                        1202.025*1e6+[-1e5  +1e5];  %  L5c   -> R3
                        1207.140*1e6+[-1e5  +1e5];  %  L5B   -> C7 , E7,
                        1227.600*1e6+[-1e5  +1e5];  %  L2a   -> G2 , J2
                        1246.000*1e6+[-7*7/16*1e6-1e5  +12*7/16*1e6+1e5];  %  L2b   -> R2
                        1268.520*1e6+[-1e5  +1e5];  %  L6a   -> C3
                        1278.750*1e6+[-1e5  +1e5];  %  L6b   -> J6 , E6
                        1561.098*1e6+[-1e5  +1e5];  %  L1a   -> C2
                        1575.42*1e6+[-1e5  +1e5];  %  L1b   -> G1, E1, J1, S1
                        1602.000*1e6+[-7*9/16*1e6-1e5  +12*9/16*1e6+1e5];  %  L1c   -> R1
                        2492.028*1e6+[-1e5  +1e5];  %  L9   -> I9
                        ];
                    wl_bin = fliplr(Core_Utils.V_LIGHT./fr_bin);
                    wl2bnd = zeros(size(this.unique_wl));
                    wl2chid = zeros(size(this.unique_wl));
                    
                    for w = 1 : length(this.unique_wl)
                        b = 1;
                        while wl2bnd(w) == 0
                            if this.unique_wl(w) > wl_bin(b,1) && this.unique_wl(w) < wl_bin(b,2)
                                wl2bnd(w) = b;
                            else
                                b = b+1;
                            end
                        end
                        wl2chid(w) = this.obs_codes_id_obs(find(this.wl_id_obs == w,1,'first'));
                    end
                    u_bnd = unique(wl2bnd);
                    n_ch_set = length(u_bnd);
                    for c = 1 : n_ch_set
                        chids = wl2chid(wl2bnd == u_bnd(c));
                        ch_set{c} = chids;
                    end
                elseif parametriz(4) == ls_parametrization.RULE
                    % ---- evaluate the rule
                    n_ch_set = 0;
                    ch_set = {};
                    for rule = opt.rule
                        cond_sig = true(size(sig_p_id));
                        parts = strsplit(rule{1},':');
                        cond = parts{1};
                        optt = parts{2};
                        conds = strsplit(cond,'&');
                        for c = conds
                            if strcmpi('PSRANGE',c)
                                cond_sig = cond_sig & ~sig2phase;
                            elseif strcmpi('PHASE',c)
                                cond_sig = cond_sig & sig2phase;
                            elseif strcmpi('GLONASS',c)
                                cond_sig = cond_sig & sig2const == 'R';
                            elseif strcmpi('NOT*GLONASS',c)
                                cond_sig = cond_sig & sig2const ~= 'R';
                            else
                                cond_sig(:) = false;
                                this.log.addWarning(sprintf('Unknown option %s in parametrization rule, ignoring parameter',c));
                            end
                        end
                        if any(cond_sig)
                            if str2num(optt) == ls_parametrization.ALL_FREQ
                                n_ch_set = n_ch_set +1;
                                ch_set{n_ch_set} = uint8(sig_p_id(cond_sig));
                            elseif str2num(optt) == ls_parametrization.SING_FREQ
                                u_wl_tmp = unique(sig2wl(cond_sig));
                                n_wl_tmp = length(u_wl_tmp);
                                for c = 1 : n_wl_tmp
                                    n_ch_set = n_ch_set +1;
                                    ch_set{n_ch_set} = uint8(sig_p_id(sig2wl == u_wl_tmp(c) & cond_sig));
                                end
                            elseif str2num(optt) == ls_parametrization.SING_TRACK
                                n_ch_set = n_ch_set +length(sig_p_id(cond_sig));
                                ch_set = [ch_set num2cell(uint8(sig_p_id(cond_sig)))];
                            elseif str2num(optt) == ls_parametrization.FREQ_CONST
                                u_wl_const = unique(sig2wl(cond_sig)*100+sig2const_num(cond_sig));
                                for c = 1 : length(u_wl_const)
                                    n_ch_set = n_ch_set +1;
                                    const = rem(u_wl_const(c),100);
                                    wl = floor(u_wl_const(c)/100);
                                    ch_set{n_ch_set} = uint8(sig_p_id(sig2const_num == const & sig2wl == wl & cond_sig));
                                end
                            end
                            
                        end
                        
                    end
                end
                this.param_ch_set{i_p_o} = ch_set;
                
                % ------- Now constructing the index, the epoch dependednce will be dealt inside the loop -------------
                for r = 1 : n_rec_set
                    rec_lid = false(size(this.A,1),1);
                    for rr = rec_set{r}
                        rec_lid = rec_lid | this.receiver_obs == rr;
                    end
                    % find an id for the receiver set to keep track of the
                    % parameters if recievr is sigle is simply the receiver
                    % progessince number, otherwise they are negative
                    % number with the index in the receiver set
                    if length(rec_set{r}) == 1
                        r_id = rec_set{r};
                    else
                        r_id = -Core_Utils.findAinB(rec_set{r},this.rec_set);
                        if r_id == 0
                            r_id = -length(this.rec_set) -1;
                            this.rec_set{end+1} = rec_set{r};
                        end
                        
                    end
                    for s = 1 : n_sat_set
                        %                         sat_lid = false(size(this.A,1),1);
                        %                         for ss = sat_set{s}
                        %                             sat_lid = sat_lid | this.satellite_obs == ss;
                        %                         end
                        sat_lid = ismember(this.satellite_obs,sat_set{s});
                        ss = sat_set{s};
                        if length(sat_set{s}) == 1
                            s_id = sat_set{s};
                        else
                            s_id = -Core_Utils.findAinB(sat_set{s},this.sat_set);
                            if s_id == 0
                                s_id = -length(this.sat_set) -1;
                                this.sat_set{end+1} = sat_set{s};
                            end
                            
                        end
                        for f = 1 : n_ch_set
                            if p ~= this.PAR_AMB || this.unique_obs_codes{ch_set{f}}(2) == 'L'
                                %                                 ch_lid(:) = false;
                                %                                 for cc = ch_set{f}
                                %                                     ch_lid = ch_lid | this.obs_codes_id_obs == cc;
                                %                                 end
                                if ~isequal(ch_set_old,ch_set{f})
                                    %ch_lid = ismemberBuiltinTypes(this.obs_codes_id_obs,ch_set{f});
                                    ch_lid = ismember(this.obs_codes_id_obs, ch_set{f});
                                    ch_set_old = ch_set{f};
                                end
                                cc = ch_set{f};
                                if length(ch_set{f}) == 1
                                    ch_id = ch_set{f};
                                else
                                    ch_id = -Core_Utils.findAinB(ch_set{f},this.ch_set);
                                    if ch_id == 0
                                        ch_id = -length(this.ch_set) -1;
                                        this.ch_set{end+1} = ch_set{f};
                                    end
                                end
                                % is a phase pr or both paramter
                                phase_code = 0;
                                u_phpr = unique(phpr_unique_obs_codes(cc));
                                if length(u_phpr) == 1
                                    phase_code = u_phpr;
                                end
                                %---- final observation id
                                obs_lid = rec_lid & sat_lid & ch_lid;
                                if any(obs_lid)
                                    % --- now dealing with epoch dependence ----
                                    cols_tmp = 0;
                                    if parametriz(1) == ls_parametrization.CONST
                                        ep_pgr_id = 1;
                                        n_prg_id = 1;
                                        time_par_tmp = [min(this.ref_time_obs(obs_lid))  max(this.ref_time_obs(obs_lid))];
                                    elseif parametriz(1) == ls_parametrization.EP_WISE
                                        ep_id = time_obs(obs_lid);
                                        u_e_tmp = unique(ep_id);
                                        [~,ep_pgr_id] = ismember(ep_id,u_e_tmp);
                                        n_prg_id = length(u_e_tmp);
                                        time_par_tmp = [u_e_tmp*obs_rate zeros(size(u_e_tmp))];
                                    elseif parametriz(1) == ls_parametrization.STEP_CONST
                                        ep_id = time_obs(obs_lid);
                                        ep_pgr_id = zeros(size(ep_id));
                                        n_prg_id = 0;
                                        time_par_tmp = [];
                                        if p == this.PAR_AMB %% ambiguity case is really unique and make sense to treat it separatly
                                            p_s = 1;
                                            if length(this.cycle_slips{rr}{ss}) >= cc
                                                cs = this.cycle_slips{rr}{ss}{cc};
                                                if ~isempty(cs)
                                                    steps = round(cs.getRefTime(time_min)/obs_rate);
                                                    time_par_tmp = zeros(length(steps),2);
                                                    for st = steps'
                                                        lid_maj = ep_id >= st;
                                                        ep_pgr_id(lid_maj) = p_s;
                                                        time_par_tmp(p_s,:) = [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]; %start of the arc
                                                        if p_s > 1
                                                            last_ep = ep_id(find(~lid_maj,1,'last'));
                                                            if ~isempty(last_ep)
                                                                time_par_tmp(p_s-1,2) = last_ep*obs_rate; % end of the arc
                                                            else
                                                                time_par_tmp(p_s-1,2) = ep_id(end)*obs_rate; % end of the arc
                                                            end
                                                        end
                                                        p_s = p_s +1;
                                                    end
                                                    n_prg_id = p_s-1;
                                                else
                                                    ep_pgr_id = 1;
                                                    n_prg_id = 1;
                                                    time_par_tmp = [ep_id(1)*obs_rate ep_id(end)*obs_rate];
                                                    
                                                end
                                            end
                                        else
                                            steps_set = opt.steps_set;
                                            
                                            if parametriz(2) == ls_parametrization.ALL_REC % you can use differents step for step wise satellite dependent paraters
                                                steps = round(steps_set{ss}.getNominalTime(obs_rate).getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    if any(lid_maj)
                                                        ep_pgr_id(lid_maj) = p_s;
                                                        time_par_tmp = [time_par_tmp; [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]]; %start of the arc
                                                        if p_s > 1
                                                            time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                        end
                                                        p_s = p_s +1;
                                                    end
                                                end
                                            elseif parametriz(3) == ls_parametrization.ALL_SAT  % you can use differents step for step wise receiver dependent paraters
                                                steps = round(steps_set{rr}.getNominalTime(obs_rate).getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    if any(lid_maj)
                                                        ep_pgr_id(lid_maj) = p_s;
                                                        time_par_tmp = [time_par_tmp; [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]]; %start of the arc
                                                        if p_s > 1
                                                            time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                        end
                                                        p_s = p_s +1;
                                                    end
                                                end
                                            end
                                            n_prg_id = p_s - 1;
                                            
                                        end
                                    elseif parametriz(1) == ls_parametrization.SPLINE_ZERO
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        u_e_tmp = unique(ep_id);
                                        [~,ep_pgr_id] = ismember(ep_id,u_e_tmp);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        n_prg_id = length(u_e_tmp);
                                    elseif parametriz(1) == ls_parametrization.SPLINE_LIN
                                        % ---- colum will be doubled -----
                                        cols_tmp = [ 0 1];
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate)/opt.spline_rate,1);
                                        u_e_tmp = unique([ep_id ep_id+1]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        a_col = this.A(obs_lid,i_col);
                                        if has_not_expanded
                                            this.A = [this.A(:,1:(i_col)) zeros(n_obs,1) this.A(:,(i_col+1):end)];
                                            this.A_idx = [this.A_idx(:,1:(i_col)) zeros(n_obs,1,'uint32') this.A_idx(:,(i_col+1):end)];
                                            this.param_class = [this.param_class(1:(i_col-1)); p*ones(2,1); this.param_class((i_col+1):end)];
                                            i_p = i_p + 1;
                                            has_not_expanded = false;
                                        end
                                        this.A = [this.A(:,1:(i_col-1)) zeros(n_obs,2) this.A(:,(i_col+1):end)];
                                        this.A(obs_lid,i_col + cols_tmp) = [a_col.*spline_v(:,1) a_col.*spline_v(:,2)];
                                        n_prg_id = length(u_e_tmp);
                                    elseif parametriz(1) == ls_parametrization.SPLINE_CUB
                                        % ---- colum will be quadrupled -----
                                        cols_tmp = [ 0 1 2 3];
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate)/opt.spline_rate,3);
                                        u_e_tmp = unique([ep_id ep_id+1 ep_id+2 ep_id+3]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(:,i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        
                                        if has_not_expanded
                                            this.A = [this.A(:,1:(i_col)) zeros(n_obs,3) this.A(:,(i_col+1):end)];
                                            this.A_idx = [this.A_idx(:,1:(i_col)) zeros(n_obs,3,'uint32') this.A_idx(:,(i_col+1):end)];
                                            this.param_class = [this.param_class(1:(i_col-1)); p*ones(4,1); this.param_class((i_col+1):end)];
                                            i_p = i_p + 3;
                                            has_not_expanded = false;
                                        end
                                        a_col = this.A(obs_lid,i_col);
                                        this.A(obs_lid,i_col + cols_tmp) = [a_col.*spline_v(:,1) a_col.*spline_v(:,2) a_col.*spline_v(:,3) a_col.*spline_v(:,4)];
                                        n_prg_id = length(u_e_tmp);
                                    end
                                    this.A_idx(obs_lid, i_col + cols_tmp) = cumulative_idx + ep_pgr_id;
                                    %[u_new_par] = unique(cumulative_idx + ep_pgr_id(:));
                                    
                                    this.time_par(cumulative_idx+(1:n_prg_id),:) =  uint32(time_par_tmp);% time of the parameter
                                    this.param_par(cumulative_idx+(1:n_prg_id),:) = repmat(uint8(parametriz),n_prg_id,1);% time of the parameter
                                    this.rec_par(cumulative_idx+(1:n_prg_id)) = r_id*ones(n_prg_id,1,'int8');  % receiver of the parameters
                                    this.sat_par(cumulative_idx+(1:n_prg_id)) = s_id*ones(n_prg_id,1,'int8');  % receiver of the parameters
                                    this.class_par(cumulative_idx+(1:n_prg_id)) =  p*ones(n_prg_id,1,'uint8');  % class of the parameter
                                    this.obs_codes_id_par(cumulative_idx+(1:n_prg_id)) = int8(ch_id)*ones(n_prg_id,1,'int8');  % obs_code id parameters
                                    this.phase_par(cumulative_idx+(1:n_prg_id)) = uint8(phase_code)*ones(n_prg_id,1,'uint8');  % phase code or both id parameters
                                    
                                    if ch_id > 0
                                        this.wl_id_par(cumulative_idx+(1:n_prg_id)) = uint8(sig2wl(ch_id))*ones(n_prg_id,1,'uint8');  % wavelength id parameters
                                    else
                                        u_wl = unique(sig2wl(this.ch_set{-ch_id})); % check if at least they belong to the same frequency
                                        if length(u_wl) == 1
                                            this.wl_id_par(cumulative_idx+(1:n_prg_id)) = uint8(u_wl)*ones(n_prg_id,1,'uint8');  % wavelength id parameters
                                        else
                                            this.wl_id_par(cumulative_idx+(1:n_prg_id)) = zeros(n_prg_id,1,'uint8');  % wavelength id parameters
                                        end
                                    end
                                    cumulative_idx = cumulative_idx + n_prg_id;
                                    col_incr = max(col_incr,length(cols_tmp));
                                end
                            end
                        end
                    end
                end
                i_col = i_col + col_incr;
                i_p = i_p + 1;
                i_p_o = i_p_o +1;
            end
            % free allocated space in excess
            this.time_par((cumulative_idx+1):end,:) =  [];%
            this.param_par((cumulative_idx+1):end,:) =  [];%
            this.rec_par((cumulative_idx+1):end) =  [];%
            this.sat_par((cumulative_idx+1):end) =  [];%
            this.class_par((cumulative_idx+1):end) =  [];%
            this.obs_codes_id_par((cumulative_idx+1):end) =  [];%
            this.wl_id_par((cumulative_idx+1):end) =  [];%
            this.phase_par((cumulative_idx+1):end) =  [];%
            this.out_par =  false(size(this.wl_id_par));%
            
            %
            %             for this,
            %             for o = 1 : length(this.unique_obs_codes)
            %                 this.unique_obs_codes_sys_c(o) = this.unique_obs_codes{o}(1);
            %                 this.unique_obs_codes_band(o) = this.unique_obs_codes{o}(3);
            %
            %             end
            %
            %
            %
        end
        
        function generateOutliedParIdx(this)
            % considering outlier generate an index of outlied paramters (i.e paramters that are observed only by outlier observation)
            %
            % SYNTAX
            %   this.generateOutliedParIdx();
            if isempty(this.outlier_obs)
                this.outlier_obs = false(size(this.obs));
            end
            out_par = unique(serialize(this.A_idx(this.outlier_obs > 0,:)));
            not_out_par = unique([serialize(this.A_idx(this.outlier_obs == 0,:)); this.A_idx_pseudo(:)]);
            this.out_par = Core_Utils.ordinal2logical(noZero(setdiff(out_par, not_out_par)),max(max(this.A_idx)));
            
        end
        
        function markShortArcs(this, arc_length)
            % mark arc equal or shorter than arc_length as outlier
            %
            % SYNTAX:
            %    this.markShortArcs(arc_length)
            if sum(this.param_class == this.PAR_AMB) > 0
                amb = double(this.A_idx(~this.outlier_obs,this.param_class == this.PAR_AMB));
                [occ,u_amb]=hist(amb,unique(amb));
                to_mark = find(occ <= arc_length);
                for m = to_mark
                    idx_o = amb == u_amb(m);
                    this.outlier_obs(idx_o) = true;
                end
            end
        end
        
        function markSingledObs(this)
            % mark single satellite or receiver epoch as outlier
            %
            % SYNTAX:
            %    this.markSingledObs()
            if sum(this.param_class == this.PAR_IONO)> 0 & this.ls_parametrization.iono(2) == LS_Parametrization.SING_REC % remove phase obe that does not have two pseudorange fro different frequency and phase of a seconde frequency too
                idx_valid_obs = find(~this.outlier_obs);
                idx_valid_ph = idx_valid_obs(this.phase_obs(idx_valid_obs) > 0);
                idx_valid_pr = idx_valid_obs(this.phase_obs(idx_valid_obs) == 0);
                for r = 1 : length(this.unique_rec_name)
                    idx_valid_rec_ph = idx_valid_ph(this.receiver_obs(idx_valid_ph) == r);
                    idx_valid_rec_pr = idx_valid_pr(this.receiver_obs(idx_valid_pr) == r);
                    idx_valid_rec_obs = idx_valid_obs(this.receiver_obs(idx_valid_obs) == r);


                    
                    u_go_id_r = unique(this.satellite_obs(idx_valid_rec_obs));
                    n_sat = length(u_go_id_r);
                    
                    u_oc_ph = unique(this.obs_codes_id_obs(idx_valid_rec_ph));
                    n_ch_ph = length(u_oc_ph);
                    fr_ph = zeros(size(u_oc_ph));
                    for i = 1 : length(u_oc_ph)
                        fr_ph(i) = this.unique_obs_codes{u_oc_ph(i)}(3);
                    end
                    u_fr_ph = unique(fr_ph);
                    n_fr_ph = length(u_fr_ph);
                    
                    u_oc_pr = unique(this.obs_codes_id_obs(idx_valid_rec_pr));
                    n_ch_pr = length(u_oc_pr);
                    fr_pr = zeros(size(u_oc_pr));
                    for i = 1 : length(u_oc_pr)
                        fr_pr(i) = this.unique_obs_codes{u_oc_pr(i)}(3);
                    end
                    u_fr_pr = unique(fr_pr);
                    n_fr_pr = length(u_fr_pr);
                    
                    
                    time_res = min(this.ref_time_obs(idx_valid_rec_obs));
                    duration = max(this.ref_time_obs(idx_valid_rec_obs)) - time_res;
                    time_res = time_res + (0:this.obs_rate:duration);
                    % build phase matrix
                    
                    ph_pres = false(length(time_res),n_sat,n_fr_ph);
                    ph_id = zeros(length(time_res),n_sat,n_ch_ph,'uint32');
                    for c = 1 : n_ch_ph
                        ch = u_oc_ph(c);
                        f = find(fr_ph(c) == u_fr_ph);
                        idx_valid_rec_ph_ch = idx_valid_rec_ph(this.obs_codes_id_obs(idx_valid_rec_ph) == ch);
                        for s = 1 : n_sat
                            sat = u_go_id_r(s);
                            idx_ph = idx_valid_rec_ph_ch(this.satellite_obs(idx_valid_rec_ph_ch) == sat);
                            
                            if any(idx_ph)
                                [~,idx_time] = ismember(this.ref_time_obs(idx_ph),time_res);
                                ph_pres(idx_time, s,f) = true;
                                ph_id(idx_time, s,c) = find(idx_ph);
                            end
                        end
                    end
                    % build pseudorange matrix
                    pr_pres = false(length(time_res),n_sat,n_fr_pr);
                    pr_id = zeros(length(time_res),n_sat,n_ch_pr,'uint32');
                    for c = 1 : n_ch_pr
                        ch = u_oc_pr(c);
                        f = find(fr_pr(c) == u_fr_pr);
                        idx_valid_rec_pr_ch = idx_valid_rec_pr(this.obs_codes_id_obs(idx_valid_rec_pr) == ch);
                        
                        for s = 1 : n_sat
                            sat = u_go_id_r(s);
                            idx_pr = idx_valid_rec_pr_ch(this.satellite_obs(idx_valid_rec_pr_ch) == sat);
                            
                            if any(idx_pr)
                                [~,idx_time] =ismember(this.ref_time_obs(idx_pr),time_res);
                                pr_pres(idx_time, s,f) = true;
                                pr_id(idx_time, s,c) = find(idx_pr);
                            end
                        end
                    end
                    idx_out_pr = sum(pr_pres,3) == 1;
                    idx_out_ph = find((sum(pr_pres,3) < 2 & sum(ph_pres,3) > 1) | sum(ph_pres,3) == 1);

                    idx_out_pr = find(idx_out_pr);
                    
                    for c = 1 : n_ch_pr
                        o_idx = pr_id(idx_out_pr+(c-1) * n_sat*length(time_res));
                        this.outlier_obs(noZero(o_idx)) = true;
                    end
                    for c = 1 : n_ch_ph
                        o_idx = ph_id(idx_out_ph+(c-1) * n_sat*length(time_res));
                        this.outlier_obs(noZero(o_idx)) = true;
                    end
                    
                end
             
            end
            
            for par = [ this.PAR_SAT_CLK this.PAR_SAT_CLK_PR this.PAR_SAT_CLK_PH]
                idx_par_clk = find(this.class_par == par & ~this.out_par);
                id_p_clk = this.param_class == par;
                if any(id_p_clk)
                for s = this.unique_sat_goid
                    idx_par = idx_par_clk(this.sat_par(idx_par_clk) == s);
                    idx_obs = find(this.satellite_obs == s & ~this.outlier_obs);
                    idx_obs = idx_obs(this.A(idx_obs, id_p_clk) ~= 0);
                    num_rec = zeros(length(idx_par),1);
                    for r = 1 : length(this.unique_rec_name)
                        [~,idx_par_rec] = ismember(unique(this.A_idx(idx_obs(this.receiver_obs(idx_obs) == r),id_p_clk )),idx_par);
                        num_rec(noZero(idx_par_rec)) = num_rec(noZero(idx_par_rec)) + 1;
                    end
                    idx_el = find(num_rec < 2);
                    par_sat_clk = this.A_idx(idx_obs,id_p_clk);
                    for i = idx_el'
                        obs_par = par_sat_clk == idx_par(i);
                        this.outlier_obs(idx_obs(obs_par)) = true;
                    end
                end
                end
            end
            for par = [ this.PAR_REC_CLK this.PAR_REC_CLK_PR this.PAR_REC_CLK_PH]
                idx_par_clk = find(this.class_par == par & ~this.out_par);
                id_p_clk = this.param_class == par;
                                if any(id_p_clk)

                for r = 1 : length(this.unique_rec_name)
                    idx_par =  idx_par_clk(this.rec_par(idx_par_clk) == r);
                    idx_obs = find(this.receiver_obs == r & ~this.outlier_obs);
                    idx_obs = idx_obs(this.A(idx_obs, id_p_clk) ~= 0);
                    num_sat = zeros(length(idx_par),1);
                    for s = this.unique_sat_goid
                        [~,idx_par_rec] = ismember(unique(this.A_idx(idx_obs(this.satellite_obs(idx_obs) == s), id_p_clk)), idx_par);
                        num_sat(noZero(idx_par_rec)) = num_sat(noZero(idx_par_rec)) + 1;
                    end
                    idx_el = find(num_sat < 2);
                    par_rec_clk = this.A_idx(idx_obs, id_p_clk);
                    for i = idx_el'
                        obs_par = par_rec_clk == idx_par(i);
                        this.outlier_obs(idx_obs(obs_par)) = true;
                    end
                end
                                end
            end
            
            
            
        end
        
        function removeFullRankDeficency(this)
            % solve full rank deficency removing parameters from the
            % estimation
            %
            % SYNTAX:
            %    this.removeFullRankDeficency()
            
            % TODO make a log of what is taken off to solve the rank
            % deficency
            
            log = Core.getLogger;
            log.addMessage(log.indent('Solving the system'));
            cc = Core.getConstellationCollector();
            
            % remove two bias per receiver
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            this.idx_rd = []; %empty previous par choosen to solve the rank deficency
            idx_rm = [];
            u_ep = unique(this.time_par);
          
                        % remove one clock per epoch
                        if sum(this.param_class == this.PAR_REC_CLK) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                            idx_rec_par = find(this.class_par == this.PAR_REC_CLK  & ~this.out_par);
                            time_par_tmp = this.time_par(idx_rec_par,1);
                            for e = u_ep'
                                idx_par = idx_rec_par(time_par_tmp == e);
                                if length(idx_par)>1
                                    [~,idx_rm_rm] = min(this.rec_par(idx_par));
                                    idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                                end
                            end
                        end
                        if sum(this.param_class == this.PAR_REC_CLK_PH) > 0 && sum(this.param_class == this.PAR_SAT_CLK_PH) > 0
                            idx_rec_par = find(this.class_par == this.PAR_REC_CLK_PH  & ~this.out_par);
                            time_par_tmp = this.time_par(idx_rec_par,1);
                            for e = u_ep'
                                idx_par = idx_rec_par(time_par_tmp == e);
                                if length(idx_par)>1
                                    [~,idx_rm_rm] = min(this.rec_par(idx_par));
                                    idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                                end
                            end
                        end
                        if sum(this.param_class == this.PAR_REC_CLK_PR) > 0 && sum(this.param_class == this.PAR_SAT_CLK_PR) > 0
                            idx_rec_par = find(this.class_par == this.PAR_REC_CLK_PR  & ~this.out_par);
                            time_par_tmp = this.time_par(idx_rec_par,1);
                            for e = u_ep'
                                idx_par = idx_rec_par(time_par_tmp == e);
                                if length(idx_par)>1
                                    [~,idx_rm_rm] = min(this.rec_par(idx_par));
                                    idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                                end
                            end
                        end
                        
             
            
            if sum(this.param_class == this.PAR_REC_CLK) > 0 || sum(this.param_class == this.PAR_REC_CLK_PH) > 0 
                idx_rc = find(this.class_par == this.PAR_REC_CLK & ~this.out_par);
                idx_pr = find(~this.outlier_obs & ~this.phase_obs);
                for r = 1 : n_rec
                    idx_rcr = idx_rc(this.rec_par(idx_rc) == r);                    
                    idx_o_r = idx_pr(this.receiver_obs(idx_pr) == r);
                    if ~isempty(idx_o_r)
                        idx_po = unique(this.A_idx(idx_o_r,this.param_class == this.PAR_REC_CLK));
                        idx_rm = [idx_rm; setdiff(idx_rcr,idx_po)]; % clock param that are not linked to no code obsevration
                    end
                    
                end
            end
            
            if sum(this.param_class == this.PAR_SAT_CLK) > 0 || sum(this.param_class == this.PAR_SAT_CLK_PH) > 0
                idx_rc = find(this.class_par == this.PAR_SAT_CLK & ~this.out_par);
                idx_pr = find(~this.outlier_obs & ~this.phase_obs);
                for s = 1 : n_sat
                    idx_rcr = idx_rc(this.sat_par(idx_rc) == s);
                    
                    
                    idx_o_r = idx_pr(this.satellite_obs(idx_pr) == s);
                    if ~isempty(idx_o_r)
                        idx_po = unique(this.A_idx(idx_o_r,this.param_class == this.PAR_SAT_CLK));
                        idx_rm = [idx_rm; setdiff(idx_rcr,idx_po)]; % clock param that are not linked to no code obsevration
                    end
                end
            end
            
            for p = [this.PAR_REC_CLK this.PAR_REC_CLK_PH this.PAR_REC_CLK_PR ]
                if sum(this.param_class == p) > 0 & sum(this.param_class == this.PAR_SAT_EB) > 0
                    idx_par = find(this.class_par == p & this.rec_par == 1);
                    idx_rm = [idx_rm; idx_par(1)];
                end
            end
            
            % ---- (multi receiver) for each epoche remove one coordinate --------
            if (sum(this.param_class == this.PAR_REC_X) > 0 || sum(this.param_class == this.PAR_REC_Y) > 0  ||  sum(this.param_class == this.PAR_REC_Z) > 0 || sum(this.param_class == this.PAR_TROPO) > 0  || sum(this.param_class == this.PAR_TROPO_E) > 0  || sum(this.param_class == this.PAR_TROPO_N) > 0) && (sum(this.param_class == this.PAR_SAT_CLK) > 0 || sum(this.param_class == this.PAR_SAT_CLK_PH) > 0 || sum(this.param_class == this.PAR_SAT_CLK_PR) > 0)
                idx_time_x = find(this.class_par == this.PAR_REC_X & ~this.out_par);
                time_tmp_x = this.time_par(idx_time_x,:);
                idx_time_y = find(this.class_par == this.PAR_REC_Y & ~this.out_par);
                time_tmp_y = this.time_par(idx_time_y,:);
                idx_time_z = find(this.class_par == this.PAR_REC_Z & ~this.out_par);
                time_tmp_z = this.time_par(idx_time_z,:);
                if ~this.free_tropo
                    idx_time_t = find(this.class_par == this.PAR_TROPO & ~this.out_par);
                    time_tmp_t = this.time_par(idx_time_t,:);
                    idx_time_e = find(this.class_par == this.PAR_TROPO_E & ~this.out_par);
                    time_tmp_e = this.time_par(idx_time_e,:);
                    idx_time_n = find(this.class_par == this.PAR_TROPO_N & ~this.out_par);
                    time_tmp_n = this.time_par(idx_time_n,:);
                end
                for e = u_ep'
                    idx_par = idx_time_x(time_tmp_x(:,1) >= e & ( time_tmp_x(:,1)==0 | time_tmp_x(:,2) < e));
                    if any(any(idx_par))
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    
                    idx_par = idx_time_y(time_tmp_y(:,1) >= e & ( time_tmp_y(:,1)==0 | time_tmp_y(:,2) < e));
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    
                    idx_par = idx_time_z(time_tmp_z(:,1) >= e & ( time_tmp_z(:,1)==0 | time_tmp_z(:,2) < e));
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                        end
                    end
                    if ~this.free_tropo
                        idx_par = idx_time_t(time_tmp_t(:,1) >= e & ( time_tmp_t(:,1)==0 | time_tmp_t(:,2) < e));
                        if any(idx_par)
                            [~,idx_rm_rm] = min(this.rec_par(idx_par));
                            if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                                idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                            end
                        end
                        
                        idx_par = idx_time_e(time_tmp_e(:,1) >= e & ( time_tmp_e(:,1)==0 | time_tmp_e(:,2) < e));
                        if any(idx_par)
                            [~,idx_rm_rm] = min(this.rec_par(idx_par));
                            if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                                idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                            end
                        end
                        
                        idx_par = idx_time_n(time_tmp_n(:,1) >= e & ( time_tmp_n(:,1)==0 | time_tmp_n(:,2) < e));
                        if any(idx_par)
                            [~,idx_rm_rm] = min(this.rec_par(idx_par));
                            if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                                idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                            end
                        end
                    end
                end
            end
            this.idx_rd = unique(noZero(idx_rm));
        end
        
        
        function removeEstParam(this, idx)
            % reomve and estimated paramter from the system
            %
            % SYNTAX:
            %   this.removeEstParam(idx)
             A = this.A_full(1: length(this.obs),idx);
             this.obs = this.obs - A * this.x(idx);
             this.remPar(idx);
        end
        
        
        function remPar(this,idx)
            % remove paramter form the system
            %
            % SYNTAX:
            %    this.remPar(idx)
            this.x(idx) = [];
            n_par = length(this.class_par);
            this.class_par(idx) = [];
            this.time_par(idx,:) = [];
            this.param_par(idx,:) = [];
            this.rec_par(idx) = [];
            this.sat_par(idx) = [];
            this.obs_codes_id_par(idx) = [];
            this.wl_id_par(idx) = [];
            this.out_par(idx) = [];
            this.phase_par(idx) = [];
            new_par_num = zeros(n_par,1);
            new_par_num(~idx) = 1 : (n_par - sum(idx));
            n_zero = this.A_idx ~= 0;
            this.A_idx(n_zero) = new_par_num(this.A_idx(n_zero));
            n_zero = this.A_idx_pseudo ~= 0;
            this.A_idx_pseudo(n_zero) = new_par_num(this.A_idx_pseudo(n_zero));
        end
        
        function absValRegularization(this,p_class, var)
            % regularize parameters to zero (Tykhnov aka ridge aka L2)
            %
            % this.absValRegularization(param_id, var)
            par_ids = this.param_class == p_class;
            u_p_id = unique(this.A_idx(:,par_ids));
            n_par = length(u_p_id);
            A_idx_tmp = zeros(n_par, 2,'uint32'); % tykhonv regualrization are now limited to two parameters
            A_tmp = zeros(n_par, 2);
            A_tmp(:,1) = 1;
            A_idx_tmp(:,1) = u_p_id;
            this.A_pseudo = [this.A_pseudo; A_tmp];
            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
            this.variance_pseudo = [this.variance_pseudo; var*ones(n_par,1)];
            this.receiver_pseudo = [this.receiver_pseudo; this.rec_par(u_p_id)];
            if isempty(this.time_pseudo)
                this.time_pseudo = GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id));
            else
                this.time_pseudo.addEpoch(GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id)));
            end
            this.satellite_pseudo = [this.satellite_pseudo; this.sat_par(u_p_id)];
        end
        
        function timeRegularization(this, p_class, var_per_sec)
            % first order tykhonv regualrization in time
            %
            % this.timeRegularization(this, param_id, var)
            par_ids = this.param_class == p_class;
            p_idx = unique(find(this.class_par == p_class));
            rec_idx = this.rec_par(p_idx);
            u_rec = unique(rec_idx);
            sat_idx = this.sat_par(p_idx);
            u_sat = unique(sat_idx);
            ch_idx  =  this.obs_codes_id_par(p_idx);
            u_ch = unique(ch_idx);
            for r = u_rec'
                for s = u_sat'
                    for c = u_ch'
                        p_tmp = p_idx(rec_idx == r & sat_idx == s & ch_idx == c); % find the idx of the obsrevations
                        if ~isempty(p_tmp)
                            u_p_id = unique(p_tmp);
                            n_par = length(u_p_id);
                            A_tmp = [ones(n_par-1,1) -ones(n_par-1,1)];
                            A_idx_tmp = [u_p_id(1:(end-1)) u_p_id(2:end)];
                            this.A_pseudo = [this.A_pseudo; A_tmp];
                            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
                            
                            this.variance_pseudo = [this.variance_pseudo; var_per_sec * diff(double(this.time_par(u_p_id)))];
                            % taking the indices of first epoch
                            this.receiver_pseudo = [this.receiver_pseudo; this.rec_par(u_p_id(1:end-1))];
                            if isempty( this.time_pseudo)
                                this.time_pseudo = GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id(1:end-1)));
                            else
                                this.time_pseudo.addEpoch(this.time_min.getMatlabTime + this.time_par(u_p_id(1:end-1)));
                            end
                            this.satellite_pseudo = [this.satellite_pseudo; this.sat_par(u_p_id(1:end-1))];
                        end
                    end
                end
            end
        end
        
        function spatialRegularization(this, law)
            % tykhonv regualrization in space
            %
            % this.spatialRegularization(this, law)
            
            
            n_rec_tot = size(this.xyz,1);
            dist_mat = zeros(n_rec_tot);
            for i = 1 : n_rec_tot
                for j = (i+1) : n_rec_tot
                    coord1 = Coordinates();
                    coord1.xyz = this.xyz(i,:);
                    coord2 = Coordinates();
                    coord2.xyz = this.xyz(j,:);
                    dist_mat(i,j) = coord1.ellDistanceTo(coord2);
                end
            end
            var_mat = law(dis_mat);
            
            % -- construct A and indices
            par_ids = this.parameter_class == p_class;
            p_idx = unique(this.A_idx(:,par_ids));
            sat_idx = this.sat_par(p_idx);
            u_sat = unique(sat_idx);
            ch_idx  =  this.ch_par(p_idx);
            time_idx = this.time_par(p_idx);
            u_ch = unique(ch_idx);
            for s = u_sat
                for c = u_ch
                    p_tmp = p_idx(rec_idx == r & sta_idx == s & ch_idx == c); % find the idx of the obsrevations
                    u_ep = unique(this.time_par(p_tmp));
                    for e = u_ep
                        p_tmp = p_idx(rec_idx == r & sta_idx == s & ch_idx == c & time_idx == e); % find the idx of the obsrevations
                        if ~isempty(p_tmp)
                            u_p_id = unique(p_tmp);
                            n_par = length(u_p_id);
                            rec_par = this.rec_par(u_p_id);
                            n_rec = length(rec_apr);
                            n_pseudobs = n_rec*(n_rec - 1);
                            A_tmp = [ones(n_par-1,1) -ones(n_par-1,1)];
                            idx_1 = [];
                            for r = 1 : n_rec
                                idx_1 = [idx_1; (r:n_rec)'];
                            end
                            idx_2 = idx_1 + 1;
                            A_idx_tmp = [u_p_id(idx_1) u_p_id(idx_2)];
                            var = var_mat(rec_par(idx_1)+n_rec_tot*(rec_par(idx_1) -1));
                            this.A_pseudo = [this.A_pseudo; A_tmp];
                            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
                            
                            this.variance_pseudo = [this.variance_pseudo; var_per_sec * var];
                            % taking the indices of first epoch
                            this.receiver_pseudo = [this.receiver_pseudo; zeros(n_pseudo_obs,1)];
                            this.time_pseudo.addEpoch(GPS_Time(this.time_min + e*ones(n_pseudo_obs,1)));
                            this.satellite_pseudo = [this.satellite_pseudo; s*ones(n_pseudo_obs,1)];
                        end
                    end
                end
            end
        end
        
        function hemisphereRegularization(this, law)
            % tykhonv regualrization on the hemisphere
            %
            % this.spatialRegularization(this, law)
            % TBD
        end
        
        function reduceForNuisanceParameters(this, param_id)
            % reduce for the parameter
        end
        
        
        function solve(this, fix)
            % %solve the least squares
            % SYNTAX:
            %    this.solve(<fix>)
            % INPUT:
            %    fix : seek to fix integer parameters
            if nargin < 2
                fix = false;
            end
            % ------ mark short arc as outlier
            this.markShortArcs(1);
            % ------ mark single reciver or satellite epoch as outlier
            this.markSingledObs();
            % ------ remove outlier
            this.generateOutliedParIdx();
            % ------ remove full rank deficencies
            this.removeFullRankDeficency();
            % apply elevation dependent weighting
            state = Core.getCurrentSettings();
            if state.getWeigthingStrategy == 2
                this.sinElevationWeigth();
            end
            % ------ form the normal matrix
            n_obs = size(this.A,1) + size(this.A_pseudo,1);
            n_par = double(max(max(this.A_idx)));
            n_rec = size(this.rec_xyz,1);

            rows = repmat((1:size(this.A,1))',1,size(this.A,2));
            rows_pseudo = repmat(size(this.A,1)+(1:size(this.A_pseudo,1))',1,size(this.A_pseudo,2));
            rows = [rows(:); rows_pseudo(:)];
            columns = double(zero2n([this.A_idx(:); this.A_idx_pseudo(:)],1));
            values = [this.A(:); this.A_pseudo(:)];
            % Creating A' instead of A
            % A = sparse(rows, columns, values, n_obs, n_par); % <- this is A
            A = sparse(columns, rows, values, n_par, n_obs); % <- stupid trick for speed-up MATLAB traspose the sparse matrix!!!
            this.A_full = A'; % save it for use in other methods
            n_out = sum(this.outlier_obs);
            A_out = A(:, this.outlier_obs > 0)'; % <- this is 1000 time slower with not trasposed A
            A(:, this.outlier_obs > 0) = [];
            A = A'; % Re-traspose A since it have been created as trasposed
            A(:, [this.idx_rd; find(this.out_par)]) = [];
            
            class_par = this.class_par;
            class_par([this.idx_rd; find(this.out_par)]) = [];
            
            zero_pars = sum(A~=0) == 0;
            A(:,zero_pars) = [];
            class_par(zero_pars) = [];
            
            
            rec_par = this.rec_par;
            rec_par([this.idx_rd; find(this.out_par)]) = [];
            rec_par(zero_pars) = [];
            
            time_par = this.time_par;
            time_par([this.idx_rd; find(this.out_par)],:) = [];
            time_par(zero_pars,:) = [];

            valid_pars = find(~Core_Utils.ordinal2logical([this.idx_rd; find(this.out_par)],n_par)); % sometimes with splines spme paramter have a zero entry
            this.idx_rd = [this.idx_rd; valid_pars(zero_pars)];
            if isempty(this.reweight_obs)
                this.reweight_obs = ones(size(this.variance_obs));
            end
            vars = [1./this.variance_obs(~this.outlier_obs).*this.reweight_obs(~this.outlier_obs); 1./this.variance_pseudo];
            mean_vars = 1 ; %mean(vars);
            vars = vars ./ mean_vars;
            Cyy =  spdiags(vars,0,n_obs - n_out,n_obs - n_out);
            x_est = zeros(n_par -length(this.idx_rd) - sum(this.out_par),1);
            y = sparse([this.obs(~this.outlier_obs); zeros(size(this.A_pseudo,1),1)]);

            Aw = A'*Cyy;
            
            %clearvars Aw
            % ------ reduce for sat clock, rec clock and iono
            idx_reduce_sat_clk = class_par == this.PAR_SAT_CLK | class_par == this.PAR_SAT_CLK_PH | class_par == this.PAR_SAT_CLK_PR;
            idx_reduce_rec_clk = class_par == this.PAR_REC_CLK | class_par == this.PAR_REC_CLK_PH | class_par == this.PAR_REC_CLK_PR;
            idx_reduce_iono = class_par == this.PAR_IONO;
            
            idx_reduce = idx_reduce_rec_clk | idx_reduce_sat_clk | idx_reduce_iono;
            
            N = Aw(~idx_reduce,:)*A(:,~idx_reduce);
            B = Aw(~idx_reduce,:)*y;
            
            max_ep = max(this.ref_time_obs);
            ref_time_obs = this.ref_time_obs(~this.outlier_obs);
            step = 7200;
            time_par_red = time_par(idx_reduce,1);
            
            iono = sum(idx_reduce_iono) > 0;
            sat_clk = sum(idx_reduce_sat_clk) > 0;
            rec_clk = sum(idx_reduce_rec_clk) > 0;

            cross_terms = {};
            ii  = 1;
            for i = 0 : step : floor(max_ep / step) * step % sparse matrix library became very slow in case of big/huge matrix the reduction can be applyed dividing the matrices in parts
                idx_time_obs = ref_time_obs > (i-this.obs_rate/10) &  ref_time_obs < (i + step -this.obs_rate/10 );
                idx_time_par_red = time_par_red > (i-this.obs_rate/10)  &  time_par_red < (i + step -this.obs_rate/10 );
                idx_red_cycle = idx_reduce;
                idx_red_cycle(idx_red_cycle) = idx_time_par_red;
                
                cp_cycle = class_par(idx_red_cycle);
                rp_cycle = rec_par(idx_red_cycle);
     
                idx_reduce_cycle_sat_clk = cp_cycle == this.PAR_SAT_CLK | cp_cycle == this.PAR_SAT_CLK_PH | cp_cycle == this.PAR_SAT_CLK_PR;
                idx_reduce_cycle_rec_clk = cp_cycle == this.PAR_REC_CLK | cp_cycle == this.PAR_REC_CLK_PH | cp_cycle == this.PAR_REC_CLK_PR;
                idx_reduce_cycle_iono = cp_cycle == this.PAR_IONO;
                
                
                
                
                Awr_t = Aw( idx_red_cycle, idx_time_obs);
                Ar_t = A(idx_time_obs, idx_red_cycle);
                Ae_t = A(idx_time_obs, ~idx_reduce);
                y_t = y(idx_time_obs);
                
                Nr_t = Awr_t*Ar_t;
                Br_t = Awr_t*y_t;
                
                Ner_t = Awr_t*Ae_t;
                
                if iono
                    n_iono = sum(idx_reduce_cycle_iono);
                    iIono = spdiags(1./diag(Nr_t(idx_reduce_cycle_iono, idx_reduce_cycle_iono)),0,n_iono,n_iono);
                    Nx_iono = Ner_t(idx_reduce_cycle_iono, :); % cross term reduce iono
                    Nx_iono_cycle = Nr_t(~idx_reduce_cycle_iono, idx_reduce_cycle_iono); % cross term reduce iono
                    Nt = Nx_iono' * iIono;
                    Nt_cycle = Nx_iono_cycle * iIono;
                    N = N - Nt * Nx_iono;
                    Nr_t = Nr_t(~idx_reduce_cycle_iono,~idx_reduce_cycle_iono) - Nt_cycle * Nr_t(idx_reduce_cycle_iono,~idx_reduce_cycle_iono);
                    Ner_t(~idx_reduce_cycle_iono, :) = Ner_t(~idx_reduce_cycle_iono, :) - Nt_cycle*Nx_iono;
                    Ner_t(idx_reduce_cycle_iono, :) = [];

                    B_iono = Br_t(idx_reduce_cycle_iono);
                    B = B - Nt * B_iono;
                    Br_t = Br_t(~idx_reduce_cycle_iono) - Nt_cycle*B_iono;
                    
                    cross_terms_t{1} = {iIono B_iono [Nx_iono Nx_iono_cycle'] idx_reduce_cycle_iono};
                end
                
                if sat_clk
                    i_sat_clk_tmp = idx_reduce_cycle_sat_clk(~idx_reduce_cycle_iono);
                    cp = cp_cycle( ~idx_reduce_cycle_iono);
                    idx_1 = cp(i_sat_clk_tmp) == this.PAR_SAT_CLK | cp(i_sat_clk_tmp) == this.PAR_SAT_CLK_PR;
                    idx_2 = cp(i_sat_clk_tmp) == this.PAR_SAT_CLK_PH;
                    if sum(idx_2) > 0 
                        iSatClk = Core_Utils.inverseByPartsDiag(Nr_t(i_sat_clk_tmp,i_sat_clk_tmp),idx_1, idx_2);%inv(N(i_sat_clk_tmp,i_sat_clk_tmp))  ;%;%spdiags(1./diag(N(i_sat_clk_tmp,i_sat_clk_tmp)),0,n_clk_sat,n_clk_sat);
                    else
                        n_sat_clk = sum(i_sat_clk_tmp);
                        iSatClk = spdiags(1./diag(Nr_t(i_sat_clk_tmp, i_sat_clk_tmp)),0,n_sat_clk,n_sat_clk);
                    end
                    Nx_satclk = Ner_t(i_sat_clk_tmp, :);
                    Nx_satclk_cyle = Nr_t(~i_sat_clk_tmp, i_sat_clk_tmp);
                    idx_full = sum(Nx_satclk~=0,1) >0;
                    Nt = Nx_satclk(:,idx_full)' * iSatClk;
                    Nt_cycle = Nx_satclk_cyle * iSatClk;
                    
                    N(idx_full,idx_full) = N(idx_full,idx_full) - sparse(full(Nt) * full(Nx_satclk(:,idx_full)));
                    Nr_t = Nr_t(~i_sat_clk_tmp,~i_sat_clk_tmp) - Nt_cycle * Nr_t(i_sat_clk_tmp, ~i_sat_clk_tmp);
                    Ner_t(~i_sat_clk_tmp, :) = Ner_t(~i_sat_clk_tmp, :) - Nt_cycle*Nx_satclk;
                    Ner_t(i_sat_clk_tmp, :) = [];
                    
                    
                    B_satclk =  Br_t(i_sat_clk_tmp);
                    B(idx_full) = B(idx_full) - Nt * B_satclk;
                    Br_t = Br_t(~i_sat_clk_tmp) - Nt_cycle * B_satclk;
                    
                    cross_terms_t{2} = {iSatClk B_satclk [Nx_satclk Nx_satclk_cyle'] idx_reduce_cycle_sat_clk};
                    
                end
                
                if rec_clk
                    i_rec_clk_tmp = idx_reduce_cycle_rec_clk(~idx_reduce_cycle_iono & ~idx_reduce_cycle_sat_clk);
                    n_rec_clk = sum(i_rec_clk_tmp);
                    if n_rec > 1 
                        if sat_clk
                        rp = rp_cycle( ~idx_reduce_cycle_sat_clk &  ~idx_reduce_cycle_iono);
                        cp = cp_cycle( ~idx_reduce_cycle_sat_clk &  ~idx_reduce_cycle_iono);
                        
                        indices = {};
                        for r = 2 : n_rec
                            indices{(r-2)*2+1} = rp(i_rec_clk_tmp) == r & (cp(i_rec_clk_tmp) == this.PAR_REC_CLK | cp(i_rec_clk_tmp) == this.PAR_REC_CLK_PR);
                            indices{(r-1)*2}   = rp(i_rec_clk_tmp) == r & cp(i_rec_clk_tmp) == this.PAR_REC_CLK_PH;
                        end
                        n_i = length(indices);
                        jj = 1;
                        while  jj <= n_i
                            if sum(indices{jj}) == 0
                                indices(jj) = [];
                                jj = jj - 1;
                                n_i = n_i - 1;
                            end
                            jj = jj +1;
                        end
                        if length(indices) > 1
                            iRecClk = Core_Utils.inverseByRecPartsDiag(Nr_t(i_rec_clk_tmp,i_rec_clk_tmp),indices);
                        else
                            iRecClk = spdiags(1./diag(Nr_t(i_rec_clk_tmp, i_rec_clk_tmp)),0,n_rec_clk,n_rec_clk);
                        end
                        else
                            cp = cp_cycle( ~idx_reduce_cycle_iono);
                            idx_1 = cp(i_rec_clk_tmp) == this.PAR_REC_CLK | cp(i_rec_clk_tmp) == this.PAR_REC_CLK_PR;
                            idx_2 = cp(i_rec_clk_tmp) == this.PAR_REC_CLK_PH;
                            if sum(idx_2) > 0
                                iRecClk = Core_Utils.inverseByPartsDiag(Nr_t(i_rec_clk_tmp,i_rec_clk_tmp),idx_1, idx_2);
                            else
                                iRecClk = spdiags(1./diag(Nr_t(i_rec_clk_tmp, i_rec_clk_tmp)),0,n_rec_clk,n_rec_clk);
                            end
                        end
                    else
                        iRecClk = inv(Nr_t(i_rec_clk_tmp,i_rec_clk_tmp));
                    end
                    
                    Nx_recclk = Ner_t(i_rec_clk_tmp, :);
                    idx_full = sum(Nx_recclk~=0,1) >0;

                    Nt = Nx_recclk(:,idx_full)' * iRecClk;
                    N(idx_full,idx_full) = N(idx_full,idx_full) - sparse(full(Nt) * full(Nx_recclk(:,idx_full)));
                    
                    B_recclk = Br_t(i_rec_clk_tmp);
                    B(idx_full) = B(idx_full) - Nt * B_recclk;
                    
                    cross_terms_t{3} = {iRecClk B_recclk Nx_recclk idx_reduce_cycle_rec_clk};
                end
                cross_terms{ii} = {cross_terms_t idx_red_cycle}; 
               ii = ii +1;
            end
            
            % ------- fix the ambiguities
            idx_amb = class_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) == this.PAR_AMB;
            ldl_strategy = false;
            if sum(this.param_class == this.PAR_AMB) > 0 && fix && any(idx_amb)
                    if ldl_strategy
                    
                    % reducing for everything else except ambiguities
                    idx_b = find(~idx_amb);
                    
                    [L,D,p] = ldl(N(~idx_amb, ~idx_amb),'vector');
                    tol = 1e-3;
                    rm_id_b = (diag(D) < tol); % find amb to be removed
                    iL = inv(L(~rm_id_b, ~rm_id_b));
                    iD = spdiags(1./diag(D(~rm_id_b, ~rm_id_b)),0,sum(~rm_id_b),sum(~rm_id_b));     
                    C_bb = zeros(sum(~idx_amb),sum(~idx_amb));
                    C_bb(p(~rm_id_b),p(~rm_id_b)) = full(iL'*iD)*full(iL);
%                     Nbb = N(~idx_amb, ~idx_amb);
%                     C_bb(p(~rm_id_b),p(~rm_id_b)) = inv(Nbb(p(~rm_id_b),p(~rm_id_b)));
                    %                         %Bbb =  B(~idx_amb);
                    BB = full(N(idx_amb,idx_b ))*C_bb;
                    N_amb_amb = N(idx_amb, idx_amb) - sparse(BB*full(N(idx_b, idx_amb)));
                    B_amb_amb = B(idx_amb) -  BB* B(~idx_amb);
                    
                    
                    % LDL decompisition
                    [L,D,P] = ldl(N_amb_amb);
                    % permute
                    Baa = P'*B_amb_amb;
                    
                    % compute threshold
                    tol = 10; %1e7*eps(max(abs(diag(D))));
                    
                    % find the hard rank deficency
                    rm_id = (diag(D) < tol); % find amb to be removed
                    
                    
                    % compute the VCV for the ambiguoties an their
                    % flost solution
                    iL = inv(L(~rm_id,~rm_id));
                    iD = inv(D(~rm_id,~rm_id));
                    C_amb_amb = full(iL'*iD)*full(iL);
                    Baa = Baa(~rm_id);
                    amb_float = C_amb_amb * Baa;
                    
                    
                    C_amb_amb_r = C_amb_amb;
                    
                    % try to fix
                    
                    % fix by receiver
                    ra = rec_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono);
                    ra = int16(P'*double(ra(idx_amb)));
                    ra(rm_id) = [];
                    not_rm_id = ~Core_Utils.ordinal2logical(rm_id,size(N_amb_amb,1));
                    ambs = zeros(size(N_amb_amb,1),1);
                    is_fixed_amb = false(size(ambs));
                    is_fixed_amb((P*rm_id) > 0) = true;
                    for r = n_rec : -4 : 1            
                        rec_amb_id = ra == r | ra == max(1,(r-1)) | ra == max(1,(r-2)) | ra == max(1,(r-3));
                        if sum(rec_amb_id) > 0
                            amb_rec = amb_float(rec_amb_id);
                            [amb_fixed, is_fixed, l_fixed] = Fixer.fix(amb_rec, C_amb_amb_r(rec_amb_id,rec_amb_id), 'lambda_partial' );
                            rec_amb_id_fix = rec_amb_id;
                            rec_amb_id_fix(rec_amb_id_fix) = l_fixed(:,1);
                            if is_fixed 
                                amb_rec(l_fixed(:,1)) = amb_fixed(l_fixed(:,1),1);
                                iaa = inv(C_amb_amb_r(rec_amb_id_fix,rec_amb_id_fix));
                                 C_amb_amb_r(~rec_amb_id,~rec_amb_id) = C_amb_amb_r(~rec_amb_id,~rec_amb_id) - C_amb_amb_r(~rec_amb_id,rec_amb_id_fix)*iaa*C_amb_amb_r(rec_amb_id_fix,~rec_amb_id);
                                 amb_float(~rec_amb_id) = amb_float(~rec_amb_id) - C_amb_amb_r(~rec_amb_id,rec_amb_id_fix)*iaa*(amb_float(rec_amb_id_fix) - amb_rec(l_fixed(:,1)));
                            end
                            rec_amb_id_fix = not_rm_id;
                            rec_amb_id_fix(not_rm_id) = rec_amb_id;
                            ambs(rec_amb_id_fix) = amb_rec;
                            is_fixed_amb(rec_amb_id_fix) = l_fixed(:,1);
                        end
                    end
                    
                    ambs = P*ambs;
                    
                    cp_red = class_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono);
                    idx_x = find(cp_red  == this.PAR_REC_X);
                    idx_y = find(cp_red  == this.PAR_REC_Y);
                    idx_z = find(cp_red  == this.PAR_REC_Z);
                    if sum(cp_red  == this.PAR_REC_X)> 0 & (length(idx_x) == length(idx_y))& (length(idx_x) == length(idx_z)) & this.ls_parametrization.rec_x(1) == LS_Parametrization.CONST
                        C_amb_amb = P(~rm_id,~rm_id)*C_amb_amb*P(~rm_id,~rm_id)';
                        is_fixed_amb = P*is_fixed_amb;
                        coo_vcv_B = sparse(zeros(size(N,1),length(idx_x)*3));
                        for c = 1 : length(idx_x)
                            coo_vcv_B( idx_x(c),(c-1)*3 + 1) = 1;
                            coo_vcv_B( idx_y(c),(c-1)*3 + 2) = 1;
                            coo_vcv_B( idx_z(c),(c)*3) = 1;
                        end
                        C_amb_amb = C_amb_amb(~is_fixed_amb,~is_fixed_amb);
                        nn = N(idx_b,idx_amb);
                        nn = nn(:,~is_fixed_amb);
                        N(idx_b,idx_b) = N(idx_b,idx_b) -  nn*C_amb_amb*nn';
                        results = N \ coo_vcv_B;
                        cp_red = cp_red(idx_b);
                        idx_coo = cp_red == this.PAR_REC_X | cp_red == this.PAR_REC_Y | cp_red == this.PAR_REC_Z;
                        this.coo_vcv = results(idx_coo,:);
                    end
                    B(idx_b) = B(idx_b) - N(idx_b,idx_amb)*ambs;
                    not_ambs = C_bb*B(idx_b);
                    x_reduced = zeros(size(N,1),1);
                    x_reduced(~idx_amb) = not_ambs;
                    x_reduced(idx_amb) = ambs;
                    else
                        % svd startegy
                        cod_avail =  false & exist('spqr') ;
                        % reduce all other paramter than ambiguoties
                        if cod_avail
                            cod_fact = factorization_cod_sparse(N(~idx_amb, ~idx_amb));
                            BB = cod_fact\[N(~idx_amb, idx_amb) B(~idx_amb)];
                            N_amb_amb = N(idx_amb, idx_amb) - N(idx_amb,~idx_amb )*BB(:,1:(end-1));
                            B_amb_amb = B(idx_amb) - N(idx_amb,~idx_amb )*BB(:,end);
                        else
                            [U,D,V] = svds(N(~idx_amb, ~idx_amb),sum(~idx_amb));
                            d = diag(D);
                            tol = max(size(N(~idx_amb, ~idx_amb))) * eps(norm(diag(D),inf));
                            keep_id = d > (tol);
                            C_bb = U(:,keep_id) * inv(D(keep_id,keep_id)) * V(:,keep_id)';
                            BB = full(N(idx_amb,~idx_amb ))*C_bb;
                            N_amb_amb = N(idx_amb, idx_amb) - sparse(BB*full(N(~idx_amb, idx_amb)));
                            B_amb_amb = B(idx_amb) -  BB*B(~idx_amb);
                        end
                        
                        
                        [C_amb_amb, amb_float,idx_amb_est] = LS_Manipulator_new.getEstimableAmb(N_amb_amb, B_amb_amb);
                        [amb_fixed, is_fixed, l_fixed] = Fixer.fix(amb_float, C_amb_amb, 'lambda_partial' );
                        ambs = zeros(sum(idx_amb),1);
                        ambs(idx_amb_est) = amb_fixed(:,1);

                        B(~idx_amb) = B(~idx_amb) - N(~idx_amb,idx_amb)*ambs;
                         if cod_avail
                             BB = cod_fact\[N(~idx_amb, idx_amb) B(~idx_amb)];
                         else
                        not_ambs = C_bb*B(~idx_amb);
                         end
                        x_reduced = zeros(size(N,1),1);
                        x_reduced(~idx_amb) = not_ambs;
                        x_reduced(idx_amb) = ambs;

                    end
            else
                cp_red = class_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono);
                idx_x = find(cp_red  == this.PAR_REC_X);
                    idx_y = find(cp_red  == this.PAR_REC_Y);
                    idx_z = find(cp_red  == this.PAR_REC_Z);
                if sum(  cp_red  == this.PAR_REC_X)> 0 & (length(idx_x) == length(idx_y))& (length(idx_x) == length(idx_z)) & this.ls_parametrization.rec_x(1) == LS_Parametrization.CONST && false
                    coo_vcv_B = sparse(zeros(size(N,1),length(idx_x)*3));
                    for c = 1 : length(idx_x)
                         coo_vcv_B( idx_x(c),(c-1)*3+1) = 1;
                         coo_vcv_B(  idx_y(c),(c-1)*3+2) = 1;
                         coo_vcv_B(  idx_z(c),(c)*3) = 1;
                     end
                     result  = N\[B coo_vcv_B];
                     x_reduced = result(:,1);
                     this.coo_vcv = result([idx_x; idx_y; idx_z] ,2 : end);
                else
                     x_reduced = N\B;
                end
            end
            % ------- substitute back
            ii = 1;
            x_est(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) = x_reduced;
            for i = 0 : step : floor(max_ep / step) * step
                % receiver clock
                if rec_clk
                    B_recclk = cross_terms{ii}{1}{3}{2};
                    Nx_recclk = cross_terms{ii}{1}{3}{3};
                    iRecClk = cross_terms{ii}{1}{3}{1};
                    idx_reduce_cycle_rec_clk = cross_terms{ii}{2};
                    idx_reduce_cycle_rec_clk(idx_reduce_cycle_rec_clk) = cross_terms{ii}{1}{3}{4};
                    B_recclk = B_recclk - sum(Nx_recclk * spdiags(x_reduced,0,length(x_reduced),length(x_reduced)),2);
                    x_rec_clk = iRecClk * B_recclk;
                    x_est(idx_reduce_cycle_rec_clk) = x_rec_clk;
                end
                
                % satellite clcok
                if sat_clk
                    B_satclk = cross_terms{ii}{1}{2}{2};
                    Nx_satclk = cross_terms{ii}{1}{2}{3};
                    iSatClk = cross_terms{ii}{1}{2}{1};
                    idx_reduce_cycle_sat_clk = cross_terms{ii}{2};
                    idx_reduce_cycle_sat_clk(idx_reduce_cycle_sat_clk) = cross_terms{ii}{1}{2}{4};
                    idx_est1 = ~idx_reduce ;
                    idx_est2 = idx_reduce_cycle_rec_clk;
                    n_est = sum(idx_est1) + sum(idx_est2);
                    B_satclk = B_satclk -   sum(Nx_satclk * spdiags([x_est(idx_est1); x_est(idx_est2)],0,n_est,n_est),2);
                    x_sat_clk = iSatClk * B_satclk;
                    x_est(idx_reduce_cycle_sat_clk) = x_sat_clk;
                end
                
                % iono
                if iono
                    B_iono = cross_terms{ii}{1}{1}{2};
                    Nx_iono = cross_terms{ii}{1}{1}{3};
                    iIono = cross_terms{ii}{1}{1}{1};
                    idx_reduce_cycle_iono = cross_terms{ii}{2};
                    idx_reduce_cycle_iono(idx_reduce_cycle_iono) = cross_terms{ii}{1}{1}{4};
                    idx_est1 = ~idx_reduce ;
                    if sat_clk
                        idx_est2 = idx_reduce_cycle_rec_clk | idx_reduce_cycle_sat_clk;
                    else
                         idx_est2 = idx_reduce_cycle_rec_clk;
                    end
                    n_est = sum(idx_est1) + sum(idx_est2);
                    B_iono = B_iono -   sum(Nx_iono * spdiags([x_est(idx_est1); x_est(idx_est2)],0,n_est,n_est),2);
                    x_iono = iIono * B_iono;
                    x_est(idx_reduce_cycle_iono) = x_iono;
                end
                ii = ii + 1;
            end
            
            
            x = zeros(n_par,1);
            idx_est = true(n_par,1);
            idx_est([this.idx_rd ; find(this.out_par)]) = false;
            x(idx_est) = x_est;
            res = nan(size(this.obs));
            
            % generate esatimations also for the out par (to get a residual)
            if n_out > 0 && false % to be debugged                            
                res_out = this.obs(this.outlier_obs) - A_out(:,~this.out_par & ~Core_Utils.ordinal2logical(this.idx_rd,n_par))*x_est;
                red_out = res_out;
                A_res_red = A_out(:,this.out_par);
                idx_empty = sum(A_res_red,2) == 0;
                red_out(idx_empty) = [];
                A_res_red(idx_empty,:) = [];
                x_out = A_res_red \ red_out;
                res_out(~idx_empty) = red_out - A_res_red * x_out;
                x(this.out_par) = x_out;
                res(this.outlier_obs) = res_out;
            end
            
            res(~this.outlier_obs) = this.obs(~this.outlier_obs) - A(1:sum(~this.outlier_obs),:)*x_est;
             this.res = res;
            this.x = x;
            
        end
        
        function applyWeightingStrategy()
        end
        
        function res = getResidual(this)
        end
        
        function elevationWeigth(this, fun)
            % weight observation by elevation
            %
            % SYNTAX:
            %   this.elevationWeigth( fun)
            this.variance_obs = this.variance_obs .* fun(this.elevation_obs);
        end
        
        function sinElevationWeigth(this)
            % weight observation by elevation by 1/sin
            %
            % SYNTAX:
            %   this.sinElevationWeigth()
            fun = @(el) (3./sin(el)).^2;
            this.elevationWeigth( fun)
        end
        
        function simpleSnoop(this, ph_thr, pr_thr)
            % simple threshold on residual
            %
            % this.Snoop(this, ph_thr, pr_thr)
            idx_out_ph = this.phase_obs & abs(this.res) > ph_thr;
            this.outlier_obs(idx_out_ph) = true;
            idx_out_pr = this.phase_obs == 0 & abs(this.res) > pr_thr;
            this.outlier_obs(idx_out_pr) = true;
        end
        
        function reweightHuber(this)
            threshold = 2;
            wfun = @(x) threshold ./ abs(x);
            this.weightOnResidualPh(wfun, threshold);
            this.weightOnResidualPr(wfun, threshold);
        end
        
        function reweightDanish(this)
            threshold = 2;
            wfun = @(x)  exp(-x.^2 ./threshold.^2);
            this.weightOnResidualPh(wfun, threshold);
            this.weightOnResidualPr(wfun, threshold);
        end
        
        function reweightDanishWM(this)
            threshold = 2;
            wfun = @(x)  max(0.5,exp(-x.^2 ./threshold.^2));
            this.weightOnResidualPh(wfun, threshold);
            this.weightOnResidualPr(wfun, threshold);
        end
        
        function reweightHubNoThr(this)
            wfun = @(x) 1 ./ abs(x);
            this.weightOnResidualPh(wfun);
            this.weightOnResidualPr(wfun);
        end
        
        function reweightTukey(this)
            threshold = 2;
            wfun = @(x) (1 - (x ./threshold).^2).^2;
             this.weightOnResidualPh(wfun, threshold);
            this.weightOnResidualPr(wfun, threshold);
        end
        
        function weightOnResidualPh(this, wfun, thr)
            if isempty(this.reweight_obs)
                this.reweight_obs = ones(size(this.variance_obs));
            end
            id_ph = this.phase_obs == 0;
            s0 = mean(abs(this.res(id_ph)).*this.reweight_obs(id_ph));
            res_n = this.res/s0;
            
            idx_rw = abs(res_n) > thr & id_ph;
            
            this.reweight_obs(idx_rw) =  wfun(res_n(idx_rw));
            if sum(this.reweight_obs(idx_rw) < 1e-3) > 0 % observation with weight less than 1e-3 are removed from the adjustment otherwise parameters that depend only on them may suffer numerical instability
                this.outlier_obs(this.reweight_obs < 1e-3) = true;
            end
        end
        
       function weightOnResidualPr(this, wfun, thr)
            if isempty(this.reweight_obs)
                this.reweight_obs = ones(size(this.variance_obs));
            end
            id_ph = this.phase_obs > 0;
            s0 = mean(abs(this.res(id_ph)).*this.reweight_obs(id_ph));
            res_n = this.res/s0;
            
            idx_rw = abs(res_n) > thr & id_ph;
            
            this.reweight_obs(idx_rw) =  wfun(res_n(idx_rw));
            if sum(this.reweight_obs(idx_rw) < 1e-3) > 0 % observation with weight less than 1e-3 are removed from the adjustment otherwise parameters that depend only on them may suffer numerical instability
                this.outlier_obs(this.reweight_obs < 1e-3) = true;
            end
        end
    
        
        function snoopGatt(this, ph_thr, pr_thr)
            % simple threshold on residuals
            %
            % this.Snoop(this, ph_thr, pr_thr)
            idx_out_ph = this.phase_obs & abs(this.res) > ph_thr;
            this.outlier_obs(idx_out_ph) = true;
            idx_out_pr = this.phase_obs == 0 & abs(this.res) > pr_thr;
            this.outlier_obs(idx_out_pr) = true;
            for r = 1 : size(this.rec_xyz,1)
                [res_ph] = getPhRes(this, r);
                idx_ko = Core_Utils.snoopGatt(res_ph, ph_thr, ph_thr/2);
                this.setPhFlag(r,idx_ko);
                [res_pr] = getPrRes(this, r);
                idx_ko = Core_Utils.snoopGatt(res_pr, pr_thr, pr_thr/2);
                this.setPrFlag(r,idx_ko);
            end
            
        end
        
        
        function [res_ph, sat, obs_id, res_id] = getPhRes(this, rec_num)
            % Get phase residuals
            %
            % OUPUT
            %   res_ph          matrix of phase residuals
            %   sat             go_id of the satellite
            %   obs_id          id of the array ls.unique_obs_codes indicating the constallation and tracking of the column
            %   id_res          indix of the value in the res array
            %
            % SYNTAX
            %   [res_ph, sat, obs_id, res_id] = this.getPhRes(rec_num)
            if nargin <2
                rec_num = 1;
            end
            idx_rec = this.receiver_obs == rec_num;
            u_stream = unique(1000 * uint32(this.satellite_obs(idx_rec  & this.phase_obs )) + uint32(this.obs_codes_id_obs(idx_rec  & this.phase_obs )));
            n_stream = length(u_stream);
            min_time_res = min(this.ref_time_obs(idx_rec));
            duration = max(this.ref_time_obs(idx_rec)) - min_time_res;
            time_res = (0:this.obs_rate:duration);
            res_ph = nan(length(time_res),n_stream);
            res_id = zeros(length(time_res),n_stream,'uint32');
            sat = nan(1, n_stream);
            obs_id = nan(1,n_stream);
            sat_c = 9999;
            idx_rec = find(idx_rec);
            for i = 1 : n_stream
                obs_id(i) = rem(u_stream(i), 1000);
                sat(i) = floor(u_stream(i) / 1000);
                if sat_c ~=  sat(i)
                    idx_sat = idx_rec(this.satellite_obs(idx_rec) == sat(i));
                    sat_c = sat(i);
                end
                idx_res = idx_sat(this.obs_codes_id_obs(idx_sat) == obs_id(i));
                if any(idx_res)
                    [~,idx_time] = ismember(this.ref_time_obs(idx_res) - min_time_res,time_res);
                    res_ph(idx_time, i) = this.res(idx_res);
                    res_id(idx_time, i) = find(idx_res);
                end
            end
            
        end
        
        function [res_pr, sat, obs_id] = getPrRes(this, rec)
            % get phase residuals
            %
            % SYNTAX:  [res_pr, sat, obs_id] = getPrRes(this)
            if nargin <2
                rec = 1;
            end
            idx_rec = this.receiver_obs == rec;
            u_stream = unique(1000*uint32(this.satellite_obs(idx_rec  & ~this.phase_obs )) + uint32(this.obs_codes_id_obs(idx_rec  & ~this.phase_obs )));
            n_stream = length(u_stream);
            min_time_res = min(this.ref_time_obs(idx_rec));
            duration = max(this.ref_time_obs(idx_rec)) - min_time_res;
            time_res = (0:this.obs_rate:duration);
            res_pr = nan(length(time_res),n_stream);
            sat = nan(1,n_stream);
            obs_id = nan(1,n_stream);
            sat_c = 9999;
                        idx_rec = find(idx_rec);

            for i = 1 : n_stream
                obs_id(i) = rem(u_stream(i) ,1000);
                sat(i) = floor(u_stream(i)/1000);
                if sat_c ~=  sat(i)
                    idx_sat = idx_rec(this.satellite_obs(idx_rec) == sat(i));
                    sat_c = sat(i);
                end
                idx_res = idx_sat(this.obs_codes_id_obs(idx_sat) == obs_id(i));
                if any(idx_res)
                    [~,idx_time] =  ismember(this.ref_time_obs(idx_res) - min_time_res,time_res);
                    res_pr(idx_time,i) = this.res(idx_res);
                end
            end            
        end
        
        function setPhFlag(this,rec,flag)
            % set phase outlier
            %
            % SYNTAX:  setPhFlag(this,rec,flag)
            idx_rec = this.receiver_obs == rec;
            u_stream = unique(1000*uint32(this.satellite_obs(idx_rec  & this.phase_obs )) + uint32(this.obs_codes_id_obs(idx_rec  & this.phase_obs )));
            n_stream = length(u_stream);
            min_time_res = min(this.ref_time_obs(idx_rec));
            duration = max(this.ref_time_obs(idx_rec)) - min_time_res;
            time_res = (0:this.obs_rate:duration);
            sat_c = 9999;
                                    idx_rec = find(idx_rec);

            for i = 1 : n_stream
                obs_id = rem(u_stream(i) ,1000);
                sat = floor(u_stream(i)/1000);
                if sat_c ~=  sat
                    idx_sat = idx_rec(this.satellite_obs(idx_rec) == sat);
                    sat_c = sat;
                end
                idx_res = idx_sat(this.obs_codes_id_obs(idx_sat) == obs_id);
                if any(idx_res)
                    [~,idx_time] =  ismember(this.ref_time_obs(idx_res) - min_time_res,time_res);
                    this.outlier_obs(idx_res) = flag(idx_time,i);
                end
            end
        end
        
        function setPrFlag(this,rec,flag)
            % set phase outlier
            %
            % SYNTAX:  setPhFlag(this,rec,flag)
            idx_rec = this.receiver_obs == rec;
            u_stream = unique(1000*uint32(this.satellite_obs(idx_rec  & ~this.phase_obs )) + uint32(this.obs_codes_id_obs(idx_rec  & ~this.phase_obs )));
            n_stream = length(u_stream);
            min_time_res = min(this.ref_time_obs(idx_rec));
            duration = max(this.ref_time_obs(idx_rec)) - min_time_res;
            time_res = (0:this.obs_rate:duration);
             sat_c = 9999;
                                     idx_rec = find(idx_rec);

            for i = 1 : n_stream
                obs_id = rem(u_stream(i) ,1000);
                sat = floor(u_stream(i)/1000);
                if sat_c ~=  sat
                    idx_sat = idx_rec(this.satellite_obs(idx_rec) == sat);
                    sat_c = sat;
                end
                idx_res = idx_sat(this.obs_codes_id_obs(idx_sat) == obs_id);
                if any(idx_res)
                    [~,idx_time] =  ismember(this.ref_time_obs(idx_res) - min_time_res,time_res);
                    this.outlier_obs(idx_res) = flag(idx_time,i);
                end
            end
        end
        
        
        function setUpSA(this, rec_work, id_sync, flag, param_selction, parametrization)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if strcmpi(flag,'???')
                o_tmp = rec_work.getObsSet('L??');
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
                o_tmp = rec_work.getObsSet('C??');
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
            else
                o_tmp = rec_work.getObsSet(flag);
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
            end
            if nargin < 6
                parametrization = LS_Parametrization();
            end
            this.unique_time = rec_work.time;
            this.bondParamsGenerateIdx(parametrization);
            
        end
        
        function setUpPPP(this, rec_work, sys_list, id_sync, param_selction, parametrization)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if nargin < 4 || isempty(param_selction)
                param_selction = [this.PAR_REC_X;
                    this.PAR_REC_Y;
                    this.PAR_REC_Z;
                    %  this.PAR_REC_PPB;
                    this.PAR_REC_EB;
                    this.PAR_AMB;
                    this.PAR_REC_CLK_PR;
                    this.PAR_REC_CLK_PH;
                    this.PAR_TROPO;
                    this.PAR_TROPO_N;
                    this.PAR_TROPO_E;
                    %                     repmat(this.PAR_TROPO_Z,10,1);
                    this.PAR_IONO
                    ];  %
            end
            if nargin < 5
                parametrization = LS_Parametrization();
            end
            
            this.setUpSA(rec_work, id_sync, '???', param_selction, parametrization);
            this.absValRegularization(this.PAR_IONO, 1e-2);
        end
        
        function setUpIonoFreePPP(this,rec_work, id_sync)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            param_selction = [this.PAR_REC_X;
                this.PAR_REC_Y;
                this.PAR_REC_Z;
               % this.PAR_REC_PPB;
                this.PAR_REC_EB;
                this.PAR_AMB;
                this.PAR_REC_CLK;
%                this.PAR_REC_CLK_PH;
%                 this.PAR_TROPO;
%                 this.PAR_TROPO_N;
%                 this.PAR_TROPO_E
                ];
            for sys = unique(rec_work.system)
            this.addObsEq(rec_work, rec_work.getPrefIonoFree('L',sys), param_selction);
            this.addObsEq(rec_work, rec_work.getPrefIonoFree('C',sys), param_selction);
            end
            
            ls_param = LS_Parametrization();
            this.bondParamsGenerateIdx(ls_param);
        end
        
        function setUpNET(this, sta_list, coo_rate, flag, param_selction, parametrization)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if nargin < 4
                param_selction = [this.PAR_REC_X ;
                    this.PAR_REC_Y;
                    this.PAR_REC_Z;
                    this.PAR_REC_PPB;
                    this.PAR_REC_EB;
                    this.PAR_IONO;
                    %  this.PAR_REC_EB_LIN;
                    this.PAR_AMB;
                    this.PAR_REC_CLK;
                    this.PAR_TROPO;
                    this.PAR_TROPO_N;
                    this.PAR_TROPO_E;
                    this.PAR_SAT_CLK;
                    this.PAR_SAT_PPB;
                    this.PAR_SAT_EBFR;
                    this.PAR_SAT_EB ];
            end
            if nargin < 5
                parametrization = LS_Parametrization();
            end
            % get time common at least to two receiver
            [p_time, id_sync] = Receiver_Work_Space.getSyncTimeExpanded(sta_list, coo_rate);
            id_rem = sum(~isnan(id_sync),2) <= 1;
            p_time.remEpoch(id_rem);
            id_sync(id_rem,:) = [];
            this.unique_time = p_time;
            
            
            % get observation types common to at least two reciver
            [o_codes, id_sync_o] = Receiver_Work_Space.getCommonObsCode(sta_list);
            id_rem_o = find(sum(~isnan(id_sync_o),2) <= 1)';
            %id_sync_o(id_rem_o,:) = nan;
            
            
            % add equations
            for r = 1 : length(sta_list)
                if strcmpi(flag,'???')
                    o_tmp = sta_list(r).work.getObsSet('L??');
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            %this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                        end
                    end
                    if ~o_tmp.isEmpty
                        this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                    end
                    o_tmp = sta_list(r).work.getObsSet('C??');
                    
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            %this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                            
                        end
                    end
                    if ~o_tmp.isEmpty
                        this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                    end
                else
                    o_tmp = sta_list(r).work.getObsSet(flag);
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                            
                        end
                    end
                    this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                end
                
            end
            this.bondParamsGenerateIdx(parametrization);
            this.absValRegularization(this.PAR_IONO, 1e4);
            
        end
        
        function [time_st, time_end] = getTimePar(this, idx)
            % get the parameter time as GPS_Time
            %
            % SYNTAX:
            %  [time_st, time_end] = getTimePar(this,idx)
            if nargin < 2
                idx = 1:size(this.time_par,1);
            end
            time_st = this.time_min.getCopy();
            time_st.addSeconds( double(this.time_par(idx,1)));
            time_end = this.time_min.getCopy();
            time_end.addSeconds( double(this.time_par(idx,2)));
        end
        
        function computeRefTimeObs(this)
            % compute referebnce time fro observations
            %
            % SYNTAX:
            %   this.computeRefTimeObs()
            this.time_min = this.time_obs.minimum();
            rate = this.time_obs.getRate();
            this.obs_rate = rate;
            this.ref_time_obs = this.time_obs.getNominalTime(rate).getRefTime(this.time_min.getMatlabTime);
        end
                
        function s0 = getSigma0Ph(this)
            % Get sigma0 of phase (PPP solution)
            %
            % SYNTAX:
            %  s0 = this.getSigma0Ph()
            s0 = mean(abs(this.res(this.phase_obs & ~ this.outlier_obs > 0)));
        end
        
        function computeAmbJmps(this)
            % determine if a loss of lock tracking happened simultaneusly
            % in all receiver of a satellite or in all satellite of a
            % receiver
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            n_phase = this.unique_obs_codes;
            min_time = this.time_min;
            min_time_mat = min_time.getMatlabTime;
            rate = this.time_obs.getRate();
            n_epoch = round(max(this.ref_time_obs)/rate)+1;
            cur_mast = zeros(n_sat,length(this.unique_obs_codes));
            
            for r  = 1 : n_rec
                % check number of channel per receiver
                n_p = 0;
                idx_o_r = this.receiver_obs == r;
                o_code_r = this.obs_codes_id_obs(idx_o_r);
                time_r = round(this.ref_time_obs(idx_o_r)/rate);
                sat_r = this.satellite_obs(idx_o_r);
                for s = 1 : n_sat
                    id_o_r_s = sat_r == this.unique_sat_goid(s);
                    o_code_r_s = o_code_r(id_o_r_s);
                    for c = 1 : length(this.unique_obs_codes)
                        if this.unique_obs_codes{c}(2) == 'L'
                            if sum(o_code_r_s == c) > 0
                                n_p = n_p + 1;
                            end
                        end
                    end
                end
                amb_mat = zeros(n_epoch,n_p); % initialize the ambiguity matrix
                n_p = 0;
                n_pa = 1;
                for s = 1 : n_sat
                    id_o_r_s = sat_r == this.unique_sat_goid(s);
                    o_code_r_s = o_code_r(id_o_r_s);
                    time_r_s = time_r(id_o_r_s);
                    for c = 1 : length(this.unique_obs_codes)
                        if this.unique_obs_codes{c}(2) == 'L'
                            id_o_r_s_c = o_code_r_s == c;
                            if sum(id_o_r_s_c) > 0
                                if cur_mast(s,c) == 0
                                    cur_mast(s,c)= r;
                                end
                                n_p = n_p + 1;
                                time_r_s_c = time_r_s(id_o_r_s_c);
                                %                                 amb_mat(time_r_s_c+1, n_p) = n_pa;
                                %                                 n_pa = n_pa+1;
                                %                                 css = this.cycle_slips{r}{this.unique_sat_goid(s)}{c};
                                %                                 if ~isempty(css)
                                %                                     for cs = 1 : css.length
                                %                                         amb_mat(css.getNominalTime(rate).getEpoch(cs).getRefTime(min_time_mat)/rate+1:end,n_p) = n_pa;
                                %                                         n_pa = n_pa+1;
                                %                                     end
                                %                                 end
                                %                                     amb_mat(time_s_r_c+1,n_p) = n_pa;
                                %                                     n_pa = n_pa+1;
                                css = unique([this.cycle_slips{cur_mast(s,c)}{this.unique_sat_goid(s)}{c}.getNominalTime(rate).getRefTime(min_time_mat); this.cycle_slips{r}{this.unique_sat_goid(s)}{c}.getNominalTime(rate).getRefTime(min_time_mat)]); % avoid near rank def
                                css = GPS_Time(min_time_mat, css);
                                if ~isempty(css)
                                    for cs = 1 : css.length
                                        cse1 = css.getNominalTime(rate).getEpoch(cs).getRefTime(min_time_mat)/rate;
                                        if cs ~= css.length
                                            cse2 = css.getNominalTime(rate).getEpoch(cs+1).getRefTime(min_time_mat)/rate;
                                        else
                                            cse2 = max(time_r_s_c)+1;
                                        end
                                        time_s_r_cs = time_r_s_c(time_r_s_c >= cse1 & time_r_s_c < cse2);
                                        amb_mat(time_s_r_cs+1,n_p) = n_pa;
                                        n_pa = n_pa+1;
                                    end
                                end
                            end
                        end
                    end
                end
                
                jmps = find(sum([ false(1,n_p); diff(amb_mat)> 0] | amb_mat == 0,2) == n_p & ~(sum(amb_mat,2) == 0));
                if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                    jmps = [1; jmps];
                end
                %if ~isempty(amb_mat) && sum(abs(amb_mat(end,:) )) ~= 0 %<- if last epoch is full start of the arc is not detected
                jmps = [jmps; find(sum(amb_mat~=0,2) > 0,1,'last')+1];
                %end
                this.rec_amb_jmp{r} = min_time.getCopy().addSeconds((jmps-1)*rate);
                this.ls_parametrization.rec_ppb_opt.steps_set = this.rec_amb_jmp;
                
            end
            if n_rec > 1
                for s = 1 : n_sat
                    % check number of channel per receiver
                    n_p = 0;
                    idx_o_s = this.satellite_obs == this.unique_sat_goid(s);
                    o_code_s = this.obs_codes_id_obs(idx_o_s);
                    time_s = round(this.ref_time_obs(idx_o_s)/rate);
                    rec_s = this.receiver_obs(idx_o_s);
                    for r = 1 : n_rec
                        id_o_s_r = rec_s == r;
                        o_code_s_r = o_code_s(id_o_s_r);
                        for c = 1 : length(this.unique_obs_codes)
                            if this.unique_obs_codes{c}(2) == 'L'
                                if sum(o_code_s_r == c) > 0
                                    n_p = n_p + 1;
                                end
                            end
                        end
                    end
                    amb_mat = zeros(n_epoch,n_p); % initialize the ambiguity matrix
                    n_p = 0;
                    n_pa = 1;
                    for r = 1 : n_rec
                        id_o_s_r = rec_s == r;
                        o_code_s_r = o_code_s(id_o_s_r);
                        time_s_r = time_s(id_o_s_r);
                        
                        for c = 1 : length(this.unique_obs_codes)
                            if this.unique_obs_codes{c}(2) == 'L'
                                id_o_s_r_c = o_code_s_r == c;
                                
                                if sum(id_o_s_r_c) > 0
                                    n_p = n_p + 1;
                                    time_s_r_c = time_s_r(id_o_s_r_c);
                                    %                                     amb_mat(time_s_r_c+1,n_p) = n_pa;
                                    %                                     n_pa = n_pa+1;
                                    css = this.cycle_slips{r}{this.unique_sat_goid(s)}{c};
                                    if ~isempty(css)
                                        for cs = 1 : css.length
                                            cse1 = css.getNominalTime(rate).getEpoch(cs).getRefTime(min_time_mat)/rate;
                                            if cs ~= css.length
                                                cse2 = css.getNominalTime(rate).getEpoch(cs+1).getRefTime(min_time_mat)/rate;
                                            else
                                                cse2 = max(time_s_r_c)+1;
                                            end
                                            time_s_r_cs = time_s_r_c(time_s_r_c >= cse1 & time_s_r_c < cse2);
                                            amb_mat(time_s_r_cs+1,n_p) = n_pa;
                                            n_pa = n_pa+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    jmps = find(sum([ false(1,n_p); diff(amb_mat)> 0] | amb_mat == 0,2) == n_p & ~(sum(amb_mat,2) == 0));
                    if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                        jmps = [1; jmps];
                    end
                    %if ~isempty(amb_mat) && sum(abs(amb_mat(end,:) )) ~= 0 %<- if last epoch is full start of the arc is not detected
                    jmps = [jmps; find(sum(amb_mat~=0,2) > 0,1,'last')+1];
                    %end
                    this.sat_amb_jmp{this.unique_sat_goid(s)} = min_time.getCopy().addSeconds((jmps-1)*rate);
                end
                this.ls_parametrization.sat_ppb_opt.steps_set = this.sat_amb_jmp;
            end
        end
    end
    methods (Static)
        function [Caa, a_hat, id_est_amb] = getEstimableAmb(N,B,tol)
            % get a estimable subset of ambiguoty and its variance
            % coavriance matrix
            %
            % SYNTAX:
            %  [Caa,a_hat] = LS_manipulator_new.getEstimableAmb(N,B)
           
            [L,D,P] = ldl(N);
             if nargin < 3
                tol = max(size(N)) * sqrt(eps(norm(diag(D),inf)));
            end
            keep_id = diag(D) > tol;
            
            L_red = L(keep_id,keep_id);
            iL_red = inv(L_red);
            D_red = D(keep_id,keep_id);
            id_est_amb = (P*keep_id)>0;
            P_red = P(id_est_amb,keep_id);
            N_red =  P_red*L_red*D_red* (L_red')*(P_red');
            Caa = (P_red)*iL_red'*inv(D_red)*(iL_red)*P_red';
            B_red = B(id_est_amb);
            a_hat = Caa*B_red;
        end
        
        function removemanually()
            % 

              if false
                                    if sum(this.param_class == this.PAR_REC_EB) > 0
                            for r = 1 : n_rec
                                idx_par = find(this.class_par == this.PAR_REC_EB & this.rec_par == r & ~this.out_par); % one pseudorange bias per reciever
                                if ~isempty(idx_par)
            
                                    % tell which is pseudorange from the electronic bias
                                    idx_par_psrange = false(size(idx_par));
                                    sys_c_par =  zeros(size(idx_par));
                                    for i = 1 : length(idx_par)
                                        if this.obs_codes_id_par(idx_par(i)) > 0
                                            idx_o_c = this.obs_codes_id_par(idx_par(i));
                                        else % note entirely general  -> (if you put an electronic bias common to more constellations this will not work)
                                            ch_s_o_c = this.ch_set{-this.obs_codes_id_par(idx_par(i))};
                                            idx_o_c = ch_s_o_c(1);
                                        end
                                        sys_c_par(i) =  this.unique_obs_codes{idx_o_c}(1);
                                        if this.phase_par(idx_par(i)) == 1
                                            idx_par_psrange(i) = true;
                                        end
                                    end
                                    idx_par_phase = idx_par(~idx_par_psrange);
            
                                    % system of the electronic bias
                                    sys_c_par_phase = sys_c_par(~idx_par_psrange);
                                    sys_c_par_psrange = sys_c_par(idx_par_psrange);
                                    idx_par_psrange =  idx_par(idx_par_psrange);
            
            
                                    if sum(this.param_class == this.PAR_REC_CLK) > 0 | sum(this.param_class == this.PAR_REC_CLK_PR) > 0 | sum(this.param_class == this.PAR_REC_CLK_PH) > 0
                                        if ~isempty(idx_par_psrange) % <--- remove one pseudorange because it is less complicated afterwards
                                            id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                            chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                            idx_idx_par = find(id_obs_par == chosen_id_obs);
                                            idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))]; % <- all bias of the same observation
                                            bnd_ref = this.wl_id_par(idx_par_psrange(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                            sys_c_ref = sys_c_par(idx_idx_par(1));
                                            idx_par_psrange(idx_idx_par) = [];
                                            sys_c_par_psrange(idx_idx_par) = [];
                                        else
                                            id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                            chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                            idx_idx_par = find(id_obs_par == chosen_id_obs);
                                            idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))]; % <- all bias of the same observation
                                            bnd_ref = this.wl_id_par(idx_par_phase(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                            sys_c_ref = sys_c_par(idx_idx_par(1));
                                            idx_par_phase(idx_idx_par) = [];
                                            sys_c_par_phase(idx_idx_par) = [];
                                        end
                                    end
                                    if sum(this.param_class == this.PAR_IONO) > 0 && this.ls_parametrization.iono(2) == LS_Parametrization.SING_REC
                                        if ~isempty(idx_par_psrange)
                                            for sys_c = unique(sys_c_par_psrange)' % <- each system has its own ionpspherese common to the same elctronic bias so they a new rank deficency is introduced
                                                if sys_c == sys_c_ref
                                                    idx_tmp= idx_par_psrange(sys_c_par_psrange == sys_c & bnd_ref ~= this.wl_id_par(idx_par_psrange));
                                                else
                                                    idx_tmp= idx_par_psrange(sys_c_par_psrange == sys_c);
                                                end
                                                if~isempty(idx_tmp)
                                                    id_obs_par  = this.obs_codes_id_par(idx_tmp);
                                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                                    idx_rm = [idx_rm; idx_tmp(idx_idx_par)];
                                                end
                                            end
                                        else
                                            for sys_c = unique(sys_c_par_phase)' % <- each system has its own ionpspherese common to the same elctronic bias so they a new rank deficency is introduced
                                                if sys_c == sys_c_ref
                                                    idx_tmp= idx_par_phase(sys_c_par_phase == sys_c & bnd_ref ~= this.wl_id_par(idx_par_phase));
                                                else
                                                    idx_tmp= idx_par_phase(sys_c_par_phase == sys_c);
                                                end
                                                if~isempty(idx_tmp)
                                                    id_obs_par  = this.obs_codes_id_par(idx_tmp);
                                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                                    idx_rm = [idx_rm; idx_tmp(idx_idx_par)];
                                                end
                                            end
                                        end
                                    end
                                    if sum(this.param_class == this.PAR_REC_PPB) > 0 | sum(this.param_class == this.PAR_REC_CLK_PH) > 0
                                        if ~isempty(idx_par_phase)
                                            idx_rm = [idx_rm; uint32(idx_par_phase(1))];% remove one phase to put it as reference
                                        end
                                    end
                                end
                            end
                        end
            
                        if sum(this.param_class == this.PAR_IONO) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                            % for each satellite if there is not ata least one double frequrency receiver observing the satellite remove the iono paramter, the delay is going to be absorbed by the clock
                            for s = this.unique_sat_goid
                                idx_par = find(this.class_par == this.PAR_IONO & this.sat_par == s & ~this.out_par); % two code biases not from the same frequency
                                n_freq = numel(unique(this.wl_id_obs(this.satellite_obs == s)));
                                if n_freq == 1
                                    idx_rm = [idx_rm; idx_par];
                                end
                            end
                        end
            
                        % remove ones bias per satellite
                        if  sum(this.param_class == this.PAR_SAT_EB) > 0
                            for s = this.unique_sat_goid
                                idx_par = find(this.class_par == this.PAR_SAT_EB & this.sat_par == s & ~this.out_par); % two code biases not from the same frequency
                                if any(idx_par)
                                    idx_par_psrange = false(size(idx_par));
                                    for i = 1 : length(idx_par)
                                        if this.phase_par(idx_par(i)) == 1
                                            idx_par_psrange(i) = true;
                                        end
                                    end
                                    idx_par_phase = idx_par(~idx_par_psrange);
                                    idx_par_psrange =  idx_par(idx_par_psrange);
                                    if sum(this.param_class == this.PAR_SAT_EBFR) == 0
                                        if sum(this.param_class == this.PAR_SAT_CLK) > 0 || sum(this.param_class == this.PAR_SAT_CLK_PH) > 0  || sum(this.param_class == this.PAR_SAT_CLK_PR) > 0
                                            if ~isempty(idx_par_psrange) % <--- remove one pseudorange because it is less complicated afterwards
                                                id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                                idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))]; % <- all bias of the same observation
                                                bnd_ref = this.wl_id_par(idx_par_psrange(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                                this.log.addMessage(this.log.indent(sprintf('Pseudorange %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_psrange(idx_idx_par(1)))},s)));
                                                %idx_par_psrange(idx_idx_par) = [];
                                            else
                                                id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                                idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))]; % <- all bias of the same observation
                                                bnd_ref = this.wl_id_par(idx_par_phase(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                                %sys_c_ref = sys_c_par(idx_idx_par(1));
                                                %idx_par_phase(idx_idx_par) = [];
                                                this.log.addMessage(this.log.indent(sprintf('Phase %s choosen as reference for sat %d',bnd_ref,s)));
                                            end
                                        end
                                        if sum(this.param_class == this.PAR_IONO) > 0
                                            if ~isempty(idx_par_psrange)
                                                idx_par_psrange = idx_par_psrange(bnd_ref ~= this.wl_id_par(idx_par_psrange)); %<- remove bias of the same frequency of one thta has been already removed
                                                if ~isempty(idx_par_psrange)
                                                    id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                                    idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))];
                                                    this.log.addMessage(this.log.indent(sprintf('Pseudorange %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_psrange(idx_idx_par(1)))},s)));
                                                end
                                            else
                                                idx_par_phase = idx_par_phase(bnd_ref ~= this.wl_id_par(idx_par_phase));  %<- remove bias of the same frequency of one thta has been already removed
                                                if ~isempty(idx_par_phase)
                                                    id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                                    idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))];
                                                    this.log.addMessage(this.log.indent(sprintf('Phase %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_phase(idx_idx_par(1)))},s)));
                                                end
                                            end
                                        end
                                    else
                                        % remove one bias per frequency
                                        u_wl_sat = unique(this.wl_id_par(idx_par))';
                                        for bnd = u_wl_sat
                                            idx_par_psrange_wl = idx_par_psrange(this.wl_id_par(idx_par_psrange) == bnd);
                                            idx_par_phase_wl = idx_par_phase(this.wl_id_par(idx_par_phase) == bnd);
                                            if ~isempty(idx_par_psrange_wl) % <--- remove one pseudorange because it is less complicated afterwards
                                                id_obs_par  = this.obs_codes_id_par(idx_par_psrange_wl);
                                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                                idx_rm = [idx_rm; uint32(idx_par_psrange_wl(idx_idx_par))]; % <- all bias of the same observation
                                                this.log.addMessage(this.log.indent(sprintf('Pseudorange %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_psrange_wl(idx_idx_par(1)))},s)));
                                            else
                                                id_obs_par  = this.obs_codes_id_par(idx_par_phase_wl);
                                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                                idx_rm = [idx_rm; uint32(idx_par_phase_wl(idx_idx_par))]; % <- all bias of the same observation
                                                this.log.addMessage(this.log.indent(sprintf('Phase %s choosen as reference for sat %d',bnd_ref,s)));
                                            end
                                        end
                                    end
                                    if sum(this.param_class == this.PAR_SAT_PPB) > 0 | sum(this.param_class == this.PAR_SAT_CLK_PH) > 0
                                        if ~isempty(idx_par_phase)
                                            idx_rm = [idx_rm; uint32(idx_par_phase(1))];% remove one phase to put it as reference
                                        end
                                    end
                                end
                            end
                        end
            
                        if sum(this.param_class == this.PAR_SAT_PPB) > 0 && sum(this.param_class == this.PAR_REC_PPB) > 0
                            idx_par = find(this.class_par == this.PAR_REC_PPB &  ~this.out_par); % remove ppb from obe recievr
                            idx_rm = [idx_rm; uint32(idx_par(1))];
            
                        end
            
            
            
                       
            
                        if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB) > 0
                            idx_par = find(this.class_par == this.PAR_REC_EB & ~this.out_par);
                            rec_wl_id = uint64(this.wl_id_par(idx_par)) + 100*uint64(this.rec_par(idx_par));
                            rec_wl_id_rm = uint64(this.wl_id_par(idx_rm)) + 100*uint64(this.rec_par(idx_rm));
                            [conflictig_id] = ismember(rec_wl_id,rec_wl_id_rm); % <- you can not remove from same receiver and frequency this could generate inconsistency
                            idx_par(conflictig_id) = [];
                            if ~isempty(idx_par)
                                idx_rm = [idx_rm; uint32(idx_par(1))];
                                this.log.addMessage(this.log.indent(sprintf('Receiver %d and %s choosen as reference',this.rec_par(idx_par(1)),this.unique_obs_codes{this.obs_codes_id_par(idx_par(1))})));
                            end
                        end
            
            
            
            
            
            
                        % remove one bias per signal from one recievr
                        if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB) > 0
                            for e = 1: length(this.unique_obs_codes)
                                if sum(this.obs_codes_id_par(idx_rm(this.class_par(idx_rm) == this.PAR_REC_EB)) == e) == 0% <- if it has not been previously removed
                                    idx_par = find(this.class_par == this.PAR_REC_EB & this.obs_codes_id_par == e & ~this.out_par);
                                    if ~isempty(idx_par)
                                        idx_rm = [idx_rm; uint32(idx_par(1))];
                                        this.log.addMessage(this.log.indent(sprintf('Receiver %d choosen as reference for obs %s',this.rec_par(idx_par(1)),this.unique_obs_codes{e})));
                                    end
                                end
                            end
                        end
            
                        % remove one linear trend bias per signal from one recievr
                        if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB_LIN) > 0
                            for e = 1: length(this.unique_obs_codes)
                                idx_par = find(this.class_par == this.PAR_REC_EB_LIN & this.obs_codes_id_par == e & ~this.out_par);
                                if ~isempty(idx_par)
                                    idx_rm = [idx_rm; uint32(idx_par(1))];
                                    this.log.addMessage(this.log.indent(sprintf('Receiver %d choosen as reference for obs %s',this.rec_par(idx_par(1)),this.unique_obs_codes{e})));
                                end
                            end
                        end
              end
            % phase only TBD!!
            
            
                        % for each sat and fro each contiguous set of ambiguity remove
                        % one abiguity per set of phase bias
            
                        % remember -> there is the possibility of island of
                        % disconnected tracking that you have not considered
            
                        %rule -> you can not have two full receiver removed at this
                        %stage
                        % you start from the most populus receiver than you move to the
                        % other keeping in mind that after the first no receiver can be
                        % completely removed
                        if true
                            if sum(this.class_par(~this.out_par) == this.PAR_AMB) > 0 && (sum(this.param_class == this.PAR_SAT_CLK_PH) > 0 || sum(this.param_class == this.PAR_REC_CLK_PH) > 0 || sum(this.param_class == this.PAR_REC_CLK) > 0 || sum(this.param_class == this.PAR_REC_CLK) > 0)
                                sat_eb_prmz = this.ls_parametrization.getParametrization(this.PAR_SAT_EB);
                                prmz_1 = sat_eb_prmz(4) == LS_Parametrization.SING_TRACK; % first parametrization
                                prmz_2 = sat_eb_prmz(4) == LS_Parametrization.RULE && sum(this.param_class == this.PAR_SAT_EBFR) > 0; % second parametrization
                                all_ambs = this.A_idx(:,this.param_class == this.PAR_AMB);
                                idx_amb_rm_sat = [];
                                if sum(this.param_class == this.PAR_SAT_CLK_PH) > 0 || sum(this.param_class == this.PAR_SAT_CLK) > 0 %|| sum(this.param_class == this.PAR_SAT_EB) > 0
                                    % find the elecronic bias assoictaed with each ambiguity
                                    idx_ambs = find(this.class_par == this.PAR_AMB);
                                    amb2eb = zeros(size(idx_ambs));
                                    % find to which electronic bias the ambiguity is tied
                                    for e = 1: length(idx_ambs)
                                        idx_obs_sample = find(all_ambs == idx_ambs(e),1,'first');
                                        if prmz_1
                                            ebs_tmp = this.obs_codes_id_par(this.A_idx(idx_obs_sample, this.param_class == this.PAR_SAT_EB));
                                        elseif prmz_2
                                            ebs_tmp = this.wl_id_par(this.A_idx(idx_obs_sample, this.param_class == this.PAR_SAT_EBFR));
                                        else
                                            this.log.addError('This Kind of parametrization is not dealt')
                                        end
                                        amb2eb(e) = ebs_tmp(1);
                                    end
                                    clearvars ebs_tmp
                                    sat_eb_const =  this.ls_parametrization.sat_eb(1) == LS_Parametrization.CONST;
                                    if sat_eb_const
                                        jmps_sat ={};
                                        jmps_sat_el ={}; % have been elimated an ambiguity from the block
                                        % determine all arcs jum
                                        for s = this.unique_sat_goid
                                            if ~this.sat_amb_jmp{s}.isEmpty
                                                jmps = round(this.sat_amb_jmp{s}.getNominalTime(this.obs_rate).getRefTime(this.time_min.getMatlabTime)) +1;
                                                jmps_sat{s} = jmps;
                                                jmps_sat_el{s} = false(length(jmps)-1,1);
                                            end
            
                                        end
                                    end
                                    ebs = unique(amb2eb)';
                                    for eb = ebs
                                        if sum(this.param_class == this.PAR_AMB) > 0 && (sum(this.param_class == this.PAR_SAT_CLK) > 0 || sum(this.param_class == this.PAR_SAT_CLK_PH) > 0)
                                            is_first_complete = false; %flag to know if the recievr has been removed completely
                                            idx_ambs_e = idx_ambs;
                                            idx_ambs_e(amb2eb ~= eb) = [];
                                            rec_sat_mtx = zeros(n_rec,n_sat);
                                            for aa = idx_ambs_e
                                                rec_sat_mtx(this.rec_par(aa),this.sat_par(aa)) = 1;
                                            end
                                            idx_ambs_e = Core_Utils.ordinal2logical(idx_ambs_e,length(this.class_par));
            
                                            [~, rec_preference] = sort(sum(rec_sat_mtx,2),'descend');
                                            for s = this.unique_sat_goid
            
                                                idx_par = idx_ambs_e & this.sat_par == s & ~this.out_par;
                                                idx_par = find(idx_par);
                                                if any(idx_par)
                                                    time_par = this.time_par(idx_par,:);
                                                    rec_par  = this.rec_par(idx_par,:);
                                                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
            
                                                    u_time = unique(time_par);
                                                    u_rec = unique(rec_par);
                                                    amb2arc = 1000*uint32(rec_par) + uint32(obs_codes_id_par);
                                                    u_arc = unique(amb2arc);
                                                    arc2eb = rem(u_arc,1000);
                                                    % sort ambiguity per arc length
                                                    arc_len = time_par(:,2) - time_par(:,1);
                                                    [~,idx_arc_len] = sort(arc_len,'descend');
            
            
                                                    [~,amb2u_arc] = ismember(amb2arc,u_arc);
                                                    u_eb = unique(arc2eb);
            
                                                    amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                                                    rec_amb_mat =  zeros(1,length(u_arc),'uint32');
                                                    for t = 1 : size(time_par,1)
                                                        amb_mat(time_par(t,1)+1:time_par(t,2),amb2u_arc(t)) = idx_par(t);
                                                        rec_amb_mat(amb2u_arc(t)) = floor(amb2arc(t) / 1000);
                                                    end
                                                    if ~sat_eb_const
                                                        jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                                                        if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                                            jmps = [1; jmps];
                                                        end
                                                    else
                                                        if any(~jmps_sat_el{s}) || true
                                                            jmps = jmps_sat{s};
                                                            %jmps_sat_el{s}(:) = false;
                                                        else
                                                            jmps = [1; size(amb_mat,1)];
                                                        end
                                                    end
            
            
                                                    for j = 1 : (length(jmps) -1)
                                                        first_rem = true; % one has to be removed al the times
            
                                                        jmp_s = jmps(j);
                                                        jmp_e = jmps(j+1);
                                                        ambs = amb_mat(jmp_s:min(size(amb_mat,1),jmp_e-1),:);
                                                        if any(any(ambs))
                                                            rr = 1;
                                                            not_found = true;
                                                            while rr <= length(rec_preference) && not_found
                                                                rp = rec_preference(rr);
                                                                if ~(sum(rec_sat_mtx(rec_preference(1),:) == 1) == 0 && sum(rec_sat_mtx(rp,:) == 1) <= 1)
                                                                    ambs_r = ambs(:,rec_amb_mat == rp);
                                                                    if any(any(ambs_r)) && (~sat_eb_const || ~jmps_sat_el{s}(j) ||  first_rem )
                                                                        cont_amb_while = true; % find the longest ambiguiti present
                                                                        ia = 1;
                                                                        while cont_amb_while
                                                                            if any(any(ambs_r == idx_par(idx_arc_len(ia))))
                                                                                idx_poss_amb = idx_par(idx_arc_len(ia));
                                                                                cont_amb_while = false;
                                                                            end
                                                                            ia = ia+1;
                                                                        end
                                                                        idx_amb_rm_sat = [idx_amb_rm_sat; uint32(idx_poss_amb)];
                                                                        if sat_eb_const
                                                                            jmps_sat_el{s}(j) = true;
                                                                        end
                                                                        first_rem = false;
                                                                        if ~sat_eb_const || ~any(~jmps_sat_el{s})
                                                                            rec_sat_mtx(rr,s) == 2; % two means eliminated
                                                                        end
                                                                        not_found = false;
                                                                    end
                                                                end
                                                                rr = rr +1;
                                                            end
                                                        end
                                                    end
            
                                                end
                                                idx_rm = [idx_rm; uint32(idx_amb_rm_sat)];
                                            end
                                        end
                                    end
                                end
            
                                if (sum(this.param_class == this.PAR_REC_CLK) > 0 || sum(this.param_class == this.PAR_REC_CLK_PH) > 0) &&  size(this.rec_xyz,1) < 2
                                    % find the elecronic bias assoictaed with each ambiguity
                                    idx_ambs = find(this.class_par == this.PAR_AMB);
                                    amb2eb = zeros(size(idx_ambs));
                                    % find to which electrinuc bias the ambiguity is tied
                                    for e = 1: length(idx_ambs)
                                        idx_obs_sample = find( all_ambs == idx_ambs(e),1,'first');
                                        amb2eb(e) = unique(this.obs_codes_id_par(this.A_idx(idx_obs_sample,this.param_class == this.PAR_REC_EB)));
                                    end
                                    ebs = unique(amb2eb)';
                                    rec_eb_const =  this.ls_parametrization.rec_eb(1) == LS_Parametrization.CONST;
                                    if (rec_eb_const || sum(this.param_class == this.PAR_REC_PPB) > 0) && sum(this.param_class == this.PAR_AMB) > 0
                                        jmps_rec ={};
                                        jmps_rec_el ={}; % have been elimated an ambiguity from the block
                                        % determine all arcs jum
                                        for r = 1: size(this.rec_xyz,1);
                                            jmps = round(this.rec_amb_jmp{r}.getNominalTime(this.obs_rate).getRefTime(this.time_min.getMatlabTime)) +1;
                                            jmps_rec{r} = jmps;
                                            jmps_rec_el{r} = false(length(jmps)-1,1);
                                        end
                                    end
                                    for eb = ebs
                                        % for each rec and for each contiguos set of ambiguity remove one
                                        if sum(this.param_class == this.PAR_AMB) > 0 && (sum(this.param_class == this.PAR_REC_CLK) > 0 || sum(this.param_class == this.PAR_REC_CLK_PH) > 0)
                                            o_ch = this.obs_codes_id_par(idx_amb_rm_sat);
                                            o_ch(o_ch < 0) = length(this.unique_obs_codes) - o_ch(o_ch < 0); % negative index stand for set  put them positive after the numebr of type of observation
                                            forbidden_arc = 1e6*uint32(this.rec_par(idx_amb_rm_sat)) + 1000*uint32(this.sat_par(idx_amb_rm_sat)) + uint32(this.ls_parametrization.rec_eb(4) == LS_Parametrization.SING_TRACK)*uint32(o_ch); % this ambiguities can not be elemitaing without genrating "tension" in the system  %+ uint32(this.obs_codes_id_par(idx_amb_rm_sat)) if there is a bias per tracking the condition is weaker
                                            for r = 1: size(this.rec_xyz,1);
                                                forbidden_arc_rec = rem(forbidden_arc(floor(forbidden_arc/1e6) == r),1e6);
                                                idx_ambs_e = idx_ambs;
                                                idx_ambs_e(amb2eb ~= eb) = [];
                                                idx_ambs_e = Core_Utils.ordinal2logical(idx_ambs_e,length(this.class_par));
                                                idx_par = idx_ambs_e & this.rec_par ==  r & ~this.out_par;
                                                idx_par = find(idx_par);
                                                if any(idx_par)
            
                                                    time_par = this.time_par(idx_par,:);
                                                    sat_par  = this.sat_par(idx_par,:);
                                                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                                                    u_time = unique(time_par);
                                                    amb2arc_a = 1000*uint32(sat_par) + uint32(obs_codes_id_par);
                                                    u_arc = unique(amb2arc_a);
                                                    arc2eb = rem(u_arc,1000);
                                                    % sort ambiguity per arc length
                                                    arc_len = time_par(:,2) - time_par(:,1);
                                                    [~,idx_arc_len] = sort(arc_len,'descend');
            
                                                    [~,amb2arc] = ismember(amb2arc_a,u_arc);
                                                    amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                                                    for t = 1 : size(time_par,1)
                                                        amb_mat(time_par(t,1)+1:time_par(t,2),amb2arc(t)) = idx_par(t);
                                                    end
                                                    if ~rec_eb_const
                                                        jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                                                        if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                                            jmps = [1; jmps];
                                                        end
                                                    else
                                                        if any(~jmps_rec_el{r})
                                                            jmps = jmps_rec{r};
                                                        else
                                                            jmps = [1; size(amb_mat,1)];
                                                        end
                                                    end
                                                    for j = 1 : (length(jmps) -1)
                                                        first_rem = true; % one has to be removed al the times
                                                        jmp_s = jmps(j);
                                                        jmp_e = jmps(j+1);
                                                        ambs = amb_mat(jmp_s:min(size(amb_mat,1),jmp_e-1),:);
                                                        if any(any(ambs))
                                                            cont_amb_while = true; % find the longest ambiguiti present
                                                            ia = 1;
                                                            u_ambs = unique(noZero(ambs));
                                                            while cont_amb_while
                                                                if any(any(u_ambs == idx_par(idx_arc_len(ia))))
                                                                    id_poss_rm = idx_par(idx_arc_len(ia));
                                                                    cont_amb_while = false;
                                                                end
                                                                ia = ia+1;
                                                            end
                                                            idx_start = sum(u_ambs == id_poss_rm) > 0; % i  case everything has been removed to exit the loop
                                                            while any(u_ambs) && sum(id_poss_rm ==  idx_rm) > 0 | sum(amb2arc_a(find(idx_par == id_poss_rm)) == forbidden_arc_rec) > 0 ...
                                                                    | (this.ls_parametrization.rec_eb(4) ~= LS_Parametrization.SING_TRACK && sum(floor(amb2arc_a(find(idx_par == id_poss_rm))/1000) == floor(forbidden_arc_rec/1000)) > 0)% it might be that the ambiguity was previouly removed in the satellite round
                                                                u_ambs(u_ambs == id_poss_rm) = 0;
                                                                cont_amb_while = true;  % find the longest ambiguiti present
                                                                while cont_amb_while & any(u_ambs)
                                                                    if any(u_ambs == idx_par(idx_arc_len(ia)))
                                                                        id_poss_rm = idx_par(idx_arc_len(ia));
                                                                        cont_amb_while = false;
                                                                    end
                                                                    ia = ia+1;
                                                                end
                                                            end
                                                            if id_poss_rm > 0 && (~rec_eb_const || ~jmps_rec_el{r}(j) ||  first_rem )
                                                                idx_rm = [idx_rm; uint32(id_poss_rm)];
                                                                if rec_eb_const
                                                                    jmps_rec_el{r}(j) = true;
                                                                end
                                                                first_rem = false;
            
                                                            end
            
            
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        
                            
                       
                            end
                        
            
        end
    end
end
