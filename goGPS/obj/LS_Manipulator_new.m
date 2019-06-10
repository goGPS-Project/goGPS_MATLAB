%   CLASS LS_Manipulator_new
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
%    |___/                    v 1.0 beta 3jp
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:
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
        PAR_SAT_CLK = 11;
        PAR_ANT_MP = 12;
        PAR_IONO = 13;
        PAR_TROPO_S = 14;
        PAR_SAT_X = 15;
        PAR_SAT_Y = 16;
        PAR_SAT_Z = 17;
        PAR_SAT_EB = 18;
        PAR_REC_EB_LIN = 19;
        
        
        CLASS_NAME = {'REC_X', 'REC_Y', 'REC_Z', 'REC_EB', 'AMB', 'REC_CLK', 'TROPO', 'TROPO_N', 'TROPO_E', 'SAT_CLK', 'ANT_MP', 'IONO', 'TROPO_S', 'SAT_X', 'SAT_Y', 'SAT_Z', 'SAT_EB', 'REC_EB_LIN'};
    end
    
    properties
        %%% sparse A , this in more efficient than a sparse(slow to index in matlab) because
        %%% normally all the obervation euation have the same numebr of
        %%% entry
        A % obervation
        A_idx
        obs
        param_class % class id of the colum of A
        time_obs   % epoch of the obervation (GPS_Time)
        satellite_obs % goid satellite of the observation
        receiver_obs % reciver of the obeservation
        azimuth_obs % azimuth of the observation
        elevation_obs % elevation of the observation
        variance_obs % varaince of the observations
        obs_codes_id_obs % id of the signal used
        phase_obs % logical to tle  if obs are phase or code
        wl_id_obs % id of the wavelength
        
        A_pseudo
        A_idx_pseudo
        time_pseudo   % epoch of the pseudo-observation (GPS_Time)
        satellite_pseudo % goid satellite of the pseudo-observation
        receiver_pseudo % reciver of the pseudo-observation
        variance_pseudo % varaince of the pseudo-observation
        
        unique_obs_codes % uniques ids (cell) of the signals  /since lot of combinations are possible they will be generated dynamically)
        unique_wl % set of unique wavelength
        rec_xyz % receiver coordinates to be used in
        unique_rec_name % names of the recivers
        unique_sat_goid % unique satellite goids
        cycle_slips = cell(1) % epoch of the cycle slip
        
        time_par   % time of the paramter !!!! very important the parameters (within the same class e.g. time for satellite s ) MUST be ordered in chronological order
        param_par  % paramtrization of the paramter
        time_min   % ref_time_of the parameter
        rec_par    % receiver of the paramters
        sat_par    % staellite of the paramters
        class_par  % class of the paramter
        obs_codes_id_par  % obs code id fo the paramter
        wl_id_par % wl id of the paramter
        
        rec_set % set of receivers
        sat_set % set of satellites
        ch_set  % set of observation codes
        
        N
        idx_rd % idx paramter removed to silve the rank deficency
        ls_parametrization;
        
        
        log
    end
    
    methods
        function addObsEq(this,rec,obs_set, param_selection)
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
            
            % --- check all paramters presence and their order -----
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
            par_rec_eb_lin_lid = param_selection == this.PAR_REC_EB_LIN;
            par_rec_eb_lin = sum(par_rec_eb_lin_lid) > 0;
            par_sat_eb_lid = param_selection == this.PAR_REC_EB;
            par_sat_eb = sum(par_sat_eb_lid) > 0;
            
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
            
            par_iono_lid = param_selection == this.PAR_IONO;
            par_iono = sum(par_iono_lid) > 0;
            
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
            
            
            mfw_on =  sum(param_selection == this.PAR_TROPO | param_selection == this.PAR_TROPO_E | param_selection == this.PAR_TROPO_N | param_selection == this.PAR_TROPO_V) > 0;
            % get the mapping function for tropo
            if mfw_on
                id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
            % check whivh observations are phase ones
            phase_s = obs_set.obs_code(:,2) == 'L'; % first char is the system second MUST be the bervation type (Phase pseudorance doppler snr)
            
            % initliaze the matrices
            A = zeros(n_obs, n_par);
            [obs,satellite_obs, azimuth_obs, elevation_obs, variance_obs, wl_obs] = deal(zeros(n_obs, 1));
            [ obs_codes_id_obs,  wl_obs] = deal(zeros(n_obs, 1,'uint8'));
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
            iono_const = GPS_SS.L_VEC(1)^2;
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
                    azimuth_obs(lines_stream) =  obs_set.sigma(s)^2;
                    elevation_obs(lines_stream) =  obs_set.sigma(s)^2;
                    obs_codes_id_obs(lines_stream) = s_s_id;
                    phase_obs(lines_stream) = phase_s(s);
                    wl_obs(lines_stream) = wl_id;
                    time_obs.addEpoch(obs_set.time.getEpoch(id_ok_stream).getMatlabTime);
                    
                    % ----------- save the cycle slips
                    if phase_s(s)
                        %obs_set.cycle_slip(,s) = true;
                        if any(obs_set.obs(:,s)) && ~any(obs_set.cycle_slip(:,s))
                            obs_set.cycle_slip(1,s) = true;
                        end
                        this.cycle_slips{r}{s_go_id}{s_s_id} = obs_set.time.getEpoch(find(obs_set.cycle_slip(:,s) & obs_set.obs(:,s)~=0));
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
                        A(lines_stream, par_rec_x_lid) =  los_stream(:,1);
                    end
                    
                    if par_sat_y
                        A(lines_stream, par_rec_y_lid) =  los_stream(:,2);
                    end
                    
                    if par_sat_z
                        A(lines_stream, par_rec_z_lid) =  los_stream(:,3);
                    end
                    % ----------- electronic bias ------------------
                    if par_rec_eb
                        A(lines_stream, par_rec_eb_lid) = 1;
                    end
                    if par_rec_eb_lin
                        A(lines_stream, par_rec_eb_lid) = 1/obs_set.wl(s);
                    end
                    if par_sat_eb
                        A(lines_stream, par_sat_eb_lid) = 1;
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
                    % ----------- ZTD ------------------
                    if par_tropo
                        A(lines_stream, par_tropo_lid) = mfw_stream;
                    end
                    % ----------- ZTD gradients ------------------
                    if par_tropo_n || par_tropo_e
                        cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                        if par_tropo_e
                            A(lines_stream, par_tropo_e_lid) = sin(az_stream) .* cotan_term; % east gradient
                        end
                        if par_tropo_n
                            A(lines_stream, par_tropo_n_lid) = cos(az_stream) .* cotan_term; % noth gradient
                        end
                    end
                    if par_tropo_v
                        A(lines_stream, par_tropo_v_lid) = mfw_stream*rec.h_ellips;
                    end
                    % ----------- Ionosphere delay --------------------
                    if par_iono
                        if phase_s(s)
                            A(lines_stream, par_iono_lid) = - obs_set.wl(s)^2/iono_const;
                        else
                            A(lines_stream, par_iono_lid) =   obs_set.wl(s)^2/iono_const;
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
            % bond paramters (splines or other models) and generate idx
            %
            % SYNTAX
            %    this.bondParamGenerateIdx(parametrization)
            this.log = Core.getLogger();
            this.ls_parametrization = ls_parametrization;
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            this.A_idx = zeros(size(this.A));
            time_min = min(this.time_obs.getMatlabTime);
            this.time_min = this.time_obs.minimum();
            obs_rate = this.time_obs.getRate; % TBD substitute this quantity with the obesravtion minimum rate
            time_obs = round(this.time_obs.getRefTime(time_min)/obs_rate);
            
            cumulative_idx = 0;
            i_col = 1;
            
            % ---- preallocation to speed up
            ch_lid = false(size(this.A,1),1);
            obs_lid = false(size(this.A,1),1);
            
            
            this.time_par = zeros(size(this.A,1),2,'uint32');
            this.param_par = zeros(size(this.A,1),4,'uint8'); % time of the paramter
            this.rec_par = zeros(size(this.A,1),1,'uint8'); % receiver of the paramters
            this.sat_par = zeros(size(this.A,1),1,'uint8');  % receiver of the paramters
            this.class_par = zeros(size(this.A,1),1,'uint8');  % class of the paramter
            this.obs_codes_id_par= zeros(size(this.A,1),1,'uint8');  % obs_code id paramters
            this.wl_id_par= zeros(size(this.A,1),1,'uint8');  % obs_code id paramters
            
            ch_set_old = []; % avoid repating expesnive task 
            
            for i_p = 1 : length(this.param_class)
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
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        ant_id = Core.getConstellationCollector.getAntennaId(this.satellite_obs(idx1));
                        sig2const(c) = ant_id(1);
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
                elseif parametriz(4) == ls_parametrization.RULE
                    % ---- evaluate the rule
                    n_ch_set = 0;
                    ch_set = {};
                    for rule = opt.rule
                        cond_sig = true(size(sig_p_id));
                        parts = strsplit(rule{1},':');
                        cond = parts{1};
                        opt = parts{2};
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
                            if str2num(opt) == ls_parametrization.ALL_FREQ
                                n_ch_set = n_ch_set +1;
                                ch_set{n_ch_set} = uint8(sig_p_id(cond_sig));
                            elseif str2num(opt) == ls_parametrization.SING_FREQ
                                u_wl_tmp = unique(sig2wl(cond_sig));
                                n_wl_tmp = length(u_wl_tmp);
                                for c = 1 : n_wl_tmp
                                    n_ch_set = n_ch_set +1;
                                    ch_set{n_ch_set} = uint8(sig_p_id(sig2wl == u_wl_tmp(c) & cond_sig));
                                end
                            elseif str2num(opt) == ls_parametrization.SING_TRACK
                                n_ch_set = n_ch_set +length(sig_p_id(cond_sig));
                                ch_set = [ch_set num2cell(uint8(sig_p_id(cond_sig)))];
                                
                            end
                        end
                        
                    end
                end
                
                % ------- Now constructing the index, the epoch dependednce will be dealt inside the loop -------------
                for r = 1 : n_rec_set
                    rec_lid = false(size(this.A,1),1);
                    for rr = rec_set{r}
                        rec_lid = rec_lid | this.receiver_obs == rr;
                    end
                    % find an id for the reciver set to keep track of the
                    % paramters if recievr is sigle is simply the reciver
                    % progessince number, otherwise they are negative
                    % number with the index in the receiver set
                    if length(rec_set{r}) == 1
                        r_id = rec_set{r};
                    else
                        r_id =Core_Utils.findAinB(rec_set{r},this.rec_set);
                        if r_id == 0
                            r_id = length(this.rec_set) -1;
                            this.rec_set{end+1} = rec_set{r};
                        end
                        
                    end
                    for s = 1 : n_sat_set
%                         sat_lid = false(size(this.A,1),1);
%                         for ss = sat_set{s}
%                             sat_lid = sat_lid | this.satellite_obs == ss;
%                         end
                        sat_lid = ismember(this.satellite_obs,sat_set{s});
                        ss = sat_set{end};
                        if length(sat_set{s}) == 1
                            s_id = sat_set{s};
                        else
                            s_id =Core_Utils.findAinB(sat_set{s},this.sat_set);
                            if s_id == 0
                                s_id = length(this.sat_set) -1;
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
                                    ch_lid = ismember(this.obs_codes_id_obs,ch_set{f});
                                    ch_set_old = ch_set{f};
                                end
                                cc = ch_set{f};
                                if length(ch_set{f}) == 1
                                    ch_id = ch_set{f};
                                else
                                    ch_id =Core_Utils.findAinB(ch_set{f},this.sat_set);
                                    if s_id == 0
                                        ch_id = length(this.ch_set) -1;
                                        this.ch_set{end+1} = ch_set{f};
                                    end
                                    
                                end
                                %---- final observation id
                                obs_lid = rec_lid & sat_lid & ch_lid;
                                if any(obs_lid)
                                    % --- now dealing with epoch dependence ----
                                    cols_tmp = 0;
                                    if parametriz(1) == ls_parametrization.CONST
                                        ep_pgr_id = 1;
                                        n_prg_id = 1;
                                        time_par_tmp = [this.time_obs.getEpoch(obs_lid).minimum.getRefTime(this.time_min.getMatlabTime)  this.time_obs.getEpoch(obs_lid).maximum.getRefTime(this.time_min.getMatlabTime)];
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
                                                else
                                                    ep_pgr_id = 1;
                                                end
                                            end
                                        else
                                            steps_set = opt.steps_set;
                                            
                                            if parametriz(2) == ls_parametrization.ALL_REC % you can use differents step for step wise satellite dependent paraters
                                                steps = round(steps_set{ss}.getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    ep_pgr_id(lid_maj) = p_s;
                                                    time_par_tmp = [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]; %start of the arc
                                                    if p_s > 1
                                                        time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                    end
                                                    p_s = p_s +1;
                                                end
                                            elseif parametriz(3) == ls_parametrization.ALL_SAT  % you can use differents step for step wise receiver dependent paraters
                                                steps = round(steps_set{rr}.getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    ep_pgr_id(lid_maj) = p_s;
                                                    time_par_tmp = [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]; %start of the arc
                                                    if p_s > 1
                                                        time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                    end
                                                    p_s = p_s +1;
                                                end
                                            end
                                        end
                                        n_prg_id = p_s - 1;
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
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate),1);
                                        u_e_tmp = unique([ep_id ep_id+1]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        this.A = [this.A(:,1:(i_col-1)) this.A(:,i_col).*spline_v(:,1) this.A(:,i_col).*spline_v(:,2) this.A(:,(i_col+1):end)];
                                        n_prg_id = length(u_e_tmp);
                                    elseif parametriz(1) == ls_parametrization.SPLINE_CUB
                                        % ---- colum will be quadrupled -----
                                        cols_tmp = [ 0 1 2 3];
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate),1);
                                        u_e_tmp = unique([ep_id ep_id+1 ep_id+2 ep_id+3]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        this.A = [this.A(:,1:(i_col-1)) this.A(:,i_col).*spline_v(:,1) this.A(:,i_col).*spline_v(:,2) this.A(:,i_col).*spline_v(:,3) this.A(:,i_col).*spline_v(:,4) this.A(:,(i_col+1):end)];
                                        n_prg_id = length(u_e_tmp);
                                    end
                                    this.A_idx(obs_lid,i_col + cols_tmp) = cumulative_idx + ep_pgr_id;
                                    [u_new_par] = unique(cumulative_idx + ep_pgr_id(:));
                                    
                                    this.time_par(cumulative_idx+(1:n_prg_id),:) =  uint32(time_par_tmp);% time of the paramter
                                    this.param_par(cumulative_idx+(1:n_prg_id),:) = repmat(uint8(parametriz),n_prg_id,1);% time of the paramter
                                    this.rec_par(cumulative_idx+(1:n_prg_id)) = r_id*ones(n_prg_id,1,'uint8');  % receiver of the paramters
                                    this.sat_par(cumulative_idx+(1:n_prg_id)) = s_id*ones(n_prg_id,1,'uint8');  % receiver of the paramters
                                    this.class_par(cumulative_idx+(1:n_prg_id)) =  p*ones(n_prg_id,1,'uint8');  % class of the paramter
                                    this.obs_codes_id_par(cumulative_idx+(1:n_prg_id)) = uint8(ch_id)*ones(n_prg_id,1,'uint8');  % obs_code id paramters
                                    
                                    if ch_id > 0
                                        this.wl_id_par(cumulative_idx+(1:n_prg_id)) = uint8(sig2wl(ch_id))*ones(n_prg_id,1,'uint8');  % wavelength id paramters
                                    else
                                        this.wl_id_par(cumulative_idx+(1:n_prg_id)) = zeros(n_prg_id,1,'uint8');  % wavelength id paramters
                                    end
                                    cumulative_idx = cumulative_idx + n_prg_id;
                                    col_incr = max(col_incr,length(cols_tmp));
                                end
                            end
                        end
                    end
                end
                i_col = i_col + col_incr;
            end
            % free allocated space in excess
            this.time_par((cumulative_idx+1):end,:) =  [];% 
            this.param_par((cumulative_idx+1):end,:) =  [];% 
            this.rec_par((cumulative_idx+1):end) =  [];% 
            this.sat_par((cumulative_idx+1):end) =  [];% 
            this.class_par((cumulative_idx+1):end) =  [];% 
            this.obs_codes_id_par((cumulative_idx+1):end) =  [];% 
              this.wl_id_par((cumulative_idx+1):end) =  [];% 
            
            
            
        end
        
        function removeFullRankDeficency(this)
            % solve full rank deficency removing paramters from the
            % estimation
            %
            % SYNTAX:
            %    this.removeFullRankDeficency()
            
            % remove two bias per receiver
            n_rec = size(this.rec_xyz,1);
            idx_rm = [];
            for r = 1 : n_rec 
                idx_par = find(this.class_par == this.PAR_REC_EB & this.rec_par == r); % one for phase one for code
                idx_par_psrange = false(size(idx_par));
                for i = 1 : length(idx_par)
                    if this.unique_obs_codes{this.obs_codes_id_par(idx_par(i))}(2) == 'C'
                        idx_par_psrange(i) = true;
                    end
                end
                idx_par_phase = idx_par(~idx_par_psrange);
                idx_par_psrange =  idx_par(idx_par_psrange);
                if ~isempty(idx_par_psrange)
                    idx_rm = [idx_rm; idx_par_psrange(1)];
                end
                if ~isempty(idx_par_phase)
                    idx_rm = [idx_rm; idx_par_phase(1)];
                end
            end
            % remove two bias per satellite
            for s = this.unique_sat_goid
                idx_par = find(this.class_par == this.PAR_SAT_EB & this.sat_par == s); % one for phase one for code
                 idx_par_psrange = false(size(idx_par));
                for i = 1 : length(idx_par)
                    if this.unique_obs_codes{this.obs_codes_id_par(idx_par(i))}(2) == 'C'
                        idx_par_psrange(i) = true;
                    end
                end
                idx_par_phase = idx_par(~idx_par_psrange);
                idx_par_psrange =  idx_par(idx_par_psrange);
                if ~isempty(idx_par_psrange)
                    idx_rm = [idx_rm; idx_par_psrange(1)];
                end
                if ~isempty(idx_par_phase)
                    idx_rm = [idx_rm; idx_par_phase(1)];
                end
            end
            
            % remove one clock per epoch
            u_ep = unique(this.time_par);
            if sum(this.param_class == this.PAR_REC_CLK) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                for e = u_ep'
                    idx_par = find(this.class_par == this.PAR_REC_CLK & this.time_par(:,1) == e);
                    if ~isempty(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                    end
                end
            end
            
            % for each sat and fro each contiguos set of ambiguity remove
            % one abiguity per set of phase bias
            if sum(this.param_class == this.PAR_AMB) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                for s = this.unique_sat_goid
                    idx_par = this.class_par == this.PAR_AMB & this.sat_par == s;
                    idx_par = find(idx_par);
                    time_par = this.time_par(idx_par,:);
                    rec_par  = this.rec_par(idx_par,:);
                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                    eb_id_par = zeros(size(idx_par));
                    % find to which electrinuc bias the ambiguity is tied
                    for e = 1: length(idx_par)
                        idx_obs_sample = find(this.A_idx(:,this.param_class == this.PAR_AMB) == idx_par(e),1,'first');
                        eb_id_par(e) = this.A_idx(idx_obs_sample,this.param_class == this.PAR_REC_EB);
                    end
                    u_time = unique(time_par);
                    u_rec = unique(rec_par);
                    u_obs_codes = unique(obs_codes_id_par);
                    arc = 1000*uint32(rec_par) + uint32(obs_codes_id_par);
                    u_arc = unique(arc);
                    [~,arc] = ismember(arc,u_arc);
                    eb_arc = eb_id_par(unique(arc));
                    u_eb = unique(eb_arc);
                    eb_arc_rem = false(size(u_eb));
                    
                    amb_mat = zeros(max(u_time),length(u_arc),'uint8');
                    for t = 1 : size(time_par,1)
                        amb_mat(time_par(t,1)+1:time_par(t,2),arc(t)) = idx_par(t);
                    end
                    jmps = [1; find(diff(sum(amb_mat,2) > 0) == 1); max(u_time)];
                    for j = 1 : (length(jmps) -1)
                        jmp_s = jmps(j);
                        jmp_e = jmps(j+1);
                        while sum(eb_arc_rem) < length(eb_arc_rem) % one ambiguity per electronic bias associated with the ambiguity
                            [~,   idx_ambs] = ismember(eb_arc,u_eb(~eb_arc_rem));
                            ambs = amb_mat(jmp_s:jmp_e,idx_ambs>0);
                            eb_arc_tmp = eb_arc(idx_ambs>0);
                            if any(ambs)
                                idx_poss_amb = mode(noZero(ambs(:)));
                                idx_rm = [idx_rm; idx_poss_amb];
                                eb_arc_rem(u_eb == eb_arc_tmp(sum(ambs == idx_poss_amb) > 0)) = true;
                            else
                                eb_arc_rem = true;
                            end
                            
                        end
                    end
                    
                    
                end
            end
            
            % for each rec and for each contiguos set of ambiguity remove one
            if sum(this.param_class == this.PAR_AMB) > 0 && sum(this.param_class == this.PAR_REC_CLK) > 0 
                for r = 1: size(this.rec_xyz,1);
                    idx_par = this.class_par == this.PAR_AMB & this.rec_par == r;
                    idx_par = find(idx_par);
                    time_par = this.time_par(idx_par,:);
                    sat_par  = this.sat_par(idx_par,:);
                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                    eb_id_par = zeros(size(idx_par));
                    % find to which electrinuc bias the ambiguity is tied
                    for e = 1: length(idx_par)
                        idx_obs_sample = find(this.A_idx(:,this.param_class == this.PAR_AMB) == idx_par(e),1,'first');
                        eb_id_par(e) = this.A_idx(idx_obs_sample,this.param_class == this.PAR_REC_EB);
                    end
                    
                    u_time = unique(time_par);
                    u_sat = unique(sat_par);
                    u_obs_codes = unique(obs_codes_id_par);
                    arc = 1000*uint32(sat_par) + uint32(obs_codes_id_par);
                    u_arc = unique(arc);
                    [~,arc] = ismember(arc,u_arc);
                    eb_arc = eb_id_par(unique(arc));
                    u_eb = unique(eb_arc);
                    eb_arc_rem = false(size(u_eb));
                    
                    amb_mat = zeros(max(u_time),length(u_arc),'uint8');
                    for t = 1 : size(time_par,1)
                        amb_mat(time_par(t,1)+1:time_par(t,2),arc(t)) = idx_par(t);
                    end
                    jmps = [1; find(diff(sum(amb_mat,2) > 0) == 1); max(u_time)];
                    for j = 1 : (length(jmps) -1)
                        jmp_s = jmps(j);
                        jmp_e = jmps(j+1);
                        while sum(eb_arc_rem) < length(eb_arc_rem) % one ambiguity per electronic bias associated with the ambiguity
                            [~,   idx_ambs] = ismember(eb_arc,u_eb(~eb_arc_rem));
                            ambs = amb_mat(jmp_s:jmp_e,idx_ambs>0);
                            eb_arc_tmp = eb_arc(idx_ambs>0);
                            if any(ambs)
                                id_poss_rm = mode(noZero(ambs(:)));
                                while sum(id_poss_rm ==  idx_rm) > 0 % it might be that the ambiguity was previouly removed in the satellite round
                                    ambs(ambs == id_poss_rm) = 0;
                                    id_poss_rm = mode(noZero(ambs(:)));
                                end
                                idx_rm = [idx_rm; id_poss_rm];
                                eb_arc_rem(u_eb == eb_arc_tmp(sum(ambs == id_poss_rm) > 0)) = true;
                            else
                                eb_arc_rem = true;
                            end
                            
                        end
                    end
                end
            end
            
            
            % ---- (multi receiver) for each epoche remove one coordinate --------
            if (sum(this.param_class == this.PAR_REC_X) > 0 || sum(this.param_class == this.PAR_REC_Y) > 0  ||  sum(this.param_class == this.PAR_REC_Z) > 0 ) && sum(this.param_class == this.PAR_SAT_CLK) > 0
                for e = u_ep'
                    idx_e = this.time_par(:,1) >= e & ( this.time_par(:,1)==0 | this.time_par(:,2) < e) ;
                    idx_par = this.class_par == this.PAR_REC_X & idx_e;
                    idx_par = find(idx_par);
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    idx_par = this.class_par == this.PAR_REC_Y & idx_e;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    idx_par = this.class_par == this.PAR_REC_Z & idx_e;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                end
            end
            
            this.idx_rd = sort(noZero(idx_rm));
            %             for i = length(this.idx_rd) : -1 : 1
            %                 ir = this.idx_rd(i);
            %                 this.A_idx(this.A_idx == ir) = 0;
            %                 i_mj = this.A_idx > ir;
            %                 this.A_idx(i_mj) = this.A_idx(i_mj) -1;
            %             end
        end
        
        
        
        
        function absValRegularization(this,p_class, var)
            % regularize paramters to zero (Tykhnov aka ridge aka L2)
            %
            % this.absValRegularization(param_id, var)
            par_ids = this.parameter_class == p_class;
            u_p_id = unique(this.A_idx(:,par_ids));
            n_par = length(u_p_id);
            [A_tmp, A_idx_tmp] = deal(zeros(n_par, 2)); % tykhonv regualrization are now limited to two paramters
            A_tmp(:,1) = 1;
            A_idx_tmp(:,1) = u_p_id;
            this.A_pseudo = [this.A_pseudo; A_tmp];
            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
            this.variance_pseudo = [this.variance_pseudo; var*ones(n_par,1)];
            this.receiver_pseudo = [this.receiver_pseudo; rec_par(u_p_id)];
            this.time_pseudo.addEpoch(GPS_Time(this.time_min + this.time_par(u_p_id)));
            this.satellite_pseudo = [this.satellite_pseudo; sat_par(u_p_id)];
        end
        
        function timeRegularization(this, p_class, var_per_sec)
            % first order tykhonv regualrization in time
            %
            % this.timeRegularization(this, param_id, var)
            par_ids = this.parameter_class == p_class;
            p_idx = unique(this.A_idx(:,par_ids));
            rec_idx = this.rec_par(p_idx);
            u_rec = unique(rec_idx);
            sat_idx = this.sat_par(p_idx);
            u_sat = unique(sat_idx);
            ch_idx  =  this.ch_par(p_idx);
            u_ch = unique(ch_idx);
            for r = u_rec
                for s = u_sat
                    for c = u_ch
                        p_tmp = p_idx(rec_idx == r & sta_idx == s & ch_idx == c); % find the idx of the obsrevations
                        if ~isempty(p_tmp)
                            u_p_id = unique(p_tmp);
                            n_par = length(u_p_id);
                            A_tmp = [ones(n_par-1,1) -ones(n_par-1,1)];
                            A_idx_tmp = [u_p_id(1:end-1) u_p_id(2,end)];
                            this.A_pseudo = [this.A_pseudo; A_tmp];
                            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
                            
                            this.variance_pseudo = [this.variance_pseudo; var_per_sec * diff(this.time_par(u_p_id))];
                            % taking the indices of first epoch
                            this.receiver_pseudo = [this.receiver_pseudo; rec_par(u_p_id(1:end-1))];
                            this.time_pseudo.addEpoch(GPS_Time(this.time_min + this.time_par(u_p_id(1:end-1))));
                            this.satellite_pseudo = [this.satellite_pseudo; sat_par(u_p_id(1:end-1))];
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
            % reduce for the paramter
        end
        
        
        function solve(this)
            % ------ form the normal matrix
            n_obs = size(this.A,1) + size(this.A_pseudo,1);
            n_par = max(max(this.A_idx));
            rows = repmat((1:size(this.A,1))',1,size(this.A,2));
            rows_pseudo = repmat((1:size(this.A_pseudo,1))',1,size(this.A_pseudo,2));
            A = sparse([rows(:); rows_pseudo(:)],zero2n([this.A_idx(:) this.A_idx_pseudo(:)],1),[this.A(:) this.A_pseudo(:)],n_obs,n_par);
            A(:,this.idx_rd) = [];
            class_par = this.class_par;
            class_par(this.idx_rd) = [];
            vars = [1./this.variance_obs; 1./this.variance_pseudo];
            mean_vars = mean(vars);
            vars = vars ./ mean_vars;
            Cyy =  spdiags(vars,0,n_obs,n_obs);
            x_est = zeros(n_par -length(this.idx_rd),1);
            Aw = A'*Cyy;
            N = Aw*A;
            y = sparse([this.obs; zeros(size(this.A_pseudo,1),1)]);
            B = Aw*y;
            
            clearvars Aw
            % ------ reduce for sat clock, rec clock and iono
            idx_reduce_sat_clk = class_par == this.PAR_SAT_CLK;
            idx_reduce_rec_clk = class_par == this.PAR_REC_CLK;
            idx_reduce_iono = class_par == this.PAR_IONO;
            
            iono = sum(this.param_class  == this.PAR_IONO) > 0;
            if iono
                n_iono = sum(idx_reduce_iono);
                iIono = spdiags(1./diag(N(idx_reduce_iono,idx_reduce_iono)),0,n_iono,n_iono);
                Nx_iono = N(~idx_reduce_iono,idx_reduce_iono); % cross term reduce iono
                Nt = Nx_iono * iIono;
                N = N(~idx_reduce_iono,~idx_reduce_iono) - Nt * N(idx_reduce_iono,~idx_reduce_iono);
                B_iono =  B(idx_reduce_iono);
                B = B(~idx_reduce_iono) - Nt * B_iono;
            end
            
            sat_clk = sum(this.param_class  == this.PAR_SAT_CLK) > 0;
            i_sat_clk_tmp = idx_reduce_sat_clk(~idx_reduce_iono);
            if sat_clk
                n_clk_sat = sum(i_sat_clk_tmp);
                iSatClk = spdiags(1./diag(N(i_sat_clk_tmp,i_sat_clk_tmp)),0,n_clk_sat,n_clk_sat);
                Nx_satclk = N(~i_sat_clk_tmp, i_sat_clk_tmp);
                Nt = Nx_satclk * iSatClk;
                N = N(~i_sat_clk_tmp,~i_sat_clk_tmp) - Nt * N(i_sat_clk_tmp, ~i_sat_clk_tmp);
                B_satclk =  B(i_sat_clk_tmp);
                B = B(~i_sat_clk_tmp) - Nt * B_satclk;
            end
            
            i_rec_clk_tmp = idx_reduce_rec_clk(~i_sat_clk_tmp);
            iRecClk = inv(N(i_rec_clk_tmp,i_rec_clk_tmp));
            Nx_recclk = N(~i_rec_clk_tmp, i_rec_clk_tmp);
            Nt = Nx_recclk * iRecClk;
            N = N(~i_rec_clk_tmp, ~i_rec_clk_tmp) - Nt * N(i_rec_clk_tmp, ~i_rec_clk_tmp);
            B_recclk = B(i_rec_clk_tmp);
            B = B(~i_rec_clk_tmp) - Nt * B_recclk;
            
            x_reduced = N\B;
            
            % ------- fix the ambiguities
            
            if sum(this.param_class == this.PAR_AMB) > 0 && false
                % get the ambiguity inverse matrxi
                idx_amb = find(this.class_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) == this.PAR_AMB);
                n_amb = length(idx_amb);
                amb_y  = sparse(idx_amb,1:length(idx_amb),ones(size(idx_amb)),size(N,1),numel(idx_amb));
                C_amb_amb = N \ amb_y;
                C_amb_amb(idx_amb,:) = [];
                amb_float = x_reduced(idx_amb);
                [amb_fixed, is_fixed, l_fixed] = Fixer.fixAmbiguities(amb_float, C_amb_amb, approach);
                if is_fixed
                    Nt = N(~idx_amb,idx_amb);
                    B(idx_amb) = [];
                    B = B - sum(Nt * spdiags(amb_fixed,0,n_amb,n_amb),2);
                    N(idx_amb,idx_amb) = [];
                end
                x_fixed = N\B;
                x_reduced(idx_amb) = amb_fixed;
                x_reduced(~idx_amb) = x_fixed;
            end
            
            % ------- substitute back
            x_est(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) = x_reduced;
            
            % receiver clcok
            B_recclk = B_recclk - sum(Nx_recclk' * spdiags(x_reduced,0,length(x_reduced),length(x_reduced)),2);
            x_rec_clk = iRecClk * B_recclk;
            x_est(idx_reduce_rec_clk) = x_rec_clk;
            
            % satellite clcok
            if sat_clk
            n_sat_clk = size(B_recclk,1);
            idx_est = ~idx_reduce_iono & ~idx_reduce_sat_clk ;
            B_satclk = B_satclk -   sum(Nx_satclk' * spdiags(x_est(idx_est),0,sum(idx_est),sum(idx_est)),2);
            x_sat_clk = iSatClk * B_satclk;
            x_est(idx_reduce_sat_clk) = x_sat_clk;
            end
            
            % iono
            if iono
            n_iono = size(B_iono,1);
            idx_est = ~idx_reduce_iono;
            B_iono = B_iono -   sum(Nx_iono' * spdiags(x_est(idx_est),0,sum(idx_est),sum(idx_est)),2);
            x_iono = iIono * B_iono;
            x_est(idx_reduce_iono) = x_iono;
            end
            x = zeros(n_par,1);
            idx_est = true(n_par,1);
            idx_est(this.idx_rd) = false;
            x(idx_est) = x_est;
            res = this.obs - A*x_est;
            this.x;
            
        end
        
        function applyWeightingStrategy()
        end
        
        function res = getResidual(this)
        end
        
        function setUpSA(this,rec_work,id_sync,flag,param_selction)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if strcmpi(flag,'???')
               this.addObsEq(rec_work, rec_work.getObsSet('L??'), param_selction);
               this.addObsEq(rec_work, rec_work.getObsSet('C??'), param_selction);
            else
                 this.addObsEq(rec_work, rec_work.getObsSet(flag), param_selction);
            end
            ls_param = LS_Parametrization();
            this.bondParamsGenerateIdx(ls_param);
        end
        
        function setUpPPP(this,rec_work,id_sync)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            param_selction = [this.PAR_REC_X this.PAR_REC_Y this.PAR_REC_Z this.PAR_REC_EB this.PAR_AMB  this.PAR_REC_CLK this.PAR_TROPO this.PAR_TROPO_N this.PAR_TROPO_E this.PAR_IONO]; %
            this.setUpSA(rec_work,id_sync,'???',param_selction);
        end
        
        function setUpIonoFreePPP(this,rec_work,id_sync)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            param_selction = [this.PAR_REC_X this.PAR_REC_Y this.PAR_REC_Z this.PAR_REC_EB this.PAR_AMB   this.PAR_REC_CLK this.PAR_TROPO this.PAR_TROPO_N this.PAR_TROPO_E]; 
           this.addObsEq(rec_work, rec_work.getPrefIonoFree('L','G'), param_selction);
            this.addObsEq(rec_work, rec_work.getPrefIonoFree('C','G'), param_selction);

            ls_param = LS_Parametrization();
            this.bondParamsGenerateIdx(ls_param);
        end
        
         function setUpNET(this,sta_list,flag)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            param_selction = [this.PAR_REC_X this.PAR_REC_Y this.PAR_REC_Z this.PAR_REC_EB this.PAR_REC_CLK this.PAR_IONO this.PAR_SAT_CLK]; %this.PAR_TROPO this.PAR_TROPO_N this.PAR_TROPO_E 
            for r = 1 : length(sta_list)
            if strcmpi(flag,'???')
                this.addObsEq(sta_list(r).work, sta_list(r).work.getObsSet('L??'),param_selction );
                this.addObsEq(sta_list(r).work, sta_list(r).work.getObsSet('C??'), param_selction);
            end
            end
            ls_param = LS_Parametrization();
            this.bondParamsGenerateIdx(ls_param);
        end
        
    end
    
end