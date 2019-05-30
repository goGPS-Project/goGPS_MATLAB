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
        
        unique_obs_codes % uniques ids (cell) of the signals  /since lot of combinations are possible they will be generated dynamically)
        unique_wl % set of unique wavelength
        rec_xyz % receiver coordinates to be used in
        unique_sat_goid % unique satellite goids
        cycle_slips % epoch of the cycle slip
        
        time_par   % time of the paramter
        rec_par    % receiver of the paramters
        sat_par    % staellite of the paramters
        class_par  % class of the paramter
        obs_codes_id_par  % obs code id fo the paramter
        
        rec_set % set of receivers
        sat_set % set of satellites
        ch_set  % set of observation codes
        
        N      
    end
    
    methods
        function addObsEq(this,rec,obs_set, param_selection)
            % add observation equations to the matrices
            %
            % SYNTAX:
            %    this.addObsEq(rec, obs_set)
            n_param = length(param_selection);
            
            % --- check all paramters presence and their order -----
            par_rec_x_lid = param_selection == this.PAR_REC_X;
            par_rec_x = sum(par_rec_x_lid) > 0;
            par_rec_y_lid = param_selection == this.PAR_REC_Y;
            par_rec_y = sum(par_rec_y_lid) > 0;
            par_rec_z_lid = param_selection == this.PAR_REC_Z;
            par_rec_z = sum(par_rec_z_lid) > 0;
            
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
                if Core_Utils.findAinB(u_obs_code(i),this.unique_obs_codes) > 0
                    this.unique_obs_codes{end+1} = u_obs_code{i};
                end
            end
            obs_code_id = Core_Utils.findAinB(obs_set.obs_code,this.unique_obs_codes);
            
            
            % add wavelength to the unique
            u_wl = unique(obs_set.wl);
            this.unique_wl = unique([this.unique_wl u_wl]);
            [~, obs_wl_id] = intersect(obs_set.wl, this.unique_wl);
            
            
            % get the mapping function for tropo
            if sum(param_selection == this.PAR_TROPO || param_selection == this.PAR_TROPO_E || param_selection == this.PAR_TROPO_N || param_selection == this.PAR_TROPO_V) > 0
                id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
            % check whivh observations are phase ones
            phase_s = obs_set.wl ~=  -1;
            
            % initliaze the matrices
            A = zeros(n_obs, n_par);
            [obs,satellite_obs, azimuth_obs, elevation_obs, variance_obs, wl_obs] = deals(zeros(n_obs, 1));
            [ obs_codes_id_obs,  wl_obs] = deals(zeros(n_obs, 1,'uint8'));
            phase_obs = false(n_obs,1);
            time_obs = GPS_Time();
            obs_count = 1;
            
            
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
                    time_obs.addEpoch(obs_set.time.getEpoch(id_ok_stream));
                    
                    % ----------- save the cycle slips
                    if phase_s(s)
                        this.cycle_splips{r}{s_go_id}{s_s_id} = obs_set.time.getEpoch(find(obs_set.cycle_slip(:,s)));
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
                        A(lines_stream, par_rec_sat_lid) = 1;
                    end
                    % ----------- ZTD ------------------
                    if par_tropo
                        A(lines_stream, par_rec_tropo_lid) = mfw_stream;
                    end
                    % ----------- ZTD gradients ------------------
                    if par_tropo_n || par_tropo_e
                        cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                        if par_tropo_e
                            A(lines_stream, par_rec_tropo_e_lid) = sin(az_stream) .* cotan_term; % east gradient
                        end
                        if par_tropo_n
                            A(lines_stream, par_rec_tropo_n_lid) = cos(az_stream) .* cotan_term; % noth gradient
                        end
                    end
                    if par_tropo_v
                        A(lines_stream, par_rec_tropo_v_lid) = mfw_stream*rec.h_ellips;
                    end
                    % ----------- Ionosphere delay --------------------
                    if par_iono
                        if phase_s(s)
                            A(lines_stream, par_iono_lid) = - obs_set.wl^2;
                        else
                            A(lines_stream, par_iono_lid) =   obs_set.wl^2;
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
            this.tme_obs.addEpoch(time_obs);
            this.receiver_obs = [this.receiver_obs; r*ones(size(phase_obs))];
        end
        
        function bondParamsGenerateIdx(this, ls_parametrization)
            % bond paramters (splines or other models) and generate idx
            %
            % SYNTAX
            %    this.bondParamGenerateIdx(parametrization)
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            this.A_idx = zeros(size(this.A));
            time_min = min(this.time_obs.getMatlabTime);
            obs_rate = this.time_obs.getRate; % TBD substitute this quantity with the obesravtion minimum rate
            time_obs = round(this.time_obs.getRefTime(time_min)/obs_rate);
            
            cumulative_idx = 0;
            i_col = 1;
            for i_p = 1 : length(this.param_class)
                p = this.param_class(i_p);
                [parametriz, opt] = ls_parametrization.getParametrization(p);
                  % ----------------- defining receiver sets ---------
                  if parametriz(2) == ls_parametrization.SING_REC
                      n_rec_set = size(this.rec_xyz,1);
                      rec_set = mat2cell(1:n_rec_set);
                  elseif parametriz(2) == ls_parametrization.ALL_REC
                      n_rec_set = 1;
                      rec_set = {1:n_rec_set};
                  elseif parametriz(2) == ls_parametrization.MULTI_REC
                      n_rec_set = length(opt.rec_sets);
                      rec_set = opt.rec_sets;
                  end
                  % ----------------- defining satellites sets ---------
                  if parametriz(2) == ls_parametrization.SING_SAT
                      n_sat_set = length(this.unique_sat_goid);
                      sat_set = mat2cell(this.unique_sat_goid);
                  elseif parametriz(2) == ls_parametrization.ALL_SAT
                      n_sat_set = 1;
                      sat_set = {this.unique_sat_goid};
                  elseif parametriz(2) == ls_parametrization.MULTI_SAT
                      n_sat_set = length(opt.sat_sets);
                      sat_set = opt.sat_sets;
                  end
                  end
                  % ----------------- defining channel sets ---------
                  % construct channel to frequency mapping
                  sig2wl = zeros(size(this.unique_obs_codes));
                  sig_p_id = 1:size(this.unique_obs_codes,1);
                  for c = 1 : length(this.unique_obs_codes)
                      idx1 = find(this.obs_codes_id_obs == c,1,'first');
                      sig2wl(c) = this.wl_obs(idx1);
                  end
                  % construct channel to phase mapping
                  sig2phase = zeros(size(this.unique_obs_codes));
                  for c = 1 : length(this.unique_obs_codes)
                      idx1 = find(this.obs_codes_id_obs == c,1,'first');
                      sig2phase(c) = this.phase_obs(idx1);
                  end
                  % construct channel to constellation mapping
                  sig2const = char(zeros(size(this.unique_obs_codes)));
                  for c = 1 : length(this.unique_obs_codes)
                      idx1 = find(this.obs_codes_id_obs == c,1,'first');
                      ant_id = Core.getConstellationCollector.getAntennaId(this.go_id(idx1));
                      sig2const(c) = ant_id(1);
                  end
                  if parametriz(3) == ls_parametrization.SING_TRACK
                      n_ch_set = size(this.unique_obs_codes);
                      ch_set = mat2cell(1:n_ch_set);
                  elseif parametriz(3) == ls_parametrization.ALL_FREQ
                      n_ch_set = 1;
                      ch_set = {1:n_ch_set};
                  elseif parametriz(3) == ls_parametrization.SING_FREQ
                      n_ch_set = length(this.unique_wl);
                      ch_set = {};
                      for c = 1 : n_ch_set
                          ch_set{c} = sig_p_id(sig2wl == this.unique_wl(c));
                      end
                  elseif parametriz(3) == ls_parametrization.RULE
                      % ---- evaluate the rule
                      n_ch_set = 0;
                      ch_set = {};
                      for rule = opt.eb_rule
                          cond_sig = true(size(sig_p_id));
                          parts = strsplit(rule,':');
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
                              elseif strcmpi('NOT_GLONASS',c)
                                  cond_sig = cond_sig & sig2const ~= 'R';
                              else
                                  cond_sig(:) = false;
                                  this.logWarning(sprintf('Unknown option %s in parametrization rule, ignoring parameter',c));
                              end
                          end
                          if any(cond_sig)
                              if str2num(opt) == ls_parametrization.ALL_FREQ
                                  n_ch_set = n_ch_set +1;
                                  ch_set{n_ch_set} = sig_p_id(cond_sig);
                              elseif str2num(opt) == ls_parametrization.SING_FREQ
                                  u_wl_tmp = unique(sig2wl(cond_sig));
                                  n_wl_tmp = length(u_wl_tmp);
                                  for c = 1 : n_wl_tmp
                                      n_ch_set = n_ch_set +1;
                                      ch_set{n_ch_set} = sig_p_id(sig2wl == u_wl_tmp(c) & cond_sig);
                                  end
                              elseif str2num(opt) == ls_parametrization.SING_TRACK
                                  n_ch_set = n_ch_set +length(sig_p_id(cond_sig));
                                  ch_set = [ch_set mat2cell(sig_p_id(cond_sig))];
                                  
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
                          sat_lid = false(size(this.A,1),1);
                          for ss = sat_set{s}
                              sat_lid = sat_lid | this.satellite_obs == ss;
                          end
                          if length(sat_set{r}) == 1
                              r_id = sat_set{s};
                          else
                              s_id =Core_Utils.findAinB(sat_set{s},this.sat_set);
                              if s_id == 0
                                  s_id = length(this.sat_set) -1;
                                  this.sat_set{end+1} = sat_set{s};
                              end
                              
                          end
                          for f = 1 : n_chanel_set
                              ch_lid = false(size(this.A,1),1);
                              for cc = ch_set{f}
                                  ch_lid = ch_lid | this.obs_codes_id_obs == cc;
                              end
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
                              % --- now dealing with epoch dependence ----
                              cols_tmp = 0;
                              if parametriz(1) = ls_parametrization.CONST
                                  ep_pgr_id = 1;
                                  n_prg_id = 1;
                                  time_par_tmp = this.time_obs.getEpoch(obs_lid).minimum.getMean(this.time_obs.getEpoch(obs_lid).maximum);
                              elseif parametriz(1) = ls_parametrization.EP_WISE
                                  ep_id = time_obs(obs_lid);
                                  u_e_tmp = unique(ep_id);
                                  [~,ep_pgr_id] = intersect(ep_id,u_e_tmp);
                                  n_prg_id = length(u_e_tmp);
                                  time_par_tmp = this.time_obs.getEpoch(obs_lid);
                              elseif parametriz(1) = ls_parametrization.STEP_CONST
                                  ep_id = time_obs(obs_lid);
                                  ep_pgr_id = zeors(size(ep_id));
                                  n_prg_id = 0;
                                  steps_set = opt.steps_set;
                                  if p == this.PAR_AMB %% ambiguity case is really unique and make sense to treat it separatly
                                      steps = round(this.cycle_slips{rr}{ss}{cc}.getRefTime(time_min)/obs_rate);
                                      p_s = 1;
                                      for st = steps'
                                          ep_pgr_id(time_obs >= st) = p_s;
                                          p_s = p_s +1;
                                      end
                                  elseif parametriz(2) == ls_parametrization.ALL_REC % you can use differents step for step wise satellite dependent paraters
                                      steps = round(steps_set{ss}.getRefTime(time_min)/obs_rate);
                                      p_s = 1;
                                      for st = steps'
                                          ep_pgr_id(time_obs >= st) = p_s;
                                          p_s = p_s +1;
                                      end
                                  elseif parametriz(3) == ls_parametrization.ALL_SAT  % you can use differents step for step wise receiver dependent paraters
                                      steps = round(steps_set{rr}.getRefTime(time_min)/obs_rate);
                                      p_s = 1;
                                      for st = steps'
                                          ep_pgr_id(time_obs >= st) = p_s;
                                          p_s = p_s +1;
                                      end
                                  end
                              elseif parametriz(1) = ls_parametrization.SPLINE_ZERO
                                  ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                  u_e_tmp = unique(ep_id);
                                  [~,ep_pgr_id] = intersect(ep_id,u_e_tmp);
                                  n_prg_id = length(u_e_tmp);
                              elseif parametriz(1) = ls_parametrization.SPLINE_LIN
                                  % ---- colum will be doubled -----
                                  cols_tmp = [ 0 1];
                                  ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                  spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate),1);
                                  u_e_tmp = unique([ep_id ep_id+1]);
                                  ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                  for i_o = cols_tmp;
                                      [~,ep_pgr_id(i_o+1)] = intersect(ep_id+i_o,u_e_tmp);
                                  end
                                  % ----- expand colum of the A matrix
                                  this.A = [this.A(:,1:(i_col-1)) this.A(:,i_col).*spline_v(:,1) this.A(:,i_col).*spline_v(:,2) this.A(:,(i_col+1):end)];
                                  n_prg_id = length(u_e_tmp);
                              elseif parametriz(1) = ls_parametrization.SPLINE_CUB
                                   % ---- colum will be quadrupled -----
                                  cols_tmp = [ 0 1 2 3];
                                  ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                  spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate),1);
                                  u_e_tmp = unique([ep_id ep_id+1 ep_id+2 ep_id+3]);
                                  ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                  for i_o = cols_tmp;
                                      [~,ep_pgr_id(i_o+1)] = intersect(ep_id+i_o,u_e_tmp);
                                  end
                                  % ----- expand colum of the A matrix
                                  this.A = [this.A(:,1:(i_col-1)) this.A(:,i_col).*spline_v(:,1) this.A(:,i_col).*spline_v(:,2) this.A(:,i_col).*spline_v(:,3) this.A(:,i_col).*spline_v(:,4) this.A(:,(i_col+1):end)];
                                  n_prg_id = length(u_e_tmp);
                              end
                              this.A_idx(obs_lid,i_col + cols_tmp) = cumulative_idx + ep_pgr_id;
                              [u_new_par] = unique(cumulative_idx + ep_pgr_id(:));
                              
                              this.time_par.addEpoch(time_par_tmp);  % time of the paramter
                              this.rec_par    = [this.rec_par;                r_id*ones(n_prg_id,1)];  % receiver of the paramters
                              this.sat_par    = [this.sat_par;                s_id*ones(n_prg_id,1)];  % receiver of the paramters
                              this.class_par  = [this.class_par;              p*ones(n_prg_id,1)];  % class of the paramter
                              this.obs_codes_id_par = [this.obs_codes_id_par; ch_id*ones(n_prg_id,1)];  % obs_code id paramters
                              cumulative_idx = cumulative_idx + n_prg_id;
                              i_col = i_col + length(cols_tmp);
                          end
                      end
                  end
                  
            end
            
            
        end
        
        function removeFullRankDeficency(this)
            % solve full rank deficency removing paramters from the
            % estimation
            %
            % SYNTAX:
            %    this.removeFullRankDeficency()
        end
        
        function absValRegularization(this,param_id, var)
            % regularize paramters to zero (Tykhnov aka ridge aka L2)
            % 
            % this.absValRegularization(param_id, var)
        end
        
        function timeRegularization(this, param_id, var)
            % first order tykhonv regualrization in time
            %
            % this.timeRegularization(this, param_id, var)
        end
        
        function spatialRegularization(this, law)
             % tykhonv regualrization in space
            %
            % this.spatialRegularization(this, law)
        end
        
        function hemisphereRegularization(this, law)
            % tykhonv regualrization on the hemisphere
            %
            % this.spatialRegularization(this, law)
        end
        
        function reduceForNuisanceParameters(this, param_id)
            
        end
        
        function applyWeightingStrategy()
        end
        
        function res = getResidual(this)
        end
        
    end
    
end