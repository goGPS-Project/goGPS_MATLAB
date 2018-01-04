%   CLASS Least_Square_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Efficently manipulate sparse least squares system
%
% EXAMPLE
%   LSM = Least_Square_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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
classdef Least_Squares_Manipulator < handle
    properties
        A_ep % Stacked epochwise design matrices [n_obs x n_param_per_epoch]
        A_idx % index of the paramter [n_obs x n_param_per_epoch]
        out_idx % index to tell if observation is outlier [ n_obs x 1]
        N_ep  % Stacked epochwise normal matrices [ n_param_per_epoch x n_param_per_epoch x n_obs]
        y % observations  [ n_obs x 1]
        w % observation weight [ n_obs x 1]
        epoch % epoch of the obseravtions and of the A lines [ n_obs x 1]
        sat % satellite of the obseravtions and of the A lines [ n_obs x 1]
        n_epochs
        param_class % [n_param_per_epoch x 1] each paramter can be part of a class 
        %   [ 1 : x
        %     2 : y
        %     3 : z
        %     4 : inter channel/frequency/system biases,
        %     5 : ambiguity, 
        %     6 : clock
        %     7 : tropo
        %     8 : tropo inclination north
        %     9 : tropo inclination east ]
        param_flag % 0: constant in time always same param, -1: constant in time differents param (e.g ambiguity), 1: same param changing epochwise
        time_regularization %[ param_class time_varability] simple time regularization constructed from psudo obs p_ep+1 - p_ep = 0 with given accuracy
        
    end
    properties(Access = private)
        Ncc % part of the normal matrix with costant paramters
        Nee % diagonal part of the normal matrix with epoch wise or multi epoch wise paramters
        Nce % cross element between constant and epoch varying paramters
        state 
    end
    methods
        function this = Least_Squares_Manipulator()
            this.state = Go_State.getCurrentSettings();
        end
        function regularize(this, reg_opt)
        end
        function setUpPPP(this, rec, ppp_opt)
            % get double frequency iono_free for all the systems
            obs_set = Observation_Set();
            for s = rec.cc.sys_c
                obs_set.merge(rec.getPrefIonoFree('L', s));
            end
            
            snr_to_fill = (double(obs_set.snr ~= 0) + 2 * double(obs_set.obs ~= 0)) == 2; % obs if present but snr is not
            if sum(sum(snr_to_fill));
                obs_set.snr = simpleFill1D(obs_set.snr, snr_to_fill);
            end
            obs_set.remUnderCutOff(15);
            % set up number of parametrs requires
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            diff_obs = nan2zero(zero2nan(obs_set.obs)-zero2nan(synt_obs));
            % remove not valid empty epoch or with only one satellite (too
            % bad conditioned)
            idx_valid_ep_l = sum(diff_obs ~= 0, 2) > 1;
            diff_obs(~idx_valid_ep_l, :) = [];
            xs_loc(~idx_valid_ep_l, :, :) = [];
            idx_valid_stream = sum(diff_obs, 1) ~= 0;
            diff_obs(:, ~idx_valid_stream) = [];
            xs_loc(:, ~idx_valid_stream, :) = [];
            
            obs_set.remEpochs(~idx_valid_ep_l);
            obs_set.sanitizeEmpty();
            n_epochs = size(obs_set.obs, 1);
            this.n_epochs = n_epochs;
            n_stream = size(obs_set.obs, 2);
            n_coo = 3;
            n_clocks = n_epochs;
            n_tropo = n_clocks;
            ep_p_idx = [1 : n_clocks];
            
            u_obs_code = cell2mat(unique(cellstr(obs_set.obs_code)));
            iob_idx = zeros(size(obs_set.wl));
            for c = 1:size(u_obs_code, 1)
                idx_b = idxCharLines(obs_set.obs_code, u_obs_code(c, :));
                iob_idx(idx_b) = c - 1;
            end
            iob_p_idx = iob_idx + n_coo;
            cycle_slip = obs_set.cycle_slip;
            cycle_slip(diff_obs == 0) = 0;
            amb_idx = ones(size(cycle_slip));
            for s = 1:n_stream
                if s > 1
                    amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                end
                cs = find(cycle_slip(:, s) > 0)';
                for c = cs
                    %check if cycle slip is not marked at first epoch of
                    %the stream
                    if c ~= find(diff_obs(:, s) ~= 0, 1, 'first')
                        amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                    end
                end
            end
            
            
            n_iob = size(u_obs_code, 1) - 1;
            
            n_obs = sum(sum(diff_obs ~= 0));
            n_amb = max(max(amb_idx));
            clocks_idx = n_coo + n_iob + n_amb + ep_p_idx;
            
            A = zeros(n_obs, n_coo+n_iob+2); % three coordinates, 1 clock, 1 ineter obs bias(can be zero), 1 amb, 3 tropo paramters
            epoch = zeros(n_obs, 1);
            sat = zeros(n_obs, 1);
            A_idx = zeros(n_obs, n_coo+n_iob+2); % three coordinates, 1 clock, 1 ineter obs bias(can be zero), 1 amb, 3 tropo paramters
            A_idx(:, 1:3) = repmat([1, 2, 3], n_obs, 1);
            y = zeros(n_obs, 1);
            w = zeros(n_obs, 1);
            obs_count = 1;
            for s = 1:n_stream
                vaild_ep_stream = diff_obs(:, s) ~= 0;
                
                obs_stream = diff_obs(vaild_ep_stream, s);
                snr_stream = obs_set.snr(vaild_ep_stream, s);
                el_stream = obs_set.el(vaild_ep_stream, s) / 180 * pi;
                az_stream = obs_set.az(vaild_ep_stream, s) / 180 * pi;
                xs_loc_stream = permute(xs_loc(vaild_ep_stream, s, :), [1, 3, 2]);
                los_stream = rowNormalize(xs_loc_stream);
                n_obs_stream = length(obs_stream);
                lines_stream = obs_count + (0:(n_obs_stream - 1));
                epoch(lines_stream) = ep_p_idx(vaild_ep_stream);
                sat(lines_stream) = s;
                y(lines_stream) = obs_stream;
                w(lines_stream) = snr_stream;
                A(lines_stream, 1:3) = - los_stream;
                if n_iob > 0
                    A(lines_stream, 4) = iob_idx(s) > 0;
                    A_idx(lines_stream, 4) = max(n_coo+1, iob_p_idx(s));
                end
                
                A(lines_stream, n_coo+n_iob+1) = obs_set.wl(s);
                A_idx(lines_stream, n_coo+n_iob+1) = n_coo + n_iob + amb_idx(vaild_ep_stream, s);
                A(lines_stream, n_coo+n_iob+2) = 1;
                A_idx(lines_stream, n_coo+n_iob+2) = n_coo + n_iob + n_amb + ep_p_idx(vaild_ep_stream);
                sine = sin(el_stream);
                cose = cos(el_stream);
                not_inf_factor = 0.01;
                if ppp_opt.tropo
                A(lines_stream, n_coo+n_iob+3) = 1 ./ (sine + not_inf_factor);
                A_idx(lines_stream, n_coo+n_iob+3) = n_coo + n_clocks + n_iob + n_amb + ep_p_idx(vaild_ep_stream);
                end
                if ppp_opt.tropo_g
                derivative_term = -cose ./ (sine + not_inf_factor).^2;
                A(lines_stream, n_coo+n_iob+4) = cos(az_stream) .* derivative_term;
                A(lines_stream, n_coo+n_iob+5) = sin(az_stream) .* derivative_term;
                
                A_idx(lines_stream, n_coo+n_iob+4) = n_coo + 2 * n_clocks + n_iob + n_amb + ep_p_idx(vaild_ep_stream);
                A_idx(lines_stream, n_coo+n_iob+5) = n_coo + 3 * n_clocks + n_iob + n_amb + ep_p_idx(vaild_ep_stream);
                end
                obs_count = obs_count + n_obs_stream;
            end
            % ---- Suppress weighting until solution is more stable/tested
            w(:) = 1;%0.005;%this.state.std_phase;
            %---------------------
            
            this.A_ep = A;
            this.A_idx = A_idx;
            this.w = w;
            this.y = y;
            this.epoch = epoch;
            this.sat = sat;
            this.param_flag = [0, 0, 0, -1*ones(n_iob), -1, 1, 1*ones(ppp_opt.tropo), 1*ones(ppp_opt.tropo_g), 1*ones(ppp_opt.tropo_g)];
            this.param_class = [1, 2, 3, 4*ones(n_iob), 5, 6, 7*ones(ppp_opt.tropo), 8*ones(ppp_opt.tropo_g), 9*ones(ppp_opt.tropo_g)];
        end
        function setTimeRegularization(this, param_class, time_variability)
            idx_param = this.time_regularization == param_class;
            if sum(idx_param) > 0
                this.time_regularization(idx_param, 2) = time_variability;
            else %if not prestn add it
                this.time_regularization = [this.time_regularization; [param_class, time_variability]];
            end
        end
        function Astack2Nstack(this)
            %DESCRIPTION: generate N stack A'*A
            n_obs = size(this.A_ep, 1);
            this.N_ep = zeros(size(this.A_ep, 2), size(this.A_ep, 2), n_obs);
            for i = 1:n_obs
                A_l = this.A_ep(i, :);
                w = this.w(i).^2;
                this.N_ep(:, :, i) = (w * A_l)' * A_l;
            end
        end
        function res = getResiduals(this, x)
            res_l = zeros(size(this.y));
            for o = 1:size(this.A_ep, 1)
                res_l(o) = this.y(o) - this.A_ep(o, :) * x(this.A_idx(o, :), 1);
            end
            n_epochs = max(this.epoch);
            n_sat = max(this.sat);
            res = zeros(n_epochs, n_sat);
            for i = 1:n_sat
                idx = this.sat == i;
                ep = this.epoch(idx);
                res(ep, i) = res_l(idx);
            end
        end
        function [x, res, s02, Cxx] = solve(this)
            idx_constant_l = this.param_flag == 0 | this.param_flag == -1;
            idx_constant = find(idx_constant_l);
            idx_non_constant = find(~idx_constant_l);
            n_constant = max(max(this.A_idx(:, idx_constant_l)));
            n_class = size(this.A_ep, 2);
            n_ep_wise = max(max(this.A_idx(:, ~idx_constant_l))) - n_constant;
            n_epochs = this.n_epochs;
            n_obs = size(this.A_ep, 1);
            n_ep_class = n_ep_wise / n_epochs;
            Ncc = zeros(n_constant, n_constant);
            Nce = zeros(n_ep_wise, n_constant);
            n_class_ep_wise = length(idx_non_constant);
            Ndiags = zeros(n_class_ep_wise, n_class_ep_wise, n_epochs); %permute(this.N_ep(~idx_constant_l,~idx_constant_l,:),[3,1,2]);
            B = zeros(n_constant+n_ep_wise, 1);
            for i = 1:n_obs
                p_idx = this.A_idx(i, :);
                p_idx(p_idx == 0) = 1;  % does not matter since terms are zeros
                N_ep = this.N_ep(:, :, i);
                A_ep = this.A_ep(i, :);
                w = this.w(i);
                y = this.y(i);
                e = this.epoch(i);
                p_c_idx = p_idx(idx_constant_l);
                p_e_idx = p_idx(~idx_constant_l) - n_constant;
                p_e_idx(p_e_idx <= 0) = 1;  % does not matter since terms are zeros
                
                
                % fill Ncc
                Ncc(p_c_idx, p_c_idx) = Ncc(p_c_idx, p_c_idx) + N_ep(idx_constant, idx_constant);
                % fill Nce
                Nce(p_e_idx, p_c_idx) = Nce(p_e_idx, p_c_idx) + N_ep(idx_non_constant, idx_constant);
                %fill Ndiags
                
                Ndiags(:, :, e) = Ndiags(:, :, e) + N_ep(idx_non_constant, idx_non_constant);
                %fill B
                B(p_idx) = B(p_idx) + A_ep' * w * y;
            end
            Nee = [];
            class_ep_wise = this.param_class(idx_non_constant);
            reg_diag0 = ones(n_epochs, 1) + [0; ones(n_epochs-2, 1); 0];
            reg_diag1 = -ones(n_epochs-1, 1);
            Ndiags = permute(Ndiags, [3, 1, 2]);
            for i = 1:n_ep_class
                N_col = [];
                for j = 1:n_ep_class
                    diag0 = Ndiags(:, i, j);
                    N_el = sparse(n_epochs, n_epochs);
                    if j == i
                        cur_class = class_ep_wise(i);
                        idx_c = this.time_regularization(:, 1) == cur_class;
                        w = (1 ./ this.time_regularization(idx_c, 2)).^2;
                        if sum(idx_c)
                            diag0 = diag0 + reg_diag0 * w;
                            diag1 = reg_diag1 * w;
                            N_el = spdiags([0; diag1], 1, N_el);
                            N_el = spdiags(diag1, -1, N_el);
                        end
                    end
                    N_el = spdiags(diag0, 0, N_el);
                    N_col = [N_col; N_el];
                end
                Nee = [Nee, N_col];
            end
            N = [[Ncc, Nce']; [Nce, Nee]];
            if nargout > 3
                %inverse by partitioning, taken from:
                % Mikhail, Edward M., and Friedrich E. Ackermann. "Observations and least squares." (1976). pp 447
                Ncc = sparse(Ncc);
                Nce = sparse(Nce);
                invNcc = (Ncc)^(-1);
                invNee = (Nee)^(-1);
                a22ia21 = invNee * Nce;
                invN11 = (Ncc - (Nce') * a22ia21)^(-1);
                invN12 = -invN11 * (a22ia21');
                invN21 = invN12';
                invN22 = invNee - a22ia21 * invN12;
                Cxx = [[invN11; invN21], [invN12; invN22]];
                x = Cxx * B;
            else
                x = N \ B;
            end
            x_class = zeros(size(x));
            for c = 1:length(this.param_class)
                x_class(this.A_idx(:, c)) = this.param_class(c);
            end
            if nargout > 1
                res = this.getResiduals(x);
            end
            x = [x, x_class];
            
        end
    end
end