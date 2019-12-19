%   CLASS Core Block
% =========================================================================
%
% DESCRIPTION
%   Class to manage pre processing utilities
%
% EXAMPLE
%   go_block = Core_Block();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_PP
%
% Note for the future: the class uses the current obs storage of goGPS
% -> switch to objects for rover and master observations is suggested

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

classdef Core_PP < handle
    
    properties (Constant, Access = private)
    end
    
    properties (Access = public)% Public Access
        log
        state
    end
    
    % ==================================================================================================================================================
    %  CREATOR
    % ==================================================================================================================================================
    
    methods (Static)
        function this = Core_PP()
            this.init()
        end                
    end
    
    % ==================================================================================================================================================
    %  PRE PROCESSING MAIN FUNCTIONS
    % ==================================================================================================================================================
    
    methods % Public Access
        function init(this)
            % Init singletons
            %
            % SYNTAX
            %   this.init()
            this.log = Core.getLogger();
            this.state = Core.getState();
        end
                
        function [ph] = jmpFix(this, ph, lambda)
            % after a missing data try to fix the jmp in the observations with an N * lambda correction
            cs_factor = this.state.cs_thr;
            ph = ph';
            for a = 1 : size(ph,2)
                arc = zero2nan(ph(:,a));
                x_arc = (1 : numel(arc))';
                arc_fit = zeros(size(arc));
                
                int_factor = zeros(size(arc));
                int_factor(~isnan(arc)) = lambda(a,1) * cs_factor;
                arc = arc ./ int_factor;
                if any(arc)
                    %PLOT:figure;
                    %PLOT:plot(x_arc, arc, 'LineWidth', 1); hold on;
                    lim = getOutliers(~isnan(arc));
                    for l = 1 : size(lim, 1) - 1
                        id = lim(1, 1) : lim(l, 2);
                        id = id(~isnan(arc(id)));
                        id = id((max(1, numel(id) - this.state.getMinArc+1) : end));
                        
                        p = polyfit(x_arc(id), arc(id), 3);
                        if l < size(lim, 1)
                            arc_fit1 = polyval(p, x_arc([lim(l, 2) lim(l + 1, 1)]));
                            
                            id2 = lim(l + 1, 1) : lim(l + 1, 2);
                            id2 = id2(1 : min(numel(id), this.state.getMinArc));
                            p2 = polyfit(x_arc(id2), arc(id2), 3);
                            arc_fit2 = polyval(p2, x_arc([lim(l, 2) lim(l + 1, 1)]));
                            
                            arc(lim(l + 1, 1) : end) = arc(lim(l + 1, 1) : end) + round(mean(arc_fit1 - arc_fit2));
                        end
                        %PLOT: plot(x_arc, arc);
                    end
                    %PLOT: plot(x_arc, arc, 'LineWidth', 2); hold on;
                end
                ph(:,a) = nan2zero(arc) .* int_factor;
            end
            ph = ph';
            
            %PLOT: figure; plot(zero2nan(ph));
            %PLOT: d_ph = diff(zero2nan([ph(1,:); ph])); figure; plot(d_ph - movmedian(d_ph, 3, 'omitnan'));
        end
    end
        
    % ==================================================================================================================================================
    %  STATIC LAUNCHERS
    % ==================================================================================================================================================
    
    methods (Static) % Public Access
        function [data, flagged_out] = flagRawObs(data, time_set, time_ref, thr_factor, thr_min, win_size)
            % Flag the data above max(thr_min, thr_factor) * mean(movstd))
            % SYNTAX: data = flagRawObs(data,  time_set, time_ref, thr_factor, <thr_min = 0 >)
            if nargin < 5
                thr_min = 0;
            end
            if nargin < 6
                win_size = 5;
            end
                        
            flag = 1;
            i = 0;
            % Interpolate to reference time to find outliers in the 4th derivative
            if ~isempty(time_set)
                data_interp = zeros(size(data));
                for s = 1 : size(data,2)
                    tmp = nan2zero(data(:,s));
                    if any(tmp)
                        id_interp = (conv(tmp ~= 0, ones(9,1), 'same') > 0);
                        data_interp(id_interp,s) = interp1(time_set, zero2nan(data(:,s)), time_ref(id_interp), 'spline', 'extrap');
                    end
                end
            else
                data_interp = data;
            end

            % loop for outliers rejection
            flagged_out = false(size(data));
            while (sum(flag(:)) > 0) && (i < 5)
                sensor = zero2nan(data_interp);
                sensor = bsxfun(@minus, sensor, median(sensor, 2, 'omitnan'));
                i = i + 1;
                % expand the set of observations to allow flogging at the border of the valid intervals
                sensor = movstd(sensor, win_size, 'omitnan');
                flag = sensor > max(thr_min, thr_factor * mean(serialize(zero2nan(sensor .* ~isnan(zero2nan(data)))), 'omitnan'));
                data(flag) = NaN;
                data_interp(flag) = NaN;
                flagged_out  = flagged_out | flag;
            end
        end
        
        function [data, flagged_out] = flagRawObsD4(data, time_set, time_ref, thr_factor, thr_min, win_size)
            % Flag the data above max(thr_min, thr_factor) * mean(movstd))
            %
            % SYNTAX
            %   data = flagRawObsD4(data,  time_set, time_ref, thr_factor, <thr_min = 0 >)
            if nargin < 5
                thr_min = 0;
            end
            if nargin < 6
                win_size = 5;
            end
                        
            flag = 1;
            i = 0;
            % Interpolate to reference time to find outliers in the 4th derivative
            if ~isempty(time_set)
                data_interp = zeros(size(data));
                for s = 1 : size(data,2)
                    tmp = zero2nan(data(:,s));
                    id_interp = ~isnan(tmp);
                    if sum(~isnan(tmp) > 2)
                        data_interp(id_interp,s) = interp1(time_set(id_interp), tmp(id_interp), time_ref(id_interp), 'spline', 'extrap');
                    end
                end
            else
                data_interp = data;
            end

            % loop for outliers rejection
            flagged_out = false(size(data));
            while (sum(flag(:)) > 0) && (i < 5)
                sensor = Core_Utils.diffAndPred(zero2nan(data_interp), 4);
                sensor = bsxfun(@minus, sensor, median(sensor, 2, 'omitnan'));
                i = i + 1;
                % expand the set of observations to allow flogging at the border of the valid intervals
                sensor = movstd(sensor, win_size, 'omitnan');
                flag = sensor > max(thr_min, thr_factor * mean(serialize(zero2nan(sensor .* ~isnan(zero2nan(data)))), 'omitnan'));
                data(flag) = NaN;
                data_interp(flag) = NaN;
                flagged_out  = flagged_out | flag;
            end
        end
        
        function [data, flagged_out] = flagBorders(data, win_size)
            % Flag the data above max(thr_min, thr_factor) * mean(movstd))
            % SYNTAX: [data, flagged_out] = flagBorders(data, win_size)
            if nargin < 2
                win_size = 20;
            end
            
            % Interpolate to reference time
            time_set = 1 : size(data,1);
            id_interp = flagExpand(~isnan(data), 5) & isnan(data);
            data_interp = data;
            for s = 1 : size(data,2)
                tmp = data(~isnan(data(:,s)),s);
                if numel(tmp) > 5
                    data_interp(id_interp(:,s),s) = interp1(time_set(~isnan(data(:,s))), tmp, time_set(id_interp(:,s)), 'spline', 'extrap');
                end
            end
            
            % loop for outliers rejection
            sensor = Core_Utils.diffAndPred(zero2nan(data_interp), 4);
            sensor = bsxfun(@minus, sensor, median(sensor, 2, 'omitnan'));
            % expand the set of observations to allow flagging at the border of the valid intervals
            sensor(isnan(data)) = 0;
            sensor = maxfilt_mat(movstd(sensor, ceil(win_size/2), 'omitnan'),5);
            sensor(isnan(data)) = 0;
            std_max = sensor; std_max(flagExpand(sensor == 0, win_size)) = 0;
            sensor_max = maxfilt_mat(max(std_max, [], 2), ceil(win_size/2));
            sensor_max = simpleFill1D(sensor_max, sensor_max == 0);
            % remove possible outliers over thr and small arcs < 5 epochs
            flagged_out = ~Core_PP.remShortArcs(~(sensor > sensor_max | isnan(data))', 5)' & ~isnan(data);
            data(flagged_out) = NaN;
        end
    end
    
    % ==================================================================================================================================================
    %  STATIC public misc utilities
    % ==================================================================================================================================================
    
    methods (Static) % Public Access
        function [obs, dt_dj, is_jumping] = remDtJumps(obs)
            % remove a jumps of the clock from ph/pr
            % this piece of code is very very criptic, but it seems to work
            % review this whenever possible
            %
            % SYNTAX:
            %   [obs, dt_dj, is_jumping] = remDtJumps(obs)
            obs = zero2nan(obs);
            dt_dj = 0;
            is_jumping = false;
            
            % Iterate for very particular instable cases
            for i = 1 : 2
                % Try to find the biggest jumps (experimental)
                ddt = median(Core_Utils.diffAndPred(zero2nan(obs),1), 2, 'omitnan');
                pos_jmp = abs(ddt) > 1e5;
                d3dt = median(Core_Utils.diffAndPred(zero2nan(obs),3), 2, 'omitnan');
                ddt = cumsum(cumsum(nan2zero(d3dt)));
                dt0 = 0;
                if sum(pos_jmp) > 0
                    is_jumping = true;
                    ddt = ddt - simpleFill1D(ddt, pos_jmp);
                    dt0 = cumsum(nan2zero(ddt));
                    ph_bk = obs;
                    obs = bsxfun(@minus, obs, dt0);
                    
                    % approach on 3rd and 4th derivate (classic)
                    ddt_sensor = cumsum(cumsum(nan2zero(d3dt)) - medfilt_mat(cumsum(nan2zero(d3dt)), 3));
                    % check if there is any discontinuity in the clock drift
                    clock_thresh = 1e3;
                    pos_jmp = abs(ddt_sensor - medfilt_mat(ddt_sensor, 3)) > clock_thresh;
                    pos_jmp(1:2) = 0;
                    if sum(pos_jmp) > 0
                        ddt = ddt - simpleFill1D(ddt, pos_jmp);
                        dt = cumsum(ddt);
                        ph_bk = obs;
                        obs = bsxfun(@minus, obs, dt);
                        d4dt = median(Core_Utils.diffAndPred(zero2nan(obs),4), 2, 'omitnan');
                        dobs = Core_Utils.diffAndPred(obs);
                        % find jmps on the median 4th derivate
                        jmp_candidate = flagExpand(abs(nan2zero(d4dt)) > clock_thresh, 2);
                        % detect the real jmp index
                        lim = getOutliers(jmp_candidate);
                        jmp = false(size(jmp_candidate));
                        for l = 1 : size(lim,1)
                            dtmp = dobs(lim(l,1) : lim(l,2), :);
                            [~, id_max] = max(abs(nan2zero(median(dtmp, 2 ,'omitnan')))');
                            jmp(id_max + lim(l,1) - 1) = true;
                                                        
                            % find the magnitude of the jump
                            dtmp_f = simpleFill1D(dtmp, isnan(dtmp)); % fill temp data
                            tmp_diff = diff(dtmp_f); 
                            tmp_diff(:, isnan(dtmp(id_max,:))) = nan; % do not use filled data for estimating the magnitude of jump
                            
                            d4dt(id_max + lim(l,1) - 1) = median(sign(dobs(id_max + lim(l,1) - 1,:)),'omitnan') * ...
                                median(abs(max(tmp_diff)),'omitnan');
                        end
                        % velocity offsets beteween sat
                        d4dt(~jmp) = 0;
                        obs = bsxfun(@minus, obs, cumsum(nan2zero(d4dt)));
                        dt = dt + cumsum(d4dt);
                    else
                        dt = zeros(size(ddt));
                    end
                    dt_dj = dt_dj + dt0 + dt;
                end
            end
            dt_dj = dt_dj / 299792458;
        end
        
        function [ph_out, dt] = remDtHF(ph)
            % [ !! ] EXPERIMENTAL FUNCTION
            % Remove high frequency from phases estimated by the 4th derivative
            % WARNING: this code works but it has a lot of border effects, especially if the dataset it's not continuous
            %
            % SYNTAX:
            %   [ph_out, dt] = remDtHF(ph)
            d4dt = median(Core_Utils.diffAndPred(zero2nan(ph),4), 2, 'omitnan');
                        
            % Filter low frequencies:
%             dt_hf = cumsum(cumsum(cumsum(cumsum(nan2zero(d4dt)))));
%             x = (1 : numel(dt_hf))';
%             ws = 5;
%             margin = 2 * round(ws/2);
%             xi = (1 - margin : numel(dt_hf) + margin)';
%             dt_lf = interp1(x(~isnan(dt_hf)), dt_hf(~isnan(dt_hf)), xi, 'spline', 'extrap');
%             dt_lf = splinerMat([],medfilt_mat(dt_lf, 5), 5);
%             dt_hf = dt_hf - dt_lf(margin + 1 : end - margin);
            dt = cumsum(Core_PP.highPass(nan2zero(d4dt), 313, 0.01, 101));
            dt = cumsum(Core_PP.highPass(nan2zero(dt), 177, 0.01, 101));
            dt = cumsum(Core_PP.highPass(nan2zero(dt), 143, 0.01, 101));
            dt = cumsum(Core_PP.highPass(nan2zero(dt), 107, 0.01, 101));
            ph_out = nan2zero(bsxfun(@minus, zero2nan(ph), dt));
            dt = dt / 299792458;
        end

        function [pr, ph, dt_pr, dt_ph] = correctTimeDesyncHF(time_ref, time, pr, ph)
            %   Correction of jumps in code and phase due to dtR and time de-sync
            % SYNTAX:
            %   [pr, ph, dt_pr, dt_ph] = correctTimeDesync(time_ref, time, pr, ph)

            log = Logger.getInstance();
            
            time_desync  = round((time_ref - time) * 1e7) / 1e7;
            %figure(1); clf; plot(diff(zero2nan(ph))); hold on;
            %figure(2); clf; plot(diff(zero2nan(pr))); hold on;
            
            ph_ds = bsxfun(@minus, ph, time_desync .* 299792458);
            pr_ds = bsxfun(@minus, pr, time_desync .* 299792458);
            % if adding the desync will improve the std it means that the receiver does not compensate for it
            [ph, flag] = Core_PP.testDesyncCorrection(ph, ph_ds);
            if flag
                time_desync_ph = time_desync;
                log.addMessage('Correcting phase for time desync', 100);
            else
                time_desync_ph = 0;
            end
            
            [pr, flag] = Core_PP.testDesyncCorrection(pr, pr_ds);
            if flag
                time_desync_pr = time_desync;
                log.addMessage('Correcting pseudo-ranges for time desync', 100);
            else
                time_desync_pr = 0;
            end
            clear pr_ds ph_ds;
            [ph, dt_ph_jumps] = Core_PP.remDtJumps(ph);
            [ph, dt] = Core_PP.remDtHF(ph);
            dt_ph = dt_ph_jumps + dt;
            %figure(1); plot(diff(zero2nan(ph)),'.k');

            % correct the pseudo-ranges for HF
            pr_lf = bsxfun(@minus, pr, dt .* 299792458);
            [pr, flag] = Core_PP.testDesyncCorrection(pr, pr_lf);
            if flag
                dt_pr = dt_ph;
                log.addMessage('Correcting pseudo-ranges for dt HF as estimated from phases observations', 100);
            else
                dt_pr = zeros(size(dt_ph));
            end

            % correct the pseudo-ranges for jumps
            pr_dj = bsxfun(@minus, pr, dt_ph_jumps .* 299792458);
            [pr, flag] = Core_PP.testDesyncCorrection(pr, pr_dj);
            if flag
                dt_pr = dt_pr + dt_ph_jumps;
                log.addMessage('Correcting pseudo-ranges for dt as estimated from phases observations', 100);
            end
            
            % in some receivers the pseudo-range is not continuous while the phase are ok
            [pr_ds, dt_pr_jumps] = Core_PP.remDtJumps(pr);
            [pr, flag] = Core_PP.testDesyncCorrection(pr, pr_ds);
            if flag
                dt_pr = dt_pr_jumps + dt_pr;
                log.addMessage('Correcting pseudo-ranges for dt as estimated from their observations', 100);
            end
            
            dt_ph = dt_ph + time_desync_ph;
            dt_pr = dt_pr + time_desync_pr;
            %figure(2); plot(diff(zero2nan(pr)),'.k');
        end
        
        function [pr, ph, dt_pr, dt_ph] = correctTimeDesync(time_ref, time, pr, ph)
            %   Correction of jumps in code and phase due to dtR and time de-sync
            % SYNTAX:
            %   [pr, ph, dt_pr, dt_ph] = correctTimeDesync(time_ref, time, pr, ph)

            log = Core.getLogger();
                        
            id_not_empty = time ~= 0;

            time_desync = nan(size(time));
            time_desync(id_not_empty)  = round((time(id_not_empty) - time_ref(id_not_empty)) * 1e7) / 1e7; % the rinex time has a maximum of 7 significant decimal digits
            time_desync = simpleFill1D(time_desync,isnan(time_desync),'nearest');
            
            if any(time_desync)
                [ph_dj, dt_ph_dj] = Core_PP.remDtJumps(ph);
                [pr_dj, dt_pr_dj] = Core_PP.remDtJumps(pr);
                ddt_pr = Core_Utils.diffAndPred(dt_pr_dj);
                
                %% time_desync is a introduced by the receiver to maintain the drift of the clock into a certain range
                ddt = [0; diff(time_desync)];
                ddrifting = ddt - ddt_pr;
                drifting = cumsum(ddt - ddt_pr);
                
                % Linear interpolation of ddrifting
                jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                jmp_fit = setdiff(find(abs(ddrifting) > 1e-7), jmp_reset); % points where desync interpolate the clock
                d_points = [drifting(jmp_reset); drifting(jmp_fit) - ddrifting(jmp_fit)/2];
                jmp = [jmp_reset; jmp_fit];
                drifting = interp1(jmp, d_points, (1 : numel(drifting))', 'spline');
                
                dt_ph = drifting + dt_ph_dj;
                dt_pr = drifting + dt_pr_dj;
                
                t_offset = round(mean(dt_pr(jmp) - time_desync(jmp) + ddrifting(jmp)/2) * 1e7) * 1e-7;
                dt_ph = dt_ph - t_offset;
                dt_pr = dt_pr - t_offset;
                
                ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                pr = bsxfun(@minus, pr, dt_pr .* 299792458);
            else
                [ph_dj, dt_ph_dj] = Core_PP.remDtJumps(ph);
                [pr_dj, dt_pr_dj] = Core_PP.remDtJumps(pr);
                ddt_pr = Core_Utils.diffAndPred(dt_pr_dj);
                jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                if numel(jmp_reset) > 2
                    drifting = interp1(jmp_reset, dt_pr_dj(jmp_reset), (1 : numel(ddt_pr))', 'spline');
                else
                    drifting = 0;
                end
                dt_pr = dt_pr_dj - drifting;
                t_offset = mean(dt_pr);
                dt_ph = dt_ph_dj - drifting - t_offset;
                dt_pr = dt_pr - t_offset;

                ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                pr = bsxfun(@minus, pr, dt_pr .* 299792458);
            end
            
            if any(dt_ph_dj)
                log.addMessage(log.indent('Correcting carrier phases jumps', 6),100);
            else
                log.addMessage(log.indent('Correcting carrier phases for a dt drift estimated from desync interpolation', 6),100);
            end
            if any(dt_pr_dj)
                log.addMessage(log.indent('Correcting pseudo-ranges jumps', 6),100);
            else
                log.addMessage(log.indent('Correcting pseudo-ranges for a dt drift estimated from desync interpolation', 6),100);
            end
        end
        
        function [obs_out] = remShortArcs(obs_in, min_len)
            % Removal of observation arcs shorter than given threshold.
            % SYNTAX:
            %   [obs_out] = Core_PP.remShortArcs(obs_in, threshold_gap);
            %
            % INPUT:
            %   obs_in = input observation matrix (num_sat x epochs)
            %   min_len = threshold on the number of epochs for an arc to be
            %                   considered "short"
            %
            % OUTPUT:
            %   obs_out = input observation matrix (num_sat x epochs)
            
            if (nargin < 2)
                min_len = 10;
            end
            
            obs_out = obs_in';
            
            for s = 1 : size(obs_in,1)
                
                % find intervals of zeros
                lim = getOutliers(obs_in(s,:) == 0 | isnan(obs_in(s, :)));
                
                % find the intervals of good obs_in
                lim = [[1; lim(:,2)+1] [lim(:,1)-1; size(obs_in,2)]];
                lim = [lim (lim(:,2) - lim(:,1) +1)];
                
                lim = lim(lim(:,3) <= min_len, :);
                
                for i = 1 : size(lim,1)
                    obs_out(lim(i,1) : lim(i,2), s) = 0;
                end
            end
            
            obs_out = obs_out';
        end        
    end
    
    % ==================================================================================================================================================
    %  PRIVATE FUNCTIONS called by public calls
    % ==================================================================================================================================================
    
    methods (Access = public)
    end
    
    % ==================================================================================================================================================
    %  STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function [data_out] = highPass(data, spline_base, cut_off, filt_border)
            % Compute an high pass filter on a dataset
            % Additionally apply a spline interpolation before thee fft transformation
            %
            % SYNTAX:
            %   [filtered_data] = highPass(data, spline_base, cut_off, filt_border);
            %
            % INPUT:
            %   data            dataset to filter
            %   spline_base     spline_base for spline filtering (if 0 do not do it)
            %   cut_off         cut-off frequency for the filter
            %   filt_border     filter dimension (larger -> smoother)
            %
            % OUTPUT:
            %   filtered_data   data reduced by splines and spectral filter
            
            if nargin < 4
                filt_border = 1e100;
            end
            
            spline_base = min(round(numel(data)/2), spline_base);
            
            %%% DEBUG figure(133); clf; plot(data); hold on;
            
            % if the number of data is even interpolate one observation
            is_even = false;
            if (cut_off > 0) && (mod(numel(data),2) == 0) % is even
                is_even = true;
                x = (1 : numel(data))';
                data = [interp1(x(~isnan(data)), data(~isnan(data)), 0, 'pchip', 'extrap'); data];
            end
            data = data - (0 : numel(data)-1)' * (median(data(max(1,length(data)-4):end)) - median(data(1:min(5,length(data))))) / numel(data);
            
            n_obs = length(data);
            
            % filter low frequencies with splines to limit border effects:
            if (spline_base > 0)
                x = (1 : n_obs)';
                ws = 5;
                margin = 2 * round(ws/2);
                xi = (1 - margin : n_obs + margin)';
                
                dt_lf = [median(data(1:min(5,length(data)))) * ones(margin,1); data(~isnan(data)); ones(margin,1) * median(data(max(1,length(data)-4):end))];
                %%% DEBUG figure(133); plot(xi - 1, dt_lf);
                dt_lf = splinerMat([],medfilt_mat(dt_lf, ws), spline_base);
                
                data = data - dt_lf(margin + 1 : end - margin);
            end
            %%% DEBUG figure(133); plot(data((1 + is_even) : end));
            
            if (cut_off > 0)
                % replicate signal to limit the border effects
                data = [flipud(data); data; flipud(data)];
                
                %%% DEBUG figure(110); clf; plot(zero2nan(data));
                
                % Take fft
                fftx = fftshift(fft(data));
                
                % Calculate the number of unique points (find the zero frequency - mean of the signal)
                f0_id = ceil((length(data) + 1) / 2);
                
                % This is an evenly spaced frequency vector with f0_id points.
                f = (0:f0_id-1) ./ length(data);
                
                % FFT is symmetric, throw away second half
                fftx_half = fftx(f0_id : end);
                
                % design the filter
                w = zeros(size(fftx_half));
                f_cut = cut_off;
                id_cut = find(f > f_cut, 1, 'first');
                filt_border = min(2*id_cut - 1, filt_border);
                w(1 : (id_cut - ((filt_border - 1) / 2))) = 0;
                w((id_cut + ((filt_border - 1) / 2)) : end) = 1;
                w((id_cut - ((filt_border - 1) / 2)) : (id_cut + ((filt_border - 1) / 2))) = 1 - (1 ./ (1+exp(-10:20/(filt_border-1):10)))'; % use a sigmoid filter for the design
                
                % compute power spectra - for debug
                %%% DEBUG psd = abs(fftx_half)/length(data)*sqrt(length(data)); % Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
                %%% DEBUG figure(125); clf; loglog(f, psd)
                
                % apply the filter
                fftx_half = fftx_half .* w;
                
                % compute power spectra of the filtered data - for debug
                %%% DEBUG psd_f = abs(fftx_half) / length(data) * sqrt(length(data)); % Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
                
                fftx = ifftshift([flipud(conj(fftx_half)); fftx_half(2:end)]);
                
                %%% DEBUG figure(125); hold on; loglog(f, psd_f)
                data_out = real(ifft(fftx));
                data_out = data_out(n_obs + ((1 + is_even) : n_obs));
                %%% DEBUG figure(133); plot(data_out); title('data');
            else
                data_out = data;
            end
        end        
    end
    
end
