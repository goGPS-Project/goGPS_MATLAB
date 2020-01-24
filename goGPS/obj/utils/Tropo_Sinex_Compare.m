%   CLASS Tropo_Sinex_Compare
% =========================================================================
%
% DESCRIPTION
%   Class to compare tropo result whith external solutions
% EXAMPLE
%   tsc = Tropo_Sinex_Compare();
%


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     Andrea Gatti
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef Tropo_Sinex_Compare < handle
    
    properties
        results
        log
    end
    
    methods
        function this = Tropo_Sinex_Compare()
            this.results = struct();
            this.log = Logger.getInstance();
        end
        
        function addTropoSinexFile(this, filename, result_n)
            % add sinex file any previous data will be overwritten
            %
            % SYNTAX
            %     this.addTropoSinexFile(filename, <result_n>)
            
            if nargin < 3
                result_n = 2;
            end
            if ~isfield(this.results, ['r' num2str(result_n)])
                this.results.(['r' num2str(result_n)]) = struct();
            end
            [results] = this.parseMultiStationTropoSinex(filename);
            stas = fieldnames(results);
            if isempty(stas)
                this.log.addWarning('No result present in file')
            else
                for s = 1:length(stas)
                    data = results.(stas{s});
                    if isfield(this.results.(['r' num2str(result_n)]), stas{s})
                        
                        [this.results.(['r' num2str(result_n)]).(stas{s}).time, idx1, idx2] = this.results.(['r' num2str(result_n)]).(stas{s}).time.injectBatch(data.time);
                        this.results.(['r' num2str(result_n)]).(stas{s}).ztd     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(stas{s}).ztd, data.TROTOT , idx1, idx2);
                        this.results.(['r' num2str(result_n)]).(stas{s}).tgn     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(stas{s}).tgn, data.TGNTOT, idx1, idx2);
                        this.results.(['r' num2str(result_n)]).(stas{s}).tge     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(stas{s}).tge, data.TGETOT, idx1, idx2);
                        idx_time_coord = this.results.(['r' num2str(result_n)]).(stas{s}).coord_time == data.time.getCentralTime;
                        if sum(idx_time_coord) == 0
                            this.results.(['r' num2str(result_n)]).(stas{s}).xyz(end+1,:)     = data.xyz;
                            this.results.(['r' num2str(result_n)]).(stas{s}).coord_time.append(data.time.getCentralTime);
                        else
                            this.results.(['r' num2str(result_n)]).(stas{s}).xyz(idx_time_coord,:)     = data.xyz;
                        end
                    else
                        this.results.(['r' num2str(result_n)]).(stas{s}) = struct();
                        this.results.(['r' num2str(result_n)]).(stas{s}).time = data.time.getCopy();
                        this.results.(['r' num2str(result_n)]).(stas{s}).ztd = data.TROTOT;
                        this.results.(['r' num2str(result_n)]).(stas{s}).tgn = data.TGNTOT;
                        this.results.(['r' num2str(result_n)]).(stas{s}).tge = data.TGETOT;
                        this.results.(['r' num2str(result_n)]).(stas{s}).coord_time = data.time.getCentralTime();
                        this.results.(['r' num2str(result_n)]).(stas{s}).xyz = data.xyz;
                    end
                end
            end
        end
        
        function [missing_days] = addIGSOfficialStation(this, sta_name, time)
            % add IGS offical station
            %
            % SYNTAX
            %     this.addIGSOfficialStation(sta_name, time)
            
            if ~iscell(sta_name)
                sta_name = {sta_name};
            end
            [mjd] = floor(time.getMJD);
            missing_days = [];
            fnp = File_Name_Processor();
                                   
            data_dir = fullfile(Core.getInstallDir, '..' , 'data');
            % Try to download data using aria
            flag_aria = true;
            for d = unique(mjd)'
                if flag_aria
                    i = 0;
                    file_name_lst = {};
                    f_ext_lst = {};
                    f_status_lst = false(0);
                    c_time = GPS_Time.fromMJD(d);
                    for s = 1 : numel(sta_name)
                        file_path = fnp.dateKeyRep(sprintf('%s/station/IGS_solutions/TROPO/${YYYY}/${DOY}/%s${DOY}0.${YY}zpd', data_dir, lower(sta_name{s})), c_time);
                        if exist(file_path, 'file') ~= 2
                            [out_dir, file_name, file_ext] = fileparts(file_path);
                            i = i + 1;
                            file_name_lst(i) = {fnp.dateKeyRep(sprintf('ftp://cddis.nasa.gov/pub/gps/products/troposphere/zpd/${YYYY}/${DOY}/%s${DOY}0.${YY}zpd', lower(sta_name{s})), c_time)};
                            f_ext_lst(i) = {'.gz'};
                            f_status_lst(i) = false;
                        end
                    end
                    if ~isempty(file_name_lst)
                        f_status_lst = Core_Utils.aria2cDownloadUncompress(file_name_lst, f_ext_lst, f_status_lst, [], out_dir);
                        flag_aria = ~isempty(f_status_lst);
                    end
                end
            end
            
            for d = unique(mjd)'
                c_time = GPS_Time.fromMJD(d);
                for s = 1 : numel(sta_name)
                    file_name = fnp.dateKeyRep(sprintf('%s/station/IGS_solutions/TROPO/${YYYY}/${DOY}/%s${DOY}0.${YY}zpd', data_dir, lower(sta_name{s})), c_time);
                    % Use ftp if and only if aria is not working
                    if ~flag_aria && exist(file_name, 'file') ~= 2
                        remote_file_name = fnp.dateKeyRep(sprintf('pub/gps/products/troposphere/zpd/${YYYY}/${DOY}/%s${DOY}0.${YY}zpd.gz',lower(sta_name{s})), c_time);
                        ftp_dw = FTP_Downloader('cddis.nasa.gov', 21);
                        [pathstr, name, ext] = fileparts(file_name);
                        ftp_dw.downloadUncompress(remote_file_name, pathstr);
                    end
                    if ~(exist(file_name, 'file') == 2)
                        missing_days = [missing_days d];
                        [pathstr, name, ext] = fileparts(file_name);
                        
                        this.log.addWarning(sprintf(' File %s not found',[name, ext]));
                    else
                        this.addTropoSinexFile(file_name);
                    end
                end
            end
        end
        
        
        function loadFilesFromFolder(this, dir_path, pattern,result_n)
            % add all sinex present in a folder that match the pattern
            %
            % SYNTAX:
            %      this.loadFilesFromFolder(dir, pattern)
            if nargin < 3
                pattern = 'd*';
            end
            if nargin < 4
                result_n = 2;
            end
            files = dir(dir_path);
            for f = 1:length(files)
                if ~isempty(regexp(files(f).name,pattern))
                    this.addTropoSinexFile([dir_path filesep files(f).name],result_n);
                end
            end
        end
        
        function addReceivers(this, recs, result_n)
            % add rceeiver
            %
            % SYNTAX
            %     this.addTropoSinexFile(filename, <result_n>)
            
            if nargin < 3
                result_n = 1;
            end
            this.results.(['r' num2str(result_n)]) = struct();
            for s = 1:length(recs)
                sta = recs(s).getMarkerName4Ch;
                if ~strcmpi(sta,'unkn')
                    data = recs(s).out;
                    if isfield(this.results.(['r' num2str(result_n)]), sta)
                        
                        [this.results.(['r' num2str(result_n)]).(sta).time, idx1, idx2] = this.results.(['r' num2str(result_n)]).(sta).time.injectBatch(data.time);
                        this.results.(['r' num2str(result_n)]).(sta).ztd     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(sta).ztd, data.ztd , idx1, idx2);
                        this.results.(['r' num2str(result_n)]).(sta).tgn     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(sta).tgn, data.tgn, idx1, idx2);
                        this.results.(['r' num2str(result_n)]).(sta).tge     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(sta).tge, data.tge, idx1, idx2);
                        
                        [this.results.(['r' num2str(result_n)]).(sta).coord_time, idx1, idx2] = this.results.(['r' num2str(result_n)]).(sta).coord_time.injectBatch(data.time_pos);
                        this.results.(['r' num2str(result_n)]).(sta).xyz     = Core_Utils.injectData(this.results.(['r' num2str(result_n)]).(sta).xyz, data.xyz, idx1, idx2);
                    else
                        this.results.(['r' num2str(result_n)]).(sta) = struct();
                        this.results.(['r' num2str(result_n)]).(sta).xyz = data.xyz;
                        this.results.(['r' num2str(result_n)]).(sta).coord_time = data.time_pos;
                        this.results.(['r' num2str(result_n)]).(sta).time = data.time;
                        this.results.(['r' num2str(result_n)]).(sta).ztd = data.ztd;
                        this.results.(['r' num2str(result_n)]).(sta).tgn = data.tgn;
                        this.results.(['r' num2str(result_n)]).(sta).tge = data.tge;
                        this.results.(['r' num2str(result_n)]).(sta).delta_enu = [data.parent.ant_delta_en data.parent.ant_delta_h];
                    end
                end
            end
        end
        
        function plotTropoDifference(this,mode)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            if nargin < 2
                mode = 1;
            end
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            for s1 = 1 : length(sta1)
                for s2 = 1 : length(sta2)
                    if strcmpi(sta1{s1},sta2{s2})
                        data1 = this.results.r1.(sta1{s1});
                        data2 = this.results.r2.(sta2{s2});
                        f = figure; f.Name = sprintf('%03d: %s ', f.Number, sta1{s1}); f.NumberTitle = 'off';
                        % plot ZTD series
                        subplot(3,4,1:3)
                        diffint = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.ztd), data1.time.getMatlabTime, zero2nan(data1.ztd),'interpolate');
                        if any(~isnan(diffint))
                            plot(data1.time.getMatlabTime,diffint,'.','Color','b');
                            hold on;
                            [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.ztd), data1.time.getMatlabTime, zero2nan(data1.ztd),'spline');
                            
                            plot(data2.time.getMatlabTime,diffagg,'Color','r');
                            stdd = nan_std(diffint);
                            ylim([-4*stdd 4*stdd])
                            xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                            setTimeTicks(5,'yyyy/mm/dd');
                            title('ZTD')
                            subplot(3,4,4)
                            hist(noNaN(diffagg),30)
                            rms = mean(abs(noNaN(diffagg)*1e3));
                            bias = mean(noNaN(diffagg)*1e3);
                            title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                            xlabel('Aggregated differences');
                            ylabel('#');
                            
                            if mode == 1
                                subplot(3,4,5:7)
                                diffint = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tge), data1.time.getMatlabTime, zero2nan(data1.tge),'interpolate');
                                plot(data1.time.getMatlabTime,diffint,'.','Color','b');
                                hold on;
                                [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tge), data1.time.getMatlabTime, zero2nan(data1.tge),'aggregate');
                                plot(data2.time.getMatlabTime,diffagg,'Color','r');
                                stdd = nan_std(diffint);
                                ylim([-4*stdd 4*stdd])
                                xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                                setTimeTicks(5,'yyyy/mm/dd');
                                title('East gradient')
                                subplot(3,4,8)
                                hist(noNaN(diffagg),30)
                                rms = mean(abs(noNaN(diffagg)*1e3));
                                bias = mean(noNaN(diffagg)*1e3);
                                title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                                xlabel('Aggregated differences');
                                ylabel('#');
                                % plot gn series
                                subplot(3,4,9:11)
                                diffint = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tgn), data1.time.getMatlabTime, zero2nan(data1.tgn),'interpolate');
                                plot(data1.time.getMatlabTime,diffint,'.','Color','b');
                                hold on;
                                [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tgn), data1.time.getMatlabTime, zero2nan(data1.tgn),'aggregate');
                                plot(data2.time.getMatlabTime,diffagg,'Color','r');
                                stdd = nan_std(diffint);
                                ylim([-4*stdd 4*stdd])
                                xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                                setTimeTicks(5,'yyyy/mm/dd');
                                title('North gradient')
                                subplot(3,4,12)
                                hist(noNaN(diffagg),30)
                                rms = mean(abs(noNaN(diffagg)*1e3));
                                bias = mean(noNaN(diffagg)*1e3);
                                title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                                xlabel('Aggregated differences');
                                ylabel('#');
                            else
                                % POLAR COORDINATES
                                subplot(3,4,5:7)
                                diffintE = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tge), data1.time.getMatlabTime, zero2nan(data1.tge),'interpolate');
                                [diffaggE] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tge), data1.time.getMatlabTime, zero2nan(data1.tge),'aggregate');
                                diffintN = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tgn), data1.time.getMatlabTime, zero2nan(data1.tgn),'interpolate');
                                [diffaggN] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tgn), data1.time.getMatlabTime, zero2nan(data1.tgn),'aggregate');
                                subplot(3,4,5:7)
                                diffint = sqrt(diffintE.^2 + diffintN.^2);
                                plot(data1.time.getMatlabTime,diffint,'.','Color','b');
                                hold on;
                                diffagg = sqrt(diffaggE.^2 + diffaggN.^2);
                                plot(data2.time.getMatlabTime,diffagg,'Color','r');
                                stdd = nan_std(diffint);
                                ylim([0 7*stdd])
                                xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                                setTimeTicks(5,'yyyy/mm/dd');
                                title('Gradient magnitude')
                                subplot(3,4,8)
                                hist(noNaN(diffagg),30)
                                rms = mean(abs(noNaN(diffagg)*1e3));
                                bias = mean(noNaN(diffagg)*1e3);
                                title(sprintf('RMS: %0.2f mm',rms))
                                xlabel('Aggregated differences');
                                ylabel('#');
                                % plot gn series
                                subplot(3,4,9:11)
                                diffint = atan2(diffintE, diffintN)/pi*180; diffint(abs(diffint) > 180) = sign(diffint(abs(diffint) > 180)).*(360 - abs(diffint(abs(diffint) > 180)));
                                
                                plot(data1.time.getMatlabTime,diffint,'.','Color','b');
                                hold on;
                                diffagg = atan2(diffaggE, diffaggN)/pi*180; diffagg(abs(diffagg) > 180) = sign(diffagg(abs(diffagg) > 180)).*(360 - abs(diffagg(abs(diffagg) > 180)));
                                plot(data2.time.getMatlabTime,diffagg,'Color','r');
                                stdd = nan_std(diffint);
                                ylim([-180 180])
                                xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                                setTimeTicks(5,'yyyy/mm/dd');
                                title('Gradient Azimuth')
                                magn = sqrt(data2.tge.^2 + data2.tgn.^2);
                                diffagg = diffagg(magn > 3*1e-3);
                                subplot(3,4,12)
                                hist(noNaN(diffagg),30)
                                rms = mean(abs(noNaN(diffagg)));
                                bias = mean(noNaN(diffagg));
                                title(sprintf('Bias: %0.2f degree RMS: %0.2f degree',bias,rms))
                                xlabel('Aggregated differences');
                                ylabel('#');
                            end
                        end
                    end
                end
            end
            
        end
        
        function plotTropoDifferenceSummary(this)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            [Csta,ia,ib] = intersect(lower(sta1),lower(sta2), 'stable');
            siz_col = ceil(sqrt(length(Csta)));
            siz_row = ceil(length(Csta)/siz_col);
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'ZTD Differences'); f.NumberTitle = 'off';
            xlim1 = -30;
            xlim2 = 40;
            dx = 2.5;
            for s = 1 : length(Csta)
                s1 = ia(s);
                s2 = ib(s);
                data1 = this.results.r1.(sta1{s1});
                data2 = this.results.r2.(sta2{s2});
                % plot ZTD series
                subplot(siz_col,siz_row,s)
                
                [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.ztd), data1.time.getMatlabTime, zero2nan(data1.ztd),'aggregate');
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Aggregated differences');
                ylabel('#');
            end
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'East Gradient Differences'); f.NumberTitle = 'off';
            xlim1 = -5;
            xlim2 = 5;
            dx = 0.5;
            for s = 1 : length(Csta)
                s1 = ia(s);
                s2 = ib(s);
                data1 = this.results.r1.(sta1{s1});
                data2 = this.results.r2.(sta2{s2});
                
                [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tge), data1.time.getMatlabTime, zero2nan(data1.tge),'aggregate');
                subplot(siz_col,siz_row,s)
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Aggregated differences');
                ylabel('#');
            end
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'North Gradient Differences'); f.NumberTitle = 'off';
            xlim1 = -5;
            xlim2 = 5;
            dx = 0.5;
            for s = 1 : length(Csta)
                s1 = ia(s);
                s2 = ib(s);
                data1 = this.results.r1.(sta1{s1});
                [diffagg] = timeSeriesComparison(data2.time.getMatlabTime, zero2nan(data2.tgn), data1.time.getMatlabTime, zero2nan(data1.tgn),'aggregate');
                subplot(siz_col,siz_row,s)
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Aggregated differences');
                ylabel('#');
            end
        end
        
        
        function plotTropoComparison(this)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            for s1 = 1 : length(sta1)
                for s2 = 1 : length(sta2)
                    if strcmpi(sta1{s1},sta2{s2})
                        data1 = this.results.r1.(sta1{s1});
                        data2 = this.results.r2.(sta2{s2});
                        f = figure; f.Name = sprintf('%03d: %s ', f.Number, sta1{s1}); f.NumberTitle = 'off';
                        % plot ZTD series
                        subplot(3,1,1)
                        plot(data1.time.getMatlabTime,zero2nan(data1.ztd),'.','Color','b');
                        hold on;
                        plot(data2.time.getMatlabTime,zero2nan(data2.ztd),'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        
                        
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        title('ZTD')
                        % plot ge series
                        subplot(3,1,2)
                        plot(data1.time.getMatlabTime,zero2nan(data1.tge),'.','Color','b');
                        hold on;
                        plot(data2.time.getMatlabTime,zero2nan(data2.tge),'r');
                        %                     diffagg = timeSeriesComparison(data2.time.getMatlabTime,data2.tge,data1.time.getMatlabTime,data1.tge,'aggregate');
                        %                     plot(data2.time.getMatlabTime,diffagg,'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        %                     subplot(3,4,11)
                        %                     hist(diffint);
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        title('East gradient')
                        % plot gn series
                        subplot(3,1,3)
                        plot(data1.time.getMatlabTime,zero2nan(data1.tgn),'.','Color','b');
                        hold on;
                        plot(data2.time.getMatlabTime,zero2nan(data2.tgn),'r');
                        %                     diffagg = timeSeriesComparison(data2.time.getMatlabTime,data2.tgn,data1.time.getMatlabTime,data1.tgn,'aggregate');
                        %                     plot(data2.time.getMatlabTime,diffagg,'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        %                     subplot(3,4,12)
                        %                     hist(diffint);
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        title('North gradient')
                    end
                end
            end
            
        end
        
        function plotCooComparison(this)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            for s1 = 1 : length(sta1)
                for s2 = 1 : length(sta2)
                    if strcmpi(sta1{s1},sta2{s2})
                        data1 = this.results.r1.(sta1{s1});
                        data2 = this.results.r2.(sta2{s2});
                        f = figure; f.Name = sprintf('%03d: %s ', f.Number, sta1{s1}); f.NumberTitle = 'off';
                        xyz_ref = median(data2.xyz,1);
                        enu1 = Coordinates.cart2local(xyz_ref, data1.xyz - repmat(xyz_ref, size(data1.xyz,1),1) )*1e3;
                        enu2 = Coordinates.cart2local(xyz_ref, data2.xyz  - repmat(xyz_ref, size(data2.xyz,1),1))*1e3;
                        % plot E
                        subplot(3,1,1)
                        [LIA,LocB] = ismembertol(data1.coord_time.getMatlabTime, data2.coord_time.getMatlabTime,1e-8);
                        LocB= LocB(LocB ~= 0);
                        plot(data1.coord_time.getMatlabTime,zero2nan(enu1(:,1)),'.','Color','b');
                        hold on;
                        plot(data2.coord_time.getMatlabTime,zero2nan(enu2(:,1)),'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        
                        ylabel('mm')
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        std1 = std(enu1(LIA,1));
                        std2 = std(enu2(LocB,1));
                        title(sprintf('East std1: %.2f std2: %.2f [mm]',std1,std2))
                        % plot N
                        subplot(3,1,2)
                        plot(data1.coord_time.getMatlabTime,zero2nan(enu1(:,2)),'.','Color','b');
                        hold on;
                        plot(data2.coord_time.getMatlabTime,zero2nan(enu2(:,2)),'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        ylabel('mm')
                        %                     subplot(3,4,11)
                        %                     hist(diffint);
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        std1 = std(enu1(LIA,2));
                        std2 = std(enu2(LocB,2));
                        title(sprintf('North std1: %.2f std2: %.2f [mm]',std1,std2))
                        % plot U
                        subplot(3,1,3)
                        plot(data1.coord_time.getMatlabTime,zero2nan(enu1(:,3)),'.','Color','b');
                        hold on;
                        plot(data2.coord_time.getMatlabTime,zero2nan(enu2(:,3)),'r');
                        setTimeTicks(5,'yyyy/mm/dd');
                        ylabel('mm')
                        %                     subplot(3,4,12)
                        %                     hist(diffint);
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        std1 = std(enu1(LIA,3));
                        std2 = std(enu2(LocB,3));
                        title(sprintf('Up std1: %.2f std2: %.2f [mm]',std1,std2))
                    end
                end
            end
            
        end
        
        
        function plotCooDifference(this,mode)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            if nargin < 2
                mode = 1;
            end
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            for s1 = 1 : length(sta1)
                for s2 = 1 : length(sta2)
                    if strcmpi(sta1{s1},sta2{s2})
                        data1 = this.results.r1.(sta1{s1});
                        data2 = this.results.r2.(sta2{s2});
                        xyz_ref = median(data2.xyz,1);
                        if mode == 1
                            enu1 = Coordinates.cart2local(xyz_ref, data1.xyz - repmat(xyz_ref, size(data1.xyz,1),1) );% - repmat(data1.delta_enu,size( data1.xyz,1),1);
                            enu2 = Coordinates.cart2local(xyz_ref, data2.xyz  - repmat(xyz_ref, size(data2.xyz,1),1));
                        else
                            enu1 = data1.xyz;
                            enu2 = data2.xyz;
                        end
                        [LIA,LocB] = ismembertol(data1.coord_time.getMatlabTime,data2.coord_time.getMatlabTime,1/24, 'DataScale', 1);
                        diff_enu = enu1(LIA,:) - enu2(LocB(LocB~=0),:);
                        disp(sprintf('%s : %.4f %.4f %.4f , %.4f ',sta1{s1}, mean(diff_enu,1) ,sqrt(sum(mean(diff_enu,1).^2))))
                        
                        
                        f = figure; f.Name = sprintf('%03d: %s ', f.Number, sta1{s1}); f.NumberTitle = 'off';
                        % plot ZTD series
                        subplot(3,4,1:3)
                        diffagg = diff_enu(:,1);
                        plot(data1.coord_time.getEpoch(LIA).getMatlabTime,diffagg,'.','Color','b');
                        stdd = nan_std(diffagg);
                        %ylim([-4*stdd 4*stdd])
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        setTimeTicks(5,'yyyy/mm/dd');
                        if mode == 1
                            title('East')
                        else
                            title('X')
                        end
                        subplot(3,4,4)
                        hist(noNaN(diffagg),30)
                        rms = mean(abs(noNaN(diffagg)*1e3));
                        bias = mean(noNaN(diffagg)*1e3);
                        title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                        xlabel('Differences');
                        ylabel('#');
                        
                        subplot(3,4,5:7)
                        diffagg = diff_enu(:,2);
                        plot(data1.coord_time.getEpoch(LIA).getMatlabTime,diffagg,'.','Color','b');
                        stdd = nan_std(diffagg);
                        %ylim([-4*stdd 4*stdd])
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        setTimeTicks(5,'yyyy/mm/dd');
                        setTimeTicks(5,'yyyy/mm/dd');
                        if mode == 1
                            title('North')
                        else
                            title('Y')
                        end
                        subplot(3,4,8)
                        hist(noNaN(diffagg),30)
                        rms = mean(abs(noNaN(diffagg)*1e3));
                        bias = mean(noNaN(diffagg)*1e3);
                        title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                        xlabel('Differences');
                        ylabel('#');
                        % plot gn series
                        subplot(3,4,9:11)
                        diffagg = diff_enu(:,3);
                        plot(data1.coord_time.getEpoch(LIA).getMatlabTime,diffagg,'.','Color','b');
                        stdd = nan_std(diffagg);
                        %ylim([-4*stdd 4*stdd])
                        xlim([data1.time.first.getMatlabTime data1.time.last.getMatlabTime])
                        setTimeTicks(5,'yyyy/mm/dd');
                        setTimeTicks(5,'yyyy/mm/dd');
                        if mode == 1
                            title('Up')
                        else
                            title('Z')
                        end
                        subplot(3,4,12)
                        hist(noNaN(diffagg),30)
                        rms = mean(abs(noNaN(diffagg)*1e3));
                        bias = mean(noNaN(diffagg)*1e3);
                        title(sprintf('Bias: %0.2f mm RMS: %0.2f mm',bias,rms))
                        xlabel('Differences');
                        ylabel('#');
                        
                    end
                end
            end
        end
        
        function [lons] = getLon(this)
            % get longitude from sinex results in degree
            %
            % SYNTAX
            % [lon] = this.getLon()
            stas = fieldnames(this.results.r2);
            n_rec = length(stas);
            lons = nan(n_rec,1);
            for i = 1 : n_rec
                [~,lons(i)] = Coordinates.cart2geod(this.results.r2.(stas{i}).xyz(1,:));
                lons(i) = lons(i)/pi*180;
            end
        end
        
        function [lats] = getLat(this)
            % get longitude from sinex results in degree
            %
            % SYNTAX
            % [lon] = this.getLon()
            stas = fieldnames(this.results.r2);
            n_rec = length(stas);
            lats = nan(n_rec,1);
            for i = 1 : n_rec
                [lats(i),~] = Coordinates.cart2geod(this.results.r2.(stas{i}).xyz(1,:));
                lats(i) = lats(i)/pi*180;
            end
        end
        
        function [h] = getHeightOrtho(this)
            % get longitude from sinex results
            %
            % SYNTAX
            % [lon] = this.getLon()
            stas = fieldnames(this.results.r2);
            n_rec = length(stas);
            h = nan(n_rec,1);
            for i = 1 : n_rec
                [lat,lon,h(i)] = Coordinates.cart2geod(this.results.r2.(stas{i}).xyz(1,:));
                h(i) = h(i) - Coordinates.getOrthometricCorrFromLatLon(lat, lon);
            end
        end
        
        function n = getNumberSinex(this)
            % get number of sinex results
            %
            % SYNTAX
            % [n] = this.getNumberSinex()
            n = length(fieldnames(this.results.r2));
        end
        
        function [ztd,time] = getZtdSinex(this,s)
            % get ztd and time of the sth station
            %
            % SYNTAX
            %  [ztd,time] = this.getZtdSinex(s)
            stas = fieldnames(this.results.r2);
            ztd = this.results.r2.(stas{s}).ztd;
            time =  this.results.r2.(stas{s}).time;
        end
        
        function [name] = getName(this,s)
            % get nathe of nth station
            %
            % SYNTAX
            %   [name] = this.getName(s)
            stas = fieldnames(this.results.r2);
            name = stas{s};
        end
        
        function plotCooDifferenceSummary(this)
            % print the difference between the results
            %
            % SYNTAX
            %   this.plotComparison()
            if nargin < 2
                mode = 1;
            end
            sta1 = fieldnames(this.results.r1);
            sta2 = fieldnames(this.results.r2);
            [Csta,ia,ib] = intersect(lower(sta1),lower(sta2), 'stable');
            siz_col = ceil(sqrt(length(Csta)));
            siz_row = ceil(length(Csta)/siz_col);
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'EAST coordinate differences'); f.NumberTitle = 'off';
            diffs = {};
            xlim1 = -20;
            xlim2 = 20;
            dx = 2.5;
            for s = 1 : length(Csta)
                s1 = ia(s);
                s2 = ib(s);
                data1 = this.results.r1.(sta1{s1});
                data2 = this.results.r2.(sta2{s2});
                xyz_ref = median(data2.xyz,1);
                if mode == 1
                    enu1 = Coordinates.cart2local(xyz_ref, data1.xyz - repmat(xyz_ref, size(data1.xyz,1),1));% - repmat(data1.delta_enu,size( data1.xyz,1),1);
                    enu2 = Coordinates.cart2local(xyz_ref, data2.xyz  - repmat(xyz_ref, size(data2.xyz,1),1));
                else
                    enu1 = data1.xyz;
                    enu2 = data2.xyz;
                end
                [LIA,LocB] = ismembertol(data1.coord_time.getMatlabTime,data2.coord_time.getMatlabTime,1/24, 'DataScale', 1);
                diff_enu = enu1(LIA,:) - enu2(LocB(LocB~=0),:);
                %disp(sprintf('%s : %.4f %.4f %.4f , %.4f ',sta1{s1}, mean(diff_enu,1) ,sqrt(sum(mean(diff_enu,1).^2))))
                diffs{s} = diff_enu;
                
                
                % plot ZTD series
                diffagg = diff_enu(:,1);
                subplot(siz_row,siz_col,s)
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Diff [mm]');
                ylabel('#');
                xlim([xlim1 xlim2]);
            end
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'NORTH coordinate differences'); f.NumberTitle = 'off';
            for s = 1 : length(Csta)
                % plot ZTD series
                diff_enu = diffs{s};
                s1 = ia(s);
                diffagg = diff_enu(:,2);
                subplot(siz_row,siz_col,s)
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Differences [mm]');
                ylabel('#');
                xlim([xlim1 xlim2]);
            end
            f = figure; f.Name = sprintf('%03d: %s ', f.Number, 'UP coordinate differences'); f.NumberTitle = 'off';
            for s = 1 : length(Csta)
                % plot ZTD series
                diff_enu = diffs{s};
                s1 = ia(s);
                diffagg = diff_enu(:,3);
                subplot(siz_row,siz_col,s)
                hist(noNaN(diffagg*1e3),[xlim1:dx:xlim2])
                rms = mean(abs(noNaN(diffagg)*1e3));
                bias = mean(noNaN(diffagg)*1e3);
                title(sprintf('%s \n Bias: %0.2f mm \n RMS: %0.2f mm',upper(sta1{s1}),bias,rms))
                xlabel('Differences [mm]');
                ylabel('#');
                xlim([xlim1 xlim2]);
            end
        end
        
    end
    methods (Static)
        function [results] = parseMultiStationTropoSinex(filename)
            % parse a sinex fil containing tropopsheric products (e.g. the ones produced by epn)
            %
            % SYNTAX:
            %  [results] = Tropo_Sinex_Compare.parseMultiStationTropoSinex(filename)
            fid = fopen(filename);
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            % get starting block lines
            st_idxes = find(txt(lim(:,1)) == '+');
            end_idxes = find(txt(lim(:,1)) == '-');
            for i = 1:length(st_idxes)
                st_idx = st_idxes(i);
                end_idx = end_idxes(i);
                if strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/DESCRIPTION')
                    for l = st_idx : end_idx
                        if strfind(txt(lim(l,1): lim(l,2)),' SOLUTION_FIELDS_1')
                            pars = strsplit(strtrim(txt(lim(l,1): lim(l,2))));
                            pars(1) = [];
                            n_par = length(pars);
                            for i = 1: n_par
                                pars{i} = strrep(pars{i},'#','NUM_');
                                all_res.(pars{i}) = [];
                            end
                        end
                    end
                    
                elseif strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/STA_COORDINATES')
                    n_lin = end_idx - st_idx -2;
                    sta_4char = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin ,1));
                    xyz = reshape(sscanf( (txt(repmat(lim(st_idx+2:end_idx-1,1),1,39) + repmat([16:54], n_lin,1)))','%f %f %f\n'),3,n_lin)';
                elseif strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/SOLUTION')
                    n_lin = end_idx - st_idx -2;
                    sta_4char_trp = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin,1));
                    obs_lines = st_idx+2:end_idx-1;
                    idx_no_ast = txt(lim(obs_lines,1)) ~= '*';
                    idx_no_ast = obs_lines(idx_no_ast);
                    n_lin = length(idx_no_ast);
                    sta_4char_trp = txt(repmat(lim(idx_no_ast,1),1,4) + repmat([1:4],n_lin,1));
                    end_col = min(lim(idx_no_ast,3));
                    
                    data = reshape(sscanf( (txt(repmat(lim(idx_no_ast,1),1,end_col) + repmat([0 :(end_col-1)],n_lin,1)))',['%*s %f:%f:%f' repmat(' %f',1,n_par)] ),n_par+3,n_lin)';
                    year = data(:,1);
                    idx_70 = year<70;
                    year(idx_70) = year(idx_70) + 2000;
                    year(~idx_70) = year(~idx_70) + 1900;
                    doy = data(:,2);
                    sod = data(:,3);
                    for i = 1 : n_par
                        all_res.(pars{i}) = [all_res.(pars{i}) ;  data(:,3+i)/1e3]; % -> covert in meters
                    end
                end
                
            end
            for s = 1: size(sta_4char,1)
                c_sta_4char = sta_4char(s,:);
                idx_sta  = strLineMatch(sta_4char_trp,c_sta_4char);
                if sum(idx_sta) > 0
                    for i = 1 : n_par
                        results.(c_sta_4char).(pars{i}) = all_res.(pars{i})(idx_sta);
                    end
                    results.(c_sta_4char).time = GPS_Time.fromDoySod(year(idx_sta),doy(idx_sta),sod(idx_sta));
                    results.(c_sta_4char).xyz = xyz(s,:);
                end
            end
        end
        
    end
    
end

