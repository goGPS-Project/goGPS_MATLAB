classdef Radiosonde < handle
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
    %  Written by:       Alice Bonfiglio, Stefano Barindelli
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
    
    properties
        time     % time as datetime                  datetime [n_records x 1]
        p        % pressure [hPa]                    double   [n_records x 1]
        h        % height [m]                        double   [n_records x 1]
        T        % temperature [ÿC]                  double   [n_records x 1]
        rh       % relatibe humidity [%]             double   [n_records x 1]
        pwv      % precipitable water vapor [mm]     double   [n_launches x 1]
    end
    
    properties (Access = private)
        log     % logger
    end
    
    methods
        % Creator
        function this = Radiosonde(region, plot_type, year, month, sta_numb) %valid for one station for one month
            % Core object creator
            this.log = Logger.getInstance();
            this.reset();
            if nargin < 5
                this.log.addMessage(sprintf('Error downloading the files, not enough input arguments'));
            else
                this.log.addMessage(sprintf('Downloading files'));
            end
            
            [rds_pwv, rds] = this.RADownloader(region, plot_type, year, month, sta_numb);
        end
    end
    
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        function reset(this)
            this.time = [];
            this.p  = [];
            this.h  = [];
            this.T  = [];
            this.rh  = [];
            this.pwv  = [];
        end
        
        function [rds_pwv, rds]= RADownloader(this,region, plot_type, year, month, sta_numb)
            
            y=clock;
            %sta_numb_cell=cellstr(sta_numb);
            sta_numb_struct=horzcat(repmat(['s_'],size(sta_numb,1),1),sta_numb);
            j=1;
            for t=1:size(month,1)
                tmp=1:size(sta_numb,1)*size(month,1);
                if strcmp(month(t,:),'01')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'02')
                    if leapyear(y(1,1))
                        from = '0100';
                        to = '2812';
                    else
                        from = '0100';
                        to = '2912';
                    end
                elseif strcmp(month(t,:),'03')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'04')
                    from = '0100';
                    to = '3012';
                elseif strcmp(month(t,:),'05')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'06')
                    from = '0100';
                    to = '3012';
                elseif strcmp(month(t,:),'07')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'08')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'09')
                    from = '0100';
                    to = '3012';
                elseif strcmp(month(t,:),'10')
                    from = '0100';
                    to = '3112';
                elseif strcmp(month(t,:),'11')
                    from = '0100';
                    to = '3012';
                elseif strcmp(month(t,:),'12')
                    from = '0100';
                    to = '3112';
                end
                
                options = weboptions;
                options.Timeout = 50;
                %plot_type = 'TEXT';
                address = ['http://weather.uwyo.edu/cgi-bin/sounding?region=' region '&TYPE=' plot_type '%3ALIST&YEAR=' year '&MONTH=' month(t,:) '&FROM=' from '&TO=' to '&STNM=' sta_numb(:,:)];
                char_array = webread(address);
                char_array = regexprep(char_array,'<script.*?/script>','');
                char_array = regexprep(char_array,'<style.*?/style>','');
                char_array = regexprep(char_array,'<.*?>','');
                
                if isempty(strfind(char_array,'Observations'))
                    warning('no data available :-(')
                    rds_pwv(1,j).pwv = [];
                    rds_pwv(1,j).datetime = [];
                    rds_pwv(1,j).sta_number = sta_numb;
                else
                    
                    %pwv
                    str_pw = 'Precipitable water [mm] for entire sounding:';
                    str_idx_pw = strfind(char_array,str_pw);
                    for i=1:length(str_idx_pw)
                        pw_vec(i,:) = char_array(1, str_idx_pw(i)+length(str_pw)+1 : str_idx_pw(i)+length(str_pw)+5);
                        pw_vec_d(i,:)=str2double(pw_vec(i,:));
                    end
                    
                    %time
                    str_obs_time = 'Observation time:';
                    str_idx_time = strfind(char_array,str_obs_time);
                    
                    for i=1:length(str_idx_time)
                        time_vec(i,:)=[char_array(1,str_idx_time(i)+length(str_obs_time)+1:str_idx_time(i)+length(str_obs_time)+7) char_array(1,str_idx_time(i)+length(str_obs_time)+8:str_idx_time(i)+length(str_obs_time)+11) '00'];
                    end
                    datetime_vec = datetime(time_vec,'InputFormat','yyMMdd/HHmmSS');
                    
                    str_idx_header = strfind(char_array,'-----------------------------------------------------------------------------');
                    str_idx_header=str_idx_header(2:2:end);
                    str_idx_header_fin = strfind(char_array,'Station information and sounding indices');
                    for i=1:length(str_idx_header)
                        temp=char_array(str_idx_header(i)+79:str_idx_header_fin(i)-2);
                        data(:,i) = textscan(temp,'%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%[^\n\r]','delimiter','','multipleDelimsAsOne',false,'TreatAsEmpty',{'[]'},'EmptyValue',NaN,'whitespace','','CollectOutput',true);
                        %data(:,i) = textscan(temp,'%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%[^\n\r]','delimiter','','multipleDelimsAsOne',false,'EmptyValue',NaN,'CollectOutput',true);
                        %pressure (hpa)
                        rds(t,1).(sta_numb_struct(:,:))(i,1).pres=str2double(data{1,i}(:,1));
                        %height (m)
                        rds(t,1).(sta_numb_struct(:,:))(i,1).height=str2double(data{1,i}(:,2));
                        %temperature (Celsius degree)
                        rds(t,1).(sta_numb_struct(:,:))(i,1).temp=str2double(data{1,i}(:,3));
                        %rel humidity (%)
                        rds(t,1).(sta_numb_struct(:,:))(i,1).relh=str2double(data{1,i}(:,5));
                        %radiosonde(i,1).name=sta_numb(s,:);
                        rds(t,1).(sta_numb_struct(:,:))(i,1).datetime=datetime_vec(i,:);
                        
                        this.time = [this.time; repmat(rds(t,1).(sta_numb_struct(1,:))(i,1).datetime,length(rds(t,1).(sta_numb_struct(1,:))(i,1).pres),1)];
                        this.p = [this.p; rds(t,1).(sta_numb_struct(1,:))(i,1).pres];
                        this.h = [this.h; rds(t,1).(sta_numb_struct(1,:))(i,1).height];
                        this.T = [this.T; rds(t,1).(sta_numb_struct(1,:))(i,1).temp];
                        this.rh = [this.rh; rds(t,1).(sta_numb_struct(1,:))(i,1).relh];
                    end
                    
                    rds_pwv(1,j).pwv = pw_vec_d;
                    rds_pwv(1,j).datetime = datetime_vec;
                    rds_pwv(1,j).sta_number = sta_numb(:,:);
                    this.pwv = [this.pwv; rds_pwv(t,1).pwv];
                    clear datetime_vec pw_vec_d time_vec
                end
                j=j+1;
            end
        end
    end
end