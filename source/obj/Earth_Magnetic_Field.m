% This class containing the spherical harmonics model for earth
% magnetic field

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Giulio Tagliaferro,  Andrea Gatti, Eugenio Realini
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, Andrea Gatti
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

classdef Earth_Magnetic_Field < handle
    properties
        H; % gauss coefficients
        G; % gauss coefficients
        years;
    end
    
    properties (Access = private)
        cache = false;
        state;
        c_doy;
        c_year;
        Hcache;
        Gcache;
        P;
        dP;
        P_d_length = 0.001;
    end
    
    methods
        function this = Earth_Magnetic_Field()
            this.state = Core.getCurrentSettings();
            this.importIGRFModel(this.state.getIgrfFile());
            load(File_Name_Processor.checkPath([this.state.igrf_dir filesep 'P_dP_schimdt.mat']));
            this.P = P;
            this.dP = dP;
        end
        
        function importIGRFModel(this, fname)
            % IMPORTANT TODO - tranform secular varaitions in coefficients
            fid = fopen(fname);
            
            tline = fgetl(fid);
            while ischar(tline)
                if tline(1) == '#'
                elseif  strcmp(tline(1:3),'c/s')
                    year_type = strsplit(tline,' ');
                    year_type(1:3) = [];
                    n_ep = length(year_type);
                    H = zeros(n_ep, 14, 14);
                    G = zeros(n_ep, 14, 14);
                    year_type = strcmp(year_type,'SV');
                    n_year = sum(year_type == 0);
                    n_sv = sum(year_type == 1);
                    
                elseif  strcmp(tline(1:3),'g/h')
                    tline(1:8) = [];
                    tline =strsplit(tline,' ');
                    years = cell2mat(tline(year_type == 0));
                    sv = zeros(n_sv,2);
                    sv_idx = find(year_type == 1);
                    for i = 1 : n_sv
                        str = strtrim(tline{sv_idx(i)});
                        sv(i,1) = str2num(str(1:4));
                        sv(i,2) = str2num(str([1:2 6:7]));
                    end
                elseif  strcmp(tline(1:2),'g ') || strcmp(tline(1:2),'h ')
                    n = str2num(tline(3:4));
                    m = str2num(tline(6:7));
                    vals = strsplit(tline(8:end));
                    if strcmp(vals{1},'')
                        vals(1) = [];
                    end
                    vals = str2double(vals);
                    
                    l_idx = sub2ind(size(H),[1:n_ep]',repmat(n+1,n_ep,1),repmat(m+1,n_ep,1));
                    if tline(1) == 'g'
                        G(l_idx) = vals;
                    else
                        H(l_idx) = vals;
                    end
                end
                tline = fgetl(fid);
            end
            y_idx = 1:n_ep;
            y_idx = y_idx ~= sv_idx;
            years_n = zeros(size(y_idx));
            
            years_n(y_idx) = sscanf(reshape(years,6,length(years)/6),'%6f');
            years_n(sv_idx) = sv(2);
            n_y = sv(2) -sv(1);
            idx_ref = years_n == sv(1);
            H(sv_idx,:,:) = H(idx_ref,:,:) + n_y .* H(sv_idx,:,:);
            G(sv_idx,:,:) = G(idx_ref,:,:) + n_y .* G(sv_idx,:,:);
            fclose(fid);
            this.H = permute(H,[3 2 1]);
            this.G = permute(G,[3 2 1]);
            this.years = years_n;
            
        end
        
        function V = getV(this, gps_time, r, phi, theta)
            % interpolate H at G at presetn epoch
            [y, doy] = gps_time.getDOY;
            re = GPS_SS.ELL_A/1000;
            iy = max(min(floor(y-1900)/5+1,24),1);
            int_time = ((y - this.years(iy))*365.25 + doy-1)/(365.25*5);
            H = this.H(:,:,iy) * (1 - int_time) + int_time *  this.H(:,:,iy+1);
            G = this.G(:,:,iy) * (1 - int_time) + int_time *  this.G(:,:,iy+1);
            n = size(H,1);
            P = zeros(n,n);
            % co to colatitude
            theta = pi/2 - theta;
            for i = 0:n-1;
                P(1:i+1,i+1) = legendre(i,cos(theta),'sch');
            end
            mphi = repmat((0:n-1)',1,n)*phi;
            cosm = cos(mphi); %some unnecesaary calculations
            sinm = sin(mphi); %some unnecesaary calculations
            arn = repmat((re/r).^(1:n),n,1);
            V =  re * sum(sum(arn .* (G .* cosm + H .* sinm) .* P));
        end
        
        function B = getBnum(this, gps_time, r, phi, theta)
            % X
            dtheta = 0.1/180*pi;
            V2 = this.getV( gps_time, r, phi, theta-dtheta/2);
            V1 = this.getV( gps_time, r, phi, theta+dtheta/2);
            X = (V2 - V1)/dtheta * 1 / r;
            % Y
            dphi = 0.1/180*pi;
            V2 = this.getV( gps_time, r, phi+dphi/2, theta);
            V1 = this.getV( gps_time, r, phi-dphi/2, theta);
            Y = -(V2 - V1)/dphi * 1 / (r * sin(pi/2 -theta));
            %Z
            dr = 1;
            V2 = this.getV( gps_time, r+dr/2, phi, theta);
            V1 = this.getV( gps_time, r-dr/2, phi, theta);
            Z = (V2 - V1)/dr;
            B = [X;Y;Z];
        end
        
        function B = getB(this, gps_time, r, lon, lat)
            % interpolate H at G at presetn epoch
            y = floor(gps_time/(86400*365.25));
            doy = floor((gps_time - y*(86400*365.25))/86400); %% approximation, but magnetic field varies very slowly
            y = y + 1980;
            re = GPS_SS.ELL_A / 1000;
            if this.cache
                if this.c_year == y || this.c_doy == doy
                    H = this.Hcache;
                    G = this.Gcache;
                else
                    this.cache = false;
                end
            end
            if ~this.cache
                iy = max(min(floor((y-1900)/5)+1,24),1);
                int_time = ((y - this.years(iy))*365.25 + doy-1)/(365.25*5);
                H = this.H(:,:,iy) * (1 - int_time) + int_time *  this.H(:,:,iy+1);
                G = this.G(:,:,iy) * (1 - int_time) + int_time *  this.G(:,:,iy+1);
                this.Hcache = H;
                this.Gcache = G;
                this.c_year = y;
                this.c_doy = doy;
                this.cache = true;
            end
            n = size(H,1);
            % co to colatitude
            colat = pi/2 - lat;
            P = this.interpolateP(cos(colat));
            mphi = repmat((0:n-1)',1,n)*lon;
            cosm = cos(mphi); %some unnecesaary calculations
            sinm = sin(mphi); %some unnecesaary calculations
            arn = repmat((re/r).^(1:n),n,1);
            
            
            
            
            N = repmat((0:n-1),n,1);
            M = N';
            
            % E dV/dphi
            marn = repmat((0:n-1)',1,n) .* arn;
            Ea =  -re / (r * sin(colat)) * sum(sum(marn .* (-G .* sinm + H.* cosm) .* P));
            
            % N dV/dtheta
            dP = this.interpolatedP(cos(colat));
            No =   re / r * sum(sum(arn .* (G .* cosm + H .* sinm) .* dP));
            % U dV/dr
            darn = repmat((1:n),n,1) .* arn;
            Up = - re/r * sum(sum(darn .* (G .* cosm + H .* sinm) .* P));
            
            B = [Ea;No;Up]/ 1e9; %nT to T
            B = local2globalVel2(B, lon,lat);
            %%% rotate B into cartesian
            %             R = [-sin(lon) cos(lon) 0;
            %                 -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
            %                 +cos(lat)*cos(lon) +cos(lat)*sin(lon) sin(lat)];
            %             [B] = R'*B;
        end
        
        function P = interpolateP(this, x)
            n1 = floor((x +1) / this.P_d_length) + 1;
            n2 = min(n1+1, 2/this.P_d_length + 1);
            d = ((x +1) - (n1-1) * this.P_d_length) / this.P_d_length;
            P = this.P(:,:,n1) * (1 - d) + this.P(:,:,n2) * d;
        end
        
        function dP = interpolatedP(this, x)
            n1 = max(0,floor((x + 1 - this.P_d_length) / this.P_d_length)) + 1;
            n2 = min(n1+1, 2/this.P_d_length - 1);
            d = ((x +1) - (n1) * this.P_d_length) / this.P_d_length;
            dP = this.dP(:,:,n1) * (1 - d) + this.dP(:,:,n2) * d;
        end
    end
    
end
