classdef Earth_Magnetic_Field < handle
    % This class containing the spherical harmonics model for earth
    % magnetic field
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
    %  Written by: Giulio Tagliaferro
    %  Contributors:     ...
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
        H;%gauss coefficients
        G;%gauss coefficients
    end
    properties (Access = private)
        state;
    end
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core_sky__
            
            if isempty(unique_instance_core_sky__)
                this = Earth_Magnetic_Field();
                unique_instance_core_sky__ = this;
            else
                this = unique_instance_core_sky__;
            end
        end
    end
    methods
        function this = Earth_Magnetic_Field()
            this.state = Go_State.getCurrentSettings();
            this.importIGRFModel(this.state.getIGRFFile());
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
                    H = zeros(n_ep,14,14);
                    G = zeros(n_ep,14,14);
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
                    
                    l_idx = sub2ind(size(H),[1:n_ep]',repmat(n+1,n_ep,1),repmat(m+1,n_ep,1))
                	if tline(1) == 'g'
                        G(l_idx) = vals;
                    else
                        H(l_idx) = vals;
                    end
                end
                tline = fgetl(fid);
            end
            
            fclose(fid);
            this.H = H;
            this.G = G;

        end
        function B = getB(this,gps_time,r,phi,theta)
        end
        
    end
    
end