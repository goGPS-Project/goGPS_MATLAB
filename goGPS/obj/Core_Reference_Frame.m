classdef Core_Reference_Frame < handle
    % This class contains properties and methods to manage reference frames
    % and station coordinates
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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
        state
        log
        cc
        
        station_code
        xyz
        vxvyvz
        ref_epoch
    end
     methods (Access = 'private')
        % Creator
        function this = Core_Reference_Frame()
            % Core object creator
            this.state = Global_Configuration.getCurrentSettings();
            this.log = Logger.getInstance();
            this.cc = Global_Configuration.getCurrentSettings().getConstellationCollector;
        end
    end
    
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_reference_frame
            
            if isempty(unique_instance_reference_frame)
                this = Core_Reference_Frame();
                unique_instance_reference_frame = this;
            else
                this = unique_instance_reference_frame;
            end
        end
    end
    
    methods
        function loadCrd(this, fname)
            fid = fopen([fname],'r');
            if fid == -1
                this.log.addWarning(sprintf('      File %s not found', fname));
                return
            end
            this.log.addMessage(sprintf('      Opening file %s for reading', fname));
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            if lim(end,3) < 3
                lim(end,:) = [];
            end
            header_line = find(txt(lim(:,1)) == '#');
            for i = header_line
                line = txt(lim(i,1):lim(i,2));
                idx = strfind(line,'Reference epoch');
                if ~isempty(idx)
                    this.ref_epoch = GPS_Time( sscanf(line((idx+15):end),'%f %f %f %f %f %f')');
                end
            end
            lim(header_line,:) = [];
            cols = 1:4;
            idx_sta_cd = repmat(lim(:,1)-1,1,length(cols)) + repmat(cols,size(lim,1),1);
            this.station_code = reshape(txt(idx_sta_cd),size(lim,1),length(cols));
            cols = 8:96;
            idx_vals = repmat(lim(:,1)-1,1,length(cols)) + repmat(cols,size(lim,1),1);
            vals = reshape(sscanf(txt(idx_vals)','%f'),6,size(lim,1))';
            this.xyz = vals(:,1:3);
            this.vxvyvz = vals(:,4:6);
        end
        function [xyz] = getCoo(this, sta_name, epoch)
            idx_sta = idxCharLines(this.station_code,sta_name);
            dt = epoch - this.ref_epoch;
            xyz = this.xyz(idx_sta,:) + dt/(365.25*86400)*this.vxvyvz(idx_sta,:);
        end
    end
    
end
