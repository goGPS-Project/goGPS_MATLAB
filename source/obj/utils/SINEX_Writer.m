%   CLASS SINEX_Writer
% =========================================================================
%
% DESCRIPTION
%   Collection of methods to produce a sinex file
%
% EXAMPLE
%   
%
% SEE ALSO
%  
% FOR A LIST OF CONSTANTs and METHODS use doc SINEX_Writer

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, ...
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
classdef SINEX_Writer < handle
    properties (Constant)
        TROPO_VERSION_HEADER = '%=TRO 1.00 GRD ${DATEPRODUCED} GRD ${STARTOFDATA} ${ENDOFDATA} P ${MARKERNAME}'
        
        
        DESCRIPTION = 'GReD Geomatic Research and Development, Lomazzo Italy';
        OUTPUT_TROPO = 'Total Troposphere Zenith Path Delay Product';
        CONTACT = 'support@gogps-project.org';
        SOFTWARE = 'goGPS Software';
        INPUT = 'IGS final GPS orbit and clock solutions, site RINEX files.';
        ACKNOWLEDGMENTS = 'International GNSS Service (IGS)';
        
        TROPO_HEADER = '*SITE EPOCH_______ TROTOT STDEV  TGNTOT  STDEV  TGETOT  STDEV';
        SINEX_MAPPING_FLAGS = {'WET GMF','WET GRID VMF1 2.5x2.5',' WET NIELL','WET GRID VMF3 1x1','WET GRID VMF3 5x5'};
        SUPPORTED_PARAMETERS = {'TROTOT','TGNTOT','TGETOT','TROWET','PWV','PRESS','TEMDRY','HUMREL'};
    end
    properties
        state;
        fid;
        fname;
        tropo_fields;
        std_fields;
    end
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
    end
    methods
        function this = SINEX_Writer(fname)
            this.state = Core.getState();
            this.fname = fname;
            this.fid = fopen(fname,'Wb');
        end
        function close(this)
            fclose(this.fid);
        end
        function writeTroSinexHeader(this, date_st, date_en, mark_name)
            fid = this.fid;
            gps_now = GPS_Time.now();
            exp_date = gps_now.toSinexStrDate();
            st_date = date_st.toSinexStrDate();
            en_date = date_en.toSinexStrDate();
            head_line = strrep(strrep(strrep(strrep(this.TROPO_VERSION_HEADER,'${DATEPRODUCED}',exp_date),'${ENDOFDATA}',en_date),'${STARTOFDATA}',st_date),'${MARKERNAME}',mark_name);
            fprintf(fid, ['%' head_line '\n']);
        end
        function writeFileReference(this)
            fid = this.fid;
            fprintf(fid,'+FILE/REFERENCE\n');
            fprintf(fid,' %-21s%s\n','DESCRIPTION',this.DESCRIPTION);
            fprintf(fid,' %-21s%s\n','OUTPUT',this.OUTPUT_TROPO);
            fprintf(fid,' %-21s%s\n','CONTACT',this.CONTACT);
            fprintf(fid,' %-21s%s\n','SOFTWARE',this.SOFTWARE);
            if isunix()
                try
                    [status,cmdout] = system('hostname');
                    if status == 0
                        fprintf(fid,' %-21s%s','HARDWARE',cmdout);
                    end
                catch
                end
            end
            fprintf(fid,' %-21s%s\n','INPUT',this.INPUT);
            fprintf(fid,'-FILE/REFERENCE\n\n');
        end
        function writeAcknoledgments(this)
            fid = this.fid;
            fprintf(fid,'+INPUT/ACKNOWLEDGMENTS\n');
            fprintf(fid,' %s\n',this.ACKNOWLEDGMENTS);
            fprintf(fid,'-INPUT/ACKNOWLEDGMENTS\n\n');
        end
        function writeSTACoo(this, mrk_name, x, y, z, ref_sys, rmrk)
            fid = this.fid;
            fprintf(fid,'+TROP/STA_COORDINATES\n');
            fprintf(fid,'*SITE PT SOLN T __STA_X_____ __STA_Y_____ __STA_Z_____ SYSTEM REMRK\n');
            
            for i = size(x,1)
                fprintf(fid,' %-4s%-3s%-5s P %-13.3f%-13.3f%-13.3f%-7s%-4s \n',mrk_name(i,:),'A','1',x(i),y(i),z(i),ref_sys(i,:),rmrk(i,:));
            end
            fprintf(fid,'-TROP/STA_COORDINATES\n\n');
        end
        function writeTropoDescription(this, cut_off, sampl_obs, smpl_tropo, map_fun,sol_fields, std_fields)
            %INPUT : cut_off (num)
            %        sampl_obs (num)
            %        sampl_tropo (num)
            %        map_fun (string)
            %        sol_fields (cell)
            %        std_fields (bool) size of sol_fields
            fid = this.fid;
            fprintf(fid,'+TROP/DESCRIPTION\n');
            fprintf(fid,'*_________KEYWORD_____________ __VALUE(S)______________________\n');
            fprintf(fid,' %-32s%-20d\n','SAMPLING INTERVAL',sampl_obs);
            fprintf(fid,' %-32s%-20d\n','SAMPLING TROP',smpl_tropo);
            fprintf(fid,' %-32s%-20d\n','ELEVATION CUTOFF ANGLE',cut_off);
            fprintf(fid,' %-32s%s\n','TROP MAPPING FUNCTION',map_fun);
            this.tropo_fields = sol_fields;
            this.std_fields = std_fields;
            sol_field_string = '';
            for i = 1 : length(this.tropo_fields)
                sol_field_string = [sol_field_string ' ' this.tropo_fields{i}];
                if std_fields(i)
                    sol_field_string = [sol_field_string ' STDDEV'];
                end
                
            end
            fprintf(fid,' %-32s%s\n','SOLUTION_FIELDS_1',sol_field_string);
            fprintf(fid,'-TROP/DESCRIPTION\n\n');
        end
        function writeTropoSolutionSt(this)
            fid = this.fid;
            fprintf(fid,'+TROP/SOLUTION\n');
            sol_field_string = '';
            for i = 1 : length(this.tropo_fields)
                if strfind(this.tropo_fields{i},'STDEV')
                    sol_field_string = [sol_field_string sprintf('%-5s',this.tropo_fields{i})];
                else
                    sol_field_string = [sol_field_string sprintf('%-8s',this.tropo_fields{i})];
                end
            end
            fprintf(fid,['*SITE EPOCH_______ ' sol_field_string '\n']);
        end
        function writeTropoSolutionStation(this, mrk_mame, gps_time, vals, stds, vals_flag)
            % IMPRTANT: vals_flags has to be consistent with the one passed
            % in this.writeTropoDescription()
            fid = this.fid;
            n_ep = gps_time.length;
            yy_doy_sod = gps_time.toSinexStrDate();
            mrk_name = repmat([' ' mrk_mame ' '], n_ep, 1);
            vals_mat = [];
            vals_flag = this.SUPPORTED_PARAMETERS(vals_flag>0);
            for i = 1 : size(vals,2)
                if strfind(vals_flag{i},'TROTOT')
                    vals_mat = [ vals_mat reshape(sprintf('%7.1f',vals(:,i)),7,n_ep)'];
                    if this.std_fields(i)
                        vals_mat = [ vals_mat reshape(sprintf('%5.1f',stds(:,i)),5,n_ep)'];
                    end
                elseif strfind(vals_flag{i},'TG')                    
                    vals_mat = [ vals_mat reshape(sprintf('%8.3f', max(-999,min(999,vals(:,i)))),8,n_ep)'];
                    if this.std_fields(i)
                        vals_mat = [ vals_mat reshape(sprintf('%5.3f',stds(:,i)),5,n_ep)'];
                    end
                else
                    vals_mat = [ vals_mat reshape(sprintf('%7.1f',vals(:,i)),7,n_ep)'];
                    if this.std_fields(i)
                        vals_mat = [ vals_mat reshape(sprintf('%5.1f',stds(:,i)),5,n_ep)'];
                    end
                end
            end
            vals_mat = [ mrk_name yy_doy_sod vals_mat repmat('\n',n_ep,1)];
            fprintf(fid,vals_mat');
            
        end
        function writeTropoSolutionEnd(this)
            fid = this.fid;
            fprintf(fid,'-TROP/SOLUTION\n');
        end
        function writeTroSinexEnd(this)
            fid = this.fid;
            fprintf(fid,'%%=ENDTRO');
        end
        
    end
end
