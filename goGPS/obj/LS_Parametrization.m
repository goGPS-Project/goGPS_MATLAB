%   CLASS LS_Parametrization
% =========================================================================
%
% DESCRIPTION
%   Class to handle the possible paramterization of the various GNSS
%   unknonw in a LS system
%
% EXAMPLE
%   LSP = LS_Parametrization();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 4 ION
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
classdef LS_Parametrization < handle
    properties (Constant)
        NONE = 0;
        
        
        % time paramterization
        CONST = 1;
        EP_WISE = 2;
        STEP_CONST = 3;
        SPLINE_ZERO = 4;
        SPLINE_LIN = 5;
        SPLINE_CUB = 6;
        RULE_BASED = 7;
        
        % rec paramterization
        SING_REC = 1;
        MULTI_REC = 2;
        ALL_REC = 3;
        
        % sat paramterization
        SING_SAT = 1;
        MULTI_SAT = 2;
        ALL_SAT = 3;
        
        % track paramterization paramterization
        SING_TRACK = 1;
        SING_FREQ = 2;
        FREQ_CONST = 6;
        SING_BAND = 7;
        SING_FREQ_BIN = 4;
        ALL_FREQ = 3;
        RULE = 5; % warnign only sequential & allowed
        
        % spatial paramterization 
        %TBD
    end
    
    properties
        % paramterization [(time paramterization) (rec paramterization) (sat paramterization) (tracking paramterization)]
        rec_x = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        rec_y = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        rec_z = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        
        sat_x = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        sat_y = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        sat_z = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        
        rec_eb = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.SING_TRACK];
        rec_ppb = [LS_Parametrization.STEP_CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.RULE];
        rec_eb_lin = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.RULE];

        
        sat_eb = [LS_Parametrization.CONST LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_TRACK];
        sat_ppb = [LS_Parametrization.STEP_CONST LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.RULE];
        sat_ebfr = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_BAND];
        
        amb  = [LS_Parametrization.STEP_CONST LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_TRACK];

        rec_clk = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        sat_clk = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
      
        tropo = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        tropo_n = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        tropo_e = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        tropo_s = [LS_Parametrization.SPLINE_CUB LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        tropo_v = [LS_Parametrization.CONST LS_Parametrization.ALL_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];
        tropo_z = [LS_Parametrization.SPLINE_CUB LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.ALL_FREQ];

              
        iono =  [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
             
        ant_mp = [LS_Parametrization.CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.SING_FREQ];
        
        % options to keep track of spline rate, rule based distinction,
        % cycle slips
        rec_x_opt;
        rec_y_opt;
        rec_z_opt;
        rec_eb_opt = struct('rule',1);
        rec_ppb_opt;
        sat_eb_opt = struct('rule',1);
        sat_ppb_opt;
        sat_ebfr_opt = struct('rule',1,'spline_rate',900);

        rec_eb_opt_lin = struct('rule',1);

        amb_opt;
        rec_clk_opt;
        tropo_opt;
        tropo_n_opt;
        tropo_e_opt;
        tropo_v_opt;
        sat_clk_opt;
        ant_mp_opt;
        iono_opt;
        tropo_s_opt= struct('spline_rate',900);
        tropo_z_opt= struct('spline_rate',1800);

        sat_x_opt;
        sat_y_opt;
        sat_z_opt;
    end
    
    methods
        function [this] = LS_Parametrization()
            this = this@handle();
          %  this.rec_eb_opt.rule = {['PSRANGE:' num2str(LS_Parametrization.SING_TRACK)],['PHASE&NOT*GLONASS:' num2str(LS_Parametrization.SING_TRACK)],['PHASE&GLONASS:' num2str(LS_Parametrization.SING_FREQ)]};
            %this.rec_eb_opt.rule = {['PSRANGE:' num2str(LS_Parametrization.SING_TRACK)],['PHASE:' num2str(LS_Parametrization.FREQ_CONST)]};
            %this.sat_eb_opt.rule = {['PSRANGE:' num2str(LS_Parametrization.SING_TRACK)],['PHASE:' num2str(LS_Parametrization.NONE)]};
            this.sat_ppb_opt.rule = {['PHASE:' num2str(LS_Parametrization.ALL_FREQ)],['PSRANGE:' num2str(LS_Parametrization.NONE)]};
            this.rec_ppb_opt.rule = {['PHASE:' num2str(LS_Parametrization.ALL_FREQ)],['PSRANGE:' num2str(LS_Parametrization.NONE)]};
            this.rec_eb_opt_lin.rule = {['PHASE&GLONASS:' num2str(LS_Parametrization.SING_FREQ)]};
        end
       
        function [parametriz, option] = getParametrization(this, par_class, obs_code)
            % get the parametrization and the options for the paramter
            % class p
            %
            % SYNTAX:
            %    [parametriz, option] = getParametrization(this, par_class)
            switch par_class
                case LS_Manipulator_new.PAR_REC_X
                    parametriz = this.rec_x;
                    option = this.rec_x_opt;
                case LS_Manipulator_new.PAR_REC_Y
                    parametriz = this.rec_y;
                    option = this.rec_y_opt;
                case LS_Manipulator_new.PAR_REC_Z
                    parametriz = this.rec_z;
                    option = this.rec_z_opt;
                case LS_Manipulator_new.PAR_REC_EB
                    parametriz = this.rec_eb;
                    option = this.rec_eb_opt;
                case LS_Manipulator_new.PAR_REC_EB_LIN
                    parametriz = this.rec_eb_lin;
                    option = this.rec_eb_opt_lin;
                case LS_Manipulator_new.PAR_SAT_EB
                    parametriz = this.sat_eb;
                    option = this.sat_eb_opt;
                case LS_Manipulator_new.PAR_AMB
                    parametriz = this.amb;
                    option = this.amb_opt;
                case LS_Manipulator_new.PAR_REC_CLK
                    parametriz = this.rec_clk;
                    option = this.rec_clk_opt;
                case LS_Manipulator_new.PAR_SAT_CLK
                    parametriz = this.sat_clk;
                    option = this.sat_clk_opt;
                case LS_Manipulator_new.PAR_TROPO
                    parametriz = this.tropo;
                    option = this.tropo_opt;
                case LS_Manipulator_new.PAR_TROPO_N
                    parametriz = this.tropo_n;
                    option = this.tropo_n_opt;
                case LS_Manipulator_new.PAR_TROPO_S
                    parametriz = this.tropo_s;
                    option = this.tropo_s_opt;
                case LS_Manipulator_new.PAR_TROPO_E
                    parametriz = this.tropo_e;
                    option = this.tropo_e_opt;
                case LS_Manipulator_new.PAR_TROPO_V
                    parametriz = this.tropo_v;
                    option = this.tropo_v_opt;
                case LS_Manipulator_new.PAR_TROPO_Z
                    parametriz = this.tropo_z;
                    option = this.tropo_z_opt;
                case LS_Manipulator_new.PAR_IONO
                    parametriz = this.iono;
                    option = this.iono_opt;
                case LS_Manipulator_new.PAR_ANT_MP
                    parametriz = this.ant_mp;
                    option = this.ant_mp_opt;
                case LS_Manipulator_new.PAR_REC_PPB
                    parametriz = this.rec_ppb;
                    option = this.rec_ppb_opt;
                case LS_Manipulator_new.PAR_SAT_PPB
                    parametriz = this.sat_ppb;
                    option = this.sat_ppb_opt;
                case LS_Manipulator_new.PAR_REC_EBFR
                    parametriz = this.rec_ebfr;
                    option = this.rec_ebfr_opt;
                case LS_Manipulator_new.PAR_SAT_EBFR
                    parametriz = this.sat_ebfr;
                    option = this.sat_ebfr_opt;
            end
            
            if nargin > 2 && parametriz(4) == this.RULE % if an obseravtion code is specified give the paramterization for thata specific obsservation code
                rule = option.rule;
                for r = 1 : length(rule)
                    parts = strsplit(rule{r},':');
                    cond = parts{1};
                    optt = parts{2};
                    conds = strsplit(cond,'&');
                    cond_passed = true;
                    for c = conds
                        if strcmpi('PSRANGE',c)
                            cond_passed = cond_passed && obs_code(2) == 'C';
                        elseif strcmpi('PHASE',c)
                            cond_passed = cond_passed && obs_code(2) == 'L';
                        elseif strcmpi('GLONASS',c)
                            cond_passed = cond_passed && obs_code(1) == 'R';
                        elseif strcmpi('NOT*GLONASS',c)
                            cond_passed = cond_passed && obs_code(1) ~= 'R';
                        else
                            cond_passed = false;
                        end
                    end
                    if cond_passed
                        parametriz(4) = str2num(optt);
                    end
                end
            end
        end
    end
end
