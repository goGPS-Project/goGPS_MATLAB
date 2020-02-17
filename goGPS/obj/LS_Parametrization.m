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
%    |___/                    v 1.0b6
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
        
        sat_x = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        sat_y = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        sat_z = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        
        rec_eb =     [LS_Parametrization.CONST      LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.SING_TRACK];
        rec_ppb =    [LS_Parametrization.STEP_CONST LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.RULE];
        rec_eb_lin = [LS_Parametrization.CONST      LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT LS_Parametrization.RULE];
        rec_ebfr = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_BAND];


        
        sat_eb =   [LS_Parametrization.CONST      LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_TRACK];
        sat_ppb =  [LS_Parametrization.STEP_CONST LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.RULE];
        sat_ebfr = [LS_Parametrization.SPLINE_CUB LS_Parametrization.ALL_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_BAND];
        
        amb  = [LS_Parametrization.STEP_CONST LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.SING_TRACK];

        rec_clk = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        sat_clk = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC  LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        
        rec_clk_pr = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        sat_clk_pr = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC  LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        rec_clk_ph = [LS_Parametrization.EP_WISE LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        sat_clk_ph = [LS_Parametrization.EP_WISE LS_Parametrization.ALL_REC  LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        
        tropo   = [LS_Parametrization.EP_WISE    LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        tropo_n = [LS_Parametrization.EP_WISE    LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        tropo_e = [LS_Parametrization.EP_WISE    LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        tropo_s = [LS_Parametrization.SPLINE_CUB LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
        tropo_v = [LS_Parametrization.CONST      LS_Parametrization.ALL_REC  LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];
        tropo_z = [LS_Parametrization.SPLINE_CUB LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.ALL_FREQ];

              
        iono =    [LS_Parametrization.EP_WISE    LS_Parametrization.SING_REC LS_Parametrization.SING_SAT LS_Parametrization.ALL_FREQ];
             
        ant_mp =  [LS_Parametrization.CONST      LS_Parametrization.SING_REC LS_Parametrization.ALL_SAT  LS_Parametrization.SING_FREQ];
        
        % options to keep track of spline rate, rule based distinction,
        % cycle slips
        rec_x_opt;
        rec_y_opt;
        rec_z_opt;
        rec_eb_opt = struct('rule',1);
        rec_ppb_opt;
        rec_ebfr_opt = struct('spline_rate',3600);

        sat_eb_opt = struct('rule',1);
        sat_ppb_opt;
        sat_ebfr_opt = struct('rule',1,'spline_rate',3600);


        rec_eb_opt_lin = struct('rule',1);

        amb_opt;
        rec_clk_opt;
        tropo_opt;
        tropo_n_opt;
        tropo_e_opt;
        tropo_v_opt;
        sat_clk_opt;
        rec_clk_pr_opt;
        sat_clk_pr_opt;
        rec_clk_ph_opt;
        sat_clk_ph_opt;
        ant_mp_opt;
        iono_opt;
        tropo_s_opt= struct('spline_rate',900);
        tropo_z_opt= struct('spline_rate',1800);

        sat_x_opt = struct('rule',1,'spline_rate',3600);
        sat_y_opt = struct('rule',1,'spline_rate',3600);
        sat_z_opt = struct('rule',1,'spline_rate',3600);
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
        
        function setRate(this, par_class, rate)
               % set the rate for the paramter
            % class p
            %
            % SYNTAX:
            %    [parametriz, option] = setRate(this, par_class, rate)
            switch par_class
                case LS_Manipulator_new.PAR_REC_X
                    this.rec_x_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_Y
                    this.rec_y_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_Z
                    this.rec_z_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_X
                    this.sat_x_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_Y
                    this.sat_y_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_Z
                    this.sat_z_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_EB
                    this.rec_eb_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_EB_LIN
                    this.rec_eb_lin_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_EB
                    this.sat_eb_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_AMB
                    this.amb_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_CLK
                    this.rec_clk_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_CLK
                    this.sat_clk_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_CLK_PR
                    this.rec_clk_pr_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_CLK_PR
                    this.sat_clk_pr_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_CLK_PH
                    this.rec_clk_ph_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_CLK_PH
                    this.sat_clk_ph_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO
                    this.tropo_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO_N
                    this.tropo_n_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO_S
                    this.tropo_s_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO_E
                    this.tropo_e_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO_V
                    this.tropo_v_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_TROPO_Z
                    this.tropo_z_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_IONO
                    this.iono_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_ANT_MP
                    this.ant_mp_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_PPB
                    this.rec_ppb_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_PPB
                    this.sat_ppb_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_REC_EBFR
                    this.rec_ebfr_opt.spline_rate = rate;
                case LS_Manipulator_new.PAR_SAT_EBFR
                    this.sat_ebfr_opt.spline_rate = rate;
            end
        end
        
          function setTimeParametrization(this, par_class, time_parm)
               % set the rate for the paramter
            % class p
            %
            % SYNTAX:
            %    [parametriz, option] = setRate(this, par_class, rate)
            switch par_class
                case LS_Manipulator_new.PAR_REC_X
                    this.rec_x(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_Y
                    this.rec_y(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_Z
                    this.rec_z(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_X
                    this.sat_x(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_Y
                    this.sat_y(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_Z
                    this.sat_z(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_EB
                    this.rec_eb(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_EB_LIN
                    this.rec_eb_lin(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_EB
                    this.sat_eb(1) = time_parm;
                case LS_Manipulator_new.PAR_AMB
                    this.amb(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_CLK
                    this.rec_clk(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_CLK
                    this.sat_clk(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_CLK_PR
                    this.rec_clk_pr(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_CLK_PR
                    this.sat_clk_pr(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_CLK_PH
                    this.rec_clk_ph(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_CLK_PH
                    this.sat_clk_ph(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO
                    this.tropo(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO_N
                    this.tropo_n(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO_S
                    this.tropo_s(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO_E
                    this.tropo_e(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO_V
                    this.tropo_v(1) = time_parm;
                case LS_Manipulator_new.PAR_TROPO_Z
                    this.tropo_z(1) = time_parm;
                case LS_Manipulator_new.PAR_IONO
                    this.iono(1) = time_parm;
                case LS_Manipulator_new.PAR_ANT_MP
                    this.ant_mp(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_PPB
                    this.rec_ppb(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_PPB
                    this.sat_ppb(1) = time_parm;
                case LS_Manipulator_new.PAR_REC_EBFR
                    this.rec_ebfr(1) = time_parm;
                case LS_Manipulator_new.PAR_SAT_EBFR
                    this.sat_ebfr(1) = time_parm;
            end
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
                case LS_Manipulator_new.PAR_SAT_X
                    parametriz = this.sat_x;
                    option = this.sat_x_opt;
                case LS_Manipulator_new.PAR_SAT_Y
                    parametriz = this.sat_y;
                    option = this.sat_y_opt;
                case LS_Manipulator_new.PAR_SAT_Z
                    parametriz = this.sat_z;
                    option = this.sat_z_opt;
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
                case LS_Manipulator_new.PAR_REC_CLK_PR
                    parametriz = this.rec_clk_pr;
                    option = this.rec_clk_pr_opt;
                case LS_Manipulator_new.PAR_SAT_CLK_PR
                    parametriz = this.sat_clk_pr;
                    option = this.sat_clk_pr_opt;
                case LS_Manipulator_new.PAR_REC_CLK_PH
                    parametriz = this.rec_clk_ph;
                    option = this.rec_clk_ph_opt;
                case LS_Manipulator_new.PAR_SAT_CLK_PH
                    parametriz = this.sat_clk_ph;
                    option = this.sat_clk_ph_opt;
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
