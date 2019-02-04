% INSTRUCTION:
% - Set up goGPS in PPP mode
% - Run goGPS with only PREPRO command
% - Launch the script, it will compute a kinematic solution for the first
%   receiver


% set up the PPP system
ls = LS_Manipulator();
id_sync = rec(1).work.getIdSync;
ls.setUpPPP(rec(1).work,id_sync , core.state.cut_off, true);
time = rec(1).work.time.getSubSet(id_sync);
rate = time.getRate();

% First order tykhonv regualrization of tropo paramater, varainces defined
% in the interface
if core.state.flag_tropo
         ls.setTimeRegularization(ls.PAR_TROPO, (core.state.std_tropo)^2 / 3600 * rate );
end
if core.state.flag_tropo_gradient
        ls.setTimeRegularization(ls.PAR_TROPO_N, (core.state.std_tropo_gradient)^2 / 3600 * rate );
        ls.setTimeRegularization(ls.PAR_TROPO_E, (core.state.std_tropo_gradient)^2 / 3600 * rate );
end

% First order tykhonv regualrization coordinate XYZ, the variance if the
% pseudo-observations refers to the receiver processing rate
varx = 1e-5;
vary = 1e-5;
varz = 1e-5;
ls.setTimeRegularization(ls.PAR_X,varx);
ls.setTimeRegularization(ls.PAR_Y,vary);
ls.setTimeRegularization(ls.PAR_Z,varz);
[x] = ls.solve;

%  x are orrections to the actual values in the receiver
xcoo = x(x(:,2) == ls.PAR_X) + rec(1).work.xyz(1);
ycoo = x(x(:,2) == ls.PAR_Y) + rec(1).work.xyz(2) ;
zcoo = x(x(:,2) == ls.PAR_Z) + rec(1).work.xyz(3);
clock = x(x(:,2) == ls.PAR_REC_CLK) + rec(1).work.getDt;
amb = x(x(:,2) == ls.PAR_AMB) ;
tropo = x(x(:,2) == ls.PAR_TROPO)+ rec(1).work.getZtd;
[gradient_n, gradient_e] = rec(1).work.getGradient;
gradient_e = x(x(:,2) == ls.PAR_TROPO_E) + gradient_e;
gradient_n = x(x(:,2) == ls.PAR_TROPO_N) + gradient_n;