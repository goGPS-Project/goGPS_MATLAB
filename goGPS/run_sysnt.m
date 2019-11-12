rec2_synt = GNSS_Station(); rec2_synt.synthetizeWork(rec(2).work.xyz, rec(2).work.time);
rec1_synt = GNSS_Station(); rec1_synt.synthetizeWork(rec(1).work.xyz, rec(1).work.time);
core.rec = [core.rec(1:4) rec1_synt rec2_synt];
rec = core.rec;
rec(5).work.preProcessing
rec(5).work.staticPPP
rec(6).work.preProcessing
rec(6).work.staticPPP
% rec(5).work.sat.cycle_slip_ph_by_ph = rec(5).work.sat.cycle_slip_ph_by_ph | rec(1).work.sat.cycle_slip_ph_by_ph;
% rec(6).work.sat.cycle_slip_ph_by_ph = rec(6).work.sat.cycle_slip_ph_by_ph | rec(2).work.sat.cycle_slip_ph_by_ph;
core.rec(5).marker_name = 'snt1';
core.rec(6).marker_name = 'snt2';
core.exec('NET T5,6 R5 -u -i')