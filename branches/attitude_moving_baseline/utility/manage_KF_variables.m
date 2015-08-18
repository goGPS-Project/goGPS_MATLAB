function KF_variables = manage_KF_variables(KF_variables, irec, modality)
      
global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old interval azR elR distR azM elM distM PDOP HDOP VDOP KPDOP KHDOP KVDOP doppler_pred_range1_R doppler_pred_range2_R ratiotest mutest succ_rate fixed_solution n_sys


if strcmp(modality,'archive') % copy from global environment to archive struct
    % copy FK results
    KF_variables(irec).Xhat_t_t = Xhat_t_t;
    KF_variables(irec).X_t1_t = X_t1_t;
    KF_variables(irec).T = T;
    KF_variables(irec).I = I;
    KF_variables(irec).Cee = Cee;
    KF_variables(irec).conf_sat = conf_sat;
    KF_variables(irec).conf_cs = conf_cs;
    KF_variables(irec).pivot = pivot;
    KF_variables(irec).pivot_old = pivot_old;
    KF_variables(irec).interval = interval;
    KF_variables(irec).azR = azR;
    KF_variables(irec).elR = elR;
    KF_variables(irec).distR = distR;
    KF_variables(irec).azM = azM;
    KF_variables(irec).elM = elM;
    KF_variables(irec).distM = distM;
    KF_variables(irec).PDOP = PDOP;
    KF_variables(irec).HDOP = HDOP;
    KF_variables(irec).VDOP = VDOP;
    KF_variables(irec).KPDOP = KPDOP;
    KF_variables(irec).KHDOP = KHDOP;
    KF_variables(irec).KVDOP = KVDOP;
    KF_variables(irec).doppler_pred_range1_R = doppler_pred_range1_R;
    KF_variables(irec).doppler_pred_range2_R = doppler_pred_range2_R;
    KF_variables(irec).ratiotest = ratiotest;
    KF_variables(irec).mutest = mutest;
    KF_variables(irec).succ_rate = succ_rate;
    KF_variables(irec).fixed_solution = fixed_solution;
    KF_variables(irec).n_sys = n_sys;
    
elseif strcmp(modality,'restore') % copy from archive struct to global environment
    Xhat_t_t = KF_variables(irec).Xhat_t_t;
    X_t1_t = KF_variables(irec).X_t1_t;
    T = KF_variables(irec).T;
    I = KF_variables(irec).I;
    Cee = KF_variables(irec).Cee;
    conf_sat = KF_variables(irec).conf_sat;
    conf_cs = KF_variables(irec).conf_cs;
    pivot = KF_variables(irec).pivot;
    pivot_old = KF_variables(irec).pivot_old;
    interval = KF_variables(irec).interval;
    azR = KF_variables(irec).azR;
    elR = KF_variables(irec).elR;
    distR = KF_variables(irec).distR;
    azM = KF_variables(irec).azM;
    elM = KF_variables(irec).elM ;
    distM = KF_variables(irec).distM;
    PDOP = KF_variables(irec).PDOP;
    HDOP = KF_variables(irec).HDOP;
    VDOP = KF_variables(irec).VDOP;
    KPDOP = KF_variables(irec).KPDOP;
    KHDOP = KF_variables(irec).KHDOP;
    KVDOP = KF_variables(irec).KVDOP;
    doppler_pred_range1_R = KF_variables(irec).doppler_pred_range1_R;
    doppler_pred_range2_R = KF_variables(irec).doppler_pred_range2_R;
    ratiotest = KF_variables(irec).ratiotest;
    mutest = KF_variables(irec).mutest;
    succ_rate = KF_variables(irec).succ_rate;
    fixed_solution = KF_variables(irec).fixed_solution;
    n_sys = KF_variables(irec).n_sys;
 
    
else
    printf('\n\n ERROR in calling ''manage_KF_variables''\n\n');
    return
end