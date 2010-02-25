/*****************************************************************\
*                            .88888.   888888ba  .d88888b         *
*                           d8'   `88  88    `8b 88.    "'        *
*         .d8888b. .d8888b. 88        a88aaaa8P' `Y88888b.        *
*         88'  `88 88'  `88 88   YP88  88              `8b        *
*         88.  .88 88.  .88 Y8.   .88  88        d8'   .8P        *
*         `8888P88 `88888P'  `88888'   dP         Y88888P         *
*              .88                                                *
*          d8888P                                                 *
\*****************************************************************/

1. Introduction
   ============
   goGPS is a software package designed to perform GPS navigation
   with low cost receivers, either in post-processing or real-time.
   It is developed in MATLAB and it is aimed at providing a tool
   useful for studying GPS navigation, implementing and testing new
   algorithms and interacting in general with GPS-related aspects.

2. Requirements
   ============
   goGPS has been developed and tested in MATLAB 7.* environments,
   running on Windows machines. The following elements are needed
   in order to use goGPS:

   - a computer with Windows operating system
   - a MATLAB 7.* installation

   For post-processing tasks:
   - RINEX observation file for the roving receiver
   - RINEX observation file for the master station
   - RINEX navigation file
     (RINEX files must have epochs in common)
     (goGPS binary data saved during a real-time session
     can be used instead of RINEX files)

   For real-time tasks:
   - 'Instrument Control Toolbox' installed on MATLAB
   - GPS receiver providing raw data on a COM port (currently only
     u-blox binary protocol is supported)
   - GPS permanent station(s) broadcasting raw data in RTCM 3.x
     format through NTRIP protocol (at least '1002' or '1004'
     messages)

3. Basic settings
   ========
   goGPS functioning mode is set in the goGPS.m file.
   Here it is possible to choose between various POST-PROCESSING
   and REAL-TIME modes, free and constrained mode and to set some
   flags. Comments on the code should be self-explanatory.

   User-defined settings are included in the file
   "global_settings.m".

   filerootIN : input path and filename prefix for loading
                goGPS-saved binary data
   filerootOUT : output path and filename prefix for saving goGPS
                 binary data
   (e.g. '../data/test' to load or save 'test_*' files in '../data'
   folder)
   ----------------------------------------------------------------
   filename_R_obs : input filename for RINEX rover observations
   filename_R_nav : input filename for RINEX rover navigation file
                    (if it is not available, insert the master one)
   filename_M_obs : input filename for RINEX master observations
   filename_M_nav : input filename for RINEX master navigation file
   ----------------------------------------------------------------
   filename_ref : name of a MAT file containing a reference path
                  described by the following two variables:
                  - 'ref_path' containing vertices in geocentric
                    cartesian coordinates
                  - 'mat_path' containing an adjacency matrix to
                    describe connections between nodes
   ----------------------------------------------------------------
   XM, YM, ZM : base station geocentric cartesian coordinates
   ----------------------------------------------------------------

4.1 HOW-TO: Post-processing
    =======================

4.2 HOW-TO: Real-time
    =================

5. File list
   =========
    amb_estimate_approx.m
    amb_estimate_observ.m
    bancroft.m
    base64encode.m
    cart2geod.m
    cart2plan.m
    check_parity.m
    check_t.m
    code_double_diff.m
    code_phase_double_diff.m
    code_SA.m
    cofactor_matrix.m
    crc24q.m
    cycle_slip_kalman.m
    cycle_slip_LS_N.m
    decode_3.m
    decode_18.m
    decode_19.m
    decode_1001.m
    decode_1002.m
    decode_1003.m
    decode_1004.m
    decode_1005.m
    decode_1006.m
    decode_1007.m
    decode_1008.m
    decode_1009.m
    decode_1010.m
    decode_1011.m
    decode_1012.m
    decode_1013.m
    decode_1014.m
    decode_1015.m
    decode_1016.m
    decode_1017.m
    decode_1019.m
    decode_1020.m
    decode_rtcm2.m
    decode_rtcm3.m
    decode_RXM_EPH.m
    decode_RXM_RAW.m
    decode_ublox.m
    dtm_bilin_interp.m
    dtm_bulk_load.m
    dtm_hdr_load.m
    e_r_corr.m
    err_iono.m
    err_tropo.m
    find_eph.m
    geoc2geod.m
    geod2cart.m
    geod2plan.m
    global2localCov.m
    global_init.m
    global_settings.m
    goGPS.m
    goGPS_gui.fig
    goGPS_gui.m
    goGPS_master_monitor.m
    goGPS_realtime.m
    goGPS_realtime_monitor.m
    goGPS_ublox_monitor.m
    gps_time.m
    input_bancroft.m
    input_kalman.m
    input_kalman_vinc.m
    julday.m
    kalman_goGPS_cod_init.m
    kalman_goGPS_cod_loop.m
    kalman_goGPS_init.m
    kalman_goGPS_loop.m
    kalman_goGPS_SA_cod_init.m
    kalman_goGPS_SA_cod_loop.m
    kalman_goGPS_vinc_init.m
    kalman_goGPS_vinc_loop.m
    KML_link_write.m
    KML_update.m
    KML_write.m
    KML_write_cov.m
    LICENSE.txt
    load_goGPSinput.m
    load_goGPSinput_old.m
    load_goGPSoutput.m
    load_goGPStime.m
    load_RINEX.m
    load_stream.m
    lorentz.m
    MQ_goGPS_loop.m
    MQ_goGPS_SA_loop.m
    NMEA_string_generator.m
    NTRIP_string_generator.m
    obs_type_find.m
    obs_type_list.m
    README.txt
    ref_2d_projection.m
    ref_3d_projection.m
    RINEX_get_epoch.m
    RINEX_get_nav.m
    RINEX_get_nav_GLO.m
    RINEX_get_obs.m
    RINEX_jump_hdr.m
    rt_find_eph.m
    rtplot_amb.m
    rtplot_googleearth.m
    rtplot_googleearth_cov.m
    rtplot_matlab.m
    rtplot_matlab_cov.m
    rtplot_skyplot.m
    rtplot_snr.m
    rttext_sat.m
    sat_corr.m
    sat_pos.m
    topocent.m
    twos_complement.m
    ublox2RINEX.m
    ublox_CFG_CFG.m
    ublox_CFG_MSG.m
    ublox_COM_find.m
    ublox_poll_message.m
    ublox_UBX_codes.m
