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

   v0.2.0beta

1. Introduction
   ============
   goGPS is a software package designed to perform GPS navigation
   with low cost receivers, either in post-processing or real-time.
   It is developed in MATLAB and it is aimed at providing a tool
   useful for studying GPS navigation, implementing and testing new
   algorithms and interacting in general with GPS-related aspects.

2. Requirements
   ============
   goGPS has been developed and tested in MATLAB 7.6+ environments,
   on both Windows and UNIX. The following elements are needed
   in order to use goGPS:

   - a computer with Windows or a UNIX-based operating system
   - a MATLAB 7.6+ installation

   For post-processing tasks:
   - RINEX observation file for the roving receiver
   - RINEX observation file for the master station
   - RINEX navigation file
     (RINEX files must have epochs in common)
     (goGPS binary data saved during a real-time session
     can be used instead of RINEX files)

   For real-time tasks:
   - 'Instrument Control Toolbox' installed on MATLAB
   - GPS receiver providing raw data on a COM port (currently
     u-blox UBX, Fastrax IT03 and SkyTraq binary protocols are
     supported) with their own drivers installed
   - GPS permanent station(s) broadcasting raw data in RTCM 3.x
     format through NTRIP protocol (at least '1002' or '1004'
     messages)

   NOTE 1: plotting on Google Earth requires it to be installed;
           if plotting error ellipses on Google Earth produces odd
           output on Windows, please switch Google Earth rendering
           engine to DirectX.
		   
   NOTE 2: if you experience problems using the graphical user
           interface under UNIX systems (Linux or Mac), please try
		   disabling the UNIX GUI (comment lines 70 and 75-80 in 
		   goGPS.m), leaving only the WINDOWS GUI enabled.
		   GUIs can be modified using the MATLAB GUIDE tool.

3. Directory tree:
   ============

goGPS
|   global_init.m
|   global_settings.m
|   goGPS.m
|   goGPS_master_monitor.m
|   goGPS_realtime.m
|   goGPS_realtime_monitor.m
|   goGPS_rover_monitor.m
|   LICENSE.txt
|   param_fastrax.m
|   param_skytraq.m
|   param_ublox.m
|   README.txt
|   update_settings.m
|   update_settings_020beta.m
|   
+---dtm
|       dtm_bulk_load.m
|       dtm_hdr_load.m
|       grid_bilin_interp.m
|       
+---gui
|   +---unix
|   |       gui_about_unix.fig
|   |       gui_about_unix.m
|   |       gui_decode_stream_unix.fig
|   |       gui_decode_stream_unix.m
|   |       gui_goGPS_unix.fig
|   |       gui_goGPS_unix.m
|   |       gui_GPS_week_unix.fig
|   |       gui_GPS_week_unix.m
|   |       gui_merge_goGPSbin_unix.fig
|   |       gui_merge_goGPSbin_unix.m
|   |       gui_polyline_simplification_unix.fig
|   |       gui_polyline_simplification_unix.m
|   |       gui_RINEX2goGPSbin_unix.fig
|   |       gui_RINEX2goGPSbin_unix.m
|   |       
|   \---win
|           gui_about.fig
|           gui_about.m
|           gui_decode_stream.fig
|           gui_decode_stream.m
|           gui_goGPS.fig
|           gui_goGPS.m
|           gui_GPS_week.fig
|           gui_GPS_week.m
|           gui_merge_goGPSbin.fig
|           gui_merge_goGPSbin.m
|           gui_polyline_simplification.fig
|           gui_polyline_simplification.m
|           gui_RINEX2goGPSbin.fig
|           gui_RINEX2goGPSbin.m
|           
+---input_output
|   |   load_stream.m
|   |   NTRIP_string_generator.m
|   |   
|   +---FTX_fastrax
|   |       checksumFTX.m
|   |       configure_fastrax.m
|   |       decode_fastrax_it03.m
|   |       decode_FTX_EPH.m
|   |       decode_FTX_PSEUDO.m
|   |       decode_FTX_TRACK.m
|   |       fastrax_check_ACK.m
|   |       FTX_TypeConv.m
|   |       italk_enable_raw.m
|   |       italk_PPS_MEAS_MS.m
|   |       italk_PPS_SYNC_TRACK.m
|   |       italk_reset_default.m
|   |       italk_TRACK_MEAS_INTERVAL.m
|   |       
|   +---goGPSbin
|   |       goGPSbinMerge.m
|   |       load_goGPSinput.m
|   |       load_goGPSoutput.m
|   |       load_observ.m
|   |       streams2goGPSbin.m
|   |       
|   +---KML
|   |       KML_link_write.m
|   |       KML_write.m
|   |       KML_write_cov.m
|   |       
|   +---NMEA
|   |       NMEA_checksum.m
|   |       NMEA_GGA_gen.m
|   |       NMEA_GSA_gen.m
|   |       NMEA_GSV_gen.m
|   |       NMEA_PGGPK_gen.m
|   |       NMEA_RMC_gen.m
|   |       
|   +---RINEX
|   |       load_RINEX.m
|   |       load_RINEX_SA.m
|   |       obs_type_find.m
|   |       RINEX2goGPSbin.m
|   |       RINEX_get_epoch.m
|   |       RINEX_get_nav.m
|   |       RINEX_get_nav_GLO.m
|   |       RINEX_get_obs.m
|   |       RINEX_parse_hdr.m
|   |       streamM2RINEX.m
|   |       streamR2RINEX.m
|   |       
|   +---RTCM
|   |       decode_1001.m
|   |       decode_1002.m
|   |       decode_1003.m
|   |       decode_1004.m
|   |       decode_1005.m
|   |       decode_1006.m
|   |       decode_1007.m
|   |       decode_1008.m
|   |       decode_1009.m
|   |       decode_1010.m
|   |       decode_1011.m
|   |       decode_1012.m
|   |       decode_1013.m
|   |       decode_1014.m
|   |       decode_1015.m
|   |       decode_1016.m
|   |       decode_1017.m
|   |       decode_1019.m
|   |       decode_1020.m
|   |       decode_1029.m
|   |       decode_18.m
|   |       decode_19.m
|   |       decode_3.m
|   |       decode_rtcm2.m
|   |       decode_rtcm3.m
|   |       
|   +---STQ_skytraq
|   |       configure_skytraq.m
|   |       decode_skytraq.m
|   |       decode_skytraq_GPS_EPH.m
|   |       decode_skytraq_MEAS_TIME.m
|   |       decode_skytraq_RAW_MEAS.m
|   |       skytraq_binary_output_rate.m
|   |       skytraq_check_ACK.m
|   |       skytraq_message_format.m
|   |       skytraq_poll_message.m
|   |       
|   \---UBX_ublox
|           configure_ublox.m
|           decode_AID_EPH.m
|           decode_AID_HUI.m
|           decode_RXM_EPH.m
|           decode_RXM_RAW.m
|           decode_RXM_SFRB.m
|           decode_ublox.m
|           ublox_CFG_CFG.m
|           ublox_CFG_MSG.m
|           ublox_CFG_RATE.m
|           ublox_check_ACK.m
|           ublox_COM_find.m
|           ublox_poll_message.m
|           ublox_UBX_codes.m
|           
+---plot
|       rtplot_amb.m
|       rtplot_googleearth.m
|       rtplot_googleearth_cov.m
|       rtplot_matlab.m
|       rtplot_matlab_cov.m
|       rtplot_matlab_cov_stopGOstop.m
|       rtplot_matlab_stopGOstop.m
|       rtplot_skyplot.m
|       rtplot_snr.m
|       rttext_sat.m
|       
+---polyline
|       LSinterp.m
|       LSsolver.m
|       polyline.m
|       polyline_arcsClustering.m
|       polyline_leastSquaresFit.m
|       polyline_nodesDetection.m
|       
+---positioning
|   |   ambiguity_init.m
|   |   ambiguity_init_SA.m
|   |   amb_estimate_observ.m
|   |   amb_estimate_observ_SA.m
|   |   clock_error.m
|   |   cofactor_matrix.m
|   |   cofactor_matrix_SA.m
|   |   cycle_slip_detection.m
|   |   cycle_slip_detection_SA.m
|   |   doppler_shift_approx.m
|   |   ecc_anomaly.m
|   |   err_iono.m
|   |   err_tropo.m
|   |   e_r_corr.m
|   |   find_eph.m
|   |   rt_find_eph.m
|   |   sat_corr.m
|   |   sat_pos.m
|   |   topocent.m
|   |   
|   +---bancroft
|   |       bancroft.m
|   |       input_bancroft.m
|   |       
|   +---kalman
|   |       input_kalman.m
|   |       input_kalman_SA.m
|   |       input_kalman_vinc.m
|   |       kalman_goGPS_cod_init.m
|   |       kalman_goGPS_cod_loop.m
|   |       kalman_goGPS_init.m
|   |       kalman_goGPS_init_model.m
|   |       kalman_goGPS_loop.m
|   |       kalman_goGPS_loop_model.m
|   |       kalman_goGPS_SA_cod_init.m
|   |       kalman_goGPS_SA_cod_loop.m
|   |       kalman_goGPS_SA_init.m
|   |       kalman_goGPS_SA_loop.m
|   |       kalman_goGPS_vinc_init.m
|   |       kalman_goGPS_vinc_loop.m
|   |       
|   \---least_squares
|           code_double_diff.m
|           code_double_diff.m~
|           code_phase_double_diff.m
|           code_phase_SA.m
|           code_SA.m
|           LS_goGPS_cod_loop.m
|           LS_goGPS_SA_cod_loop.m
|           LS_goGPS_SA_loop.m
|           
\---utility
    |   base64encode.m
    |   crc24q.m
    |   fbin2dec.m
    |   lorentz.m
    |   ref_2d_projection.m
    |   ref_3d_projection.m
    |   twos_complement.m
    |   
    +---geo
    |       cart2geod.m
    |       cart2plan.m
    |       geod2cart.m
    |       geod2plan.m
    |       global2localCov.m
    |       global2localPos.m
    |       local2globalCov.m
    |       local2globalPos.m
    |       utm2deg.m
    |       
    \---gps
            check_t.m
            decode_subframe_1.m
            decode_subframe_2.m
            decode_subframe_3.m
            decode_subframe_4.m
            gps_parity.m
            gps_time.m
            julday.m