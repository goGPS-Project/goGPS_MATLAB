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
     u-blox ANTARIS binary protocol is supported)
   - GPS permanent station(s) broadcasting raw data in RTCM 3.x
     format through NTRIP protocol (at least '1002' or '1004'
     messages)

3. File list
   =========

   amb_estimate_approx.m
   amb_estimate_approx_SA.m
   amb_estimate_observ.m
   amb_estimate_observ_SA.m
   bancroft.m
   base64encode.m
   cart2geod.m
   cart2plan.m
   check_t.m
   code_SA.m
   code_double_diff.m
   code_phase_SA.m
   code_phase_double_diff.m
   cofactor_matrix.m
   cofactor_matrix_SA.m
   crc24q.m
   cycle_slip_LS_N.m
   cycle_slip_LS_N_SA.m
   cycle_slip_kalman.m
   cycle_slip_kalman_SA.m
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
   decode_18.m
   decode_19.m
   decode_3.m
   decode_AID_EPH.m
   decode_AID_HUI.m
   decode_RXM_EPH.m
   decode_RXM_RAW.m
   decode_RXM_SFRB.m
   decode_rtcm2.m
   decode_rtcm3.m
   decode_subframe_1.m
   decode_subframe_2.m
   decode_subframe_3.m
   decode_subframe_4.m
   decode_ublox.m
   dtm_bulk_load.m
   dtm_hdr_load.m
   e_r_corr.m
   ecc_anomaly.m
   err_iono.m
   err_tropo.m
   fbin2dec.m
   find_eph.m
   geoc2geod.m
   geod2cart.m
   geod2plan.m
   global2localCov.m
   global_init.m
   global_settings.m
   goGPS.m
   goGPS_master_monitor.m
   goGPS_realtime.m
   goGPS_realtime_monitor.m
   goGPS_ublox_monitor.m
   goGPSbinMerge.m
   gps_parity.m
   gps_time.m
   grid_bilin_interp.m
   gui_GPS_week.fig
   gui_GPS_week.m
   gui_GPS_week_unix.fig
   gui_GPS_week_unix.m
   gui_RINEX2goGPSbin.fig
   gui_RINEX2goGPSbin.m
   gui_RINEX2goGPSbin_unix.fig
   gui_RINEX2goGPSbin_unix.m
   gui_decode_stream.fig
   gui_decode_stream.m
   gui_decode_stream_unix.fig
   gui_decode_stream_unix.m
   gui_goGPS.fig
   gui_goGPS.m
   gui_goGPS_unix.fig
   gui_goGPS_unix.m
   gui_merge_goGPSbin.fig
   gui_merge_goGPSbin.m
   gui_merge_goGPSbin_unix.fig
   gui_merge_goGPSbin_unix.m
   input_bancroft.m
   input_kalman.m
   input_kalman_SA.m
   input_kalman_vinc.m
   julday.m
   kalman_goGPS_SA_cod_init.m
   kalman_goGPS_SA_cod_loop.m
   kalman_goGPS_SA_init.m
   kalman_goGPS_SA_loop.m
   kalman_goGPS_cod_init.m
   kalman_goGPS_cod_loop.m
   kalman_goGPS_init.m
   kalman_goGPS_loop.m
   kalman_goGPS_vinc_init.m
   kalman_goGPS_vinc_loop.m
   load_RINEX.m
   load_goGPSinput.m
   load_goGPSoutput.m
   load_observ.m
   load_stream.m
   lorentz.m
   NMEA_GSV_gen.m
   NMEA_PGGPK_gen.m
   NMEA_RMC_gen.m
   NMEA_checksum.m
   NTRIP_string_generator.m
   obs_type_find.m
   ref_2d_projection.m
   ref_3d_projection.m
   RINEX2goGPSbin.m
   RINEX_get_epoch.m
   RINEX_get_nav.m
   RINEX_get_nav_GLO.m
   RINEX_get_obs.m
   RINEX_parse_hdr.m
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
   streamM2RINEX.m
   streamR2RINEX.m
   streams2goGPSbin.m
   topocent.m
   twos_complement.m
   ublox_CFG_CFG.m
   ublox_CFG_MSG.m
   ublox_CFG_RATE.m
   ublox_COM_find.m
   ublox_UBX_codes.m
   ublox_check_ACK.m
   ublox_poll_message.m
   update_settings.m