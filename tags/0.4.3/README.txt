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

   v0.4.3

1. Introduction
   ============
   goGPS is a software package designed to perform GPS positioning,
   either in post-processing or real-time.
   It is developed in MATLAB and it is aimed at providing a tool
   useful for studying GPS positioning, implementing and testing
   new algorithms and interacting in general with GPS-related
   aspects.

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
     u-blox UBX, Fastrax IT03, SkyTraq and NVS BINR* binary
	 protocols are supported) with their own drivers installed
   - GPS permanent station(s) broadcasting raw data in RTCM 3.x
     format through NTRIP protocol (at least '1002' or '1004'
     messages)

   NOTE: plotting on Google Earth requires it to be installed;
         if plotting error ellipses on Google Earth produces odd
         output on Windows, please switch Google Earth rendering
         engine to DirectX.
   
   *this development was supported by the JSPS Grant-in-Aid for Scientific Research
    (Issue No. 24700105)