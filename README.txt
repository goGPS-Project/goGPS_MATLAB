/*******************************************\
*               ___ ___ ___                 *
*     __ _ ___ / __| _ | __|                *
*    / _` / _ \ (_ |  _|__ \                *
*    \__, \___/\___|_| |___/                *
*    |___/                    v 1.0 beta 5  *
*                                           *
\*******************************************/

   v1.0 beta    

1. Introduction
   ============
   goGPS in it's current release it's a software package designed to 
   perform GNSS static positioning in post-processing.
   It is developed in MATLAB and it is aimed at providing a tool
   useful for studying GNSS positioning, implementing and testing
   new algorithms and interacting in general with GNSS-related
   aspects.
   
   At the moment, goGPS can perform the following main tasks:
   
   - Precise Point Positioning (PPP) static stand alone processing
     using a unified least squares approach
     
   - Network Solution (NET) with two or more static receivers
     using an undifferenced unified least squares approach
     
   Both are supported by two engines, the first works with observation 
   combinations, the second use uncombined data. While the first engine
   is a bit faster and still suggested for PPP, for network solutions
   the uncombined engine is preferrable.
     
   All the corrections for precise positioning are implemented.
   Additional resources required for the GNSS processing can be
   automatically downloaded by the software itself, these include:
    - orbits
    - satellite clock corrections
    - Earth rotation parameters
    - operational manouvers
    - atmospheric loading files
    - Vienna Mapping Functions
    - Ionospheric maps

2. Requirements
   ============
   goGPS 1.0 has been developed and tested in MATLAB 2016a+
   environments, on Linux, Windows and Mac OS.
   
   The following elements are needed in order to use goGPS:

   - a computer with Windows or a UNIX-based operating system
   - a MATLAB 2016a+ installation
   
   For post-processing tasks:
   - RINEX observation file/s for the receiver/s to process    
     (RINEX files must have epochs in common for the NET solution)
     RINEX can be both in v3 or v2 format

   Additional accepted inputs:
   - RINEX meteorological file of stations close to the GNSS receivers
   - CRD file containing apriori or fixed position for some stations
   - ocean loading (blq FES2004) files
   
   Suggested additional programs (for UNIX systems)
   - gzip command (in path)
   - aria2c software (https://aria2.github.io)
   
3. Notes
   =====
   The results of the processing can be saved in various formats, 
   however at the end of processing two main variables will appear 
   into the workspace: rec and core.
   
   Core is the main goGPS class is instantiated into the object core 
   at the beginning of goGPS execution.
   
   rec (array of GNSS_Station): it is the most important object,
   it contains the list of receivers/stations needed for computations
   and to store results
   
4. Additional information
   ======================
   Have a look at: 
   https://github.com/goGPS-Project/goGPS_MATLAB/blob/goGPS_1.0_beta/docs/coding%20goGPS.pdf
   
   
   *this development is supported by Geomatics Research & Development srl
