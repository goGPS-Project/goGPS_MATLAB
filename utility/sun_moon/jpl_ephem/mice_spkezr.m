%-Abstract
%
%   MICE_SPKEZR returns the state (position and velocity) of 
%   a target body relative to an observing body, optionally 
%   corrected for light  time (planetary aberration) and stellar 
%   aberration.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      targ      the scalar string name of a target body.
%                Optionally, you may supply the integer ID code 
%                for the object as an integer string, i.e. both 
%                'MOON' and '301' are legitimate strings that 
%                indicate the Moon is the target body. 
%
%                The target and observer define a state vector 
%                whose position component points from the observer 
%                to the target.
%
%      et        the scalar or 1XN-vector of double precision ephemeris epochs, 
%                expressed as seconds past J2000 TDB, at which the state of the 
%                target body relative to the observer is to be computed,
%                'et' refers to time at the observer's location
%
%      ref       the scalar string name of the reference frame relative
%                to which the output state vector should be
%                expressed. This may be any frame supported by the SPICE
%                system, including built-in frames (documented in the
%                Frames Required Reading) and frames defined by a loaded
%                frame kernel (FK). 
% 
%                When 'ref' designates a non-inertial frame, the 
%                orientation of the frame is evaluated at an epoch  
%                dependent on the selected aberration correction. 
%
%      abcorr    a scalar string that indicates the aberration corrections
%                to apply to the state of the target body to account 
%                for one-way light time and stellar aberration.
%
%                'abcorr' may be any of the following: 
%  
%                   'NONE'     Apply no correction. Return the  
%                              geometric state of the target   
%                              body relative to the observer.  
%  
%                The following values of 'abcorr' apply to the
%                "reception" case in which photons depart from the
%                target's location at the light-time corrected epoch
%                et-lt and *arrive* at the observer's location at
%                'et':
%  
%                   'LT'       Correct for one-way light time (also 
%                              called "planetary aberration") using a 
%                              Newtonian formulation. This correction 
%                              yields the state of the target at the 
%                              moment it emitted photons arriving at 
%                              the observer at 'et'. 
%  
%                              The light time correction uses an
%                              iterative solution of the light time 
%                              equation (see Particulars for details). 
%                              The solution invoked by the "LT" option 
%                              uses one iteration. 
%  
%                   'LT+S'     Correct for one-way light time and 
%                              stellar aberration using a Newtonian 
%                              formulation. This option modifies the 
%                              state obtained with the "LT" option to 
%                              account for the observer's velocity 
%                              relative to the solar system 
%                              barycenter. The result is the apparent 
%                              state of the target---the position and 
%                              velocity of the target as seen by the 
%                              observer. 
%  
%                   'CN'       Converged Newtonian light time 
%                              correction. In solving the light time 
%                              equation, the "CN" correction iterates 
%                              until the solution converges (three 
%                              iterations on all supported platforms). 
%  
%                              The "CN" correction typically does not 
%                              substantially improve accuracy because 
%                              the errors made by ignoring 
%                              relativistic effects may be larger than 
%                              the improvement afforded by obtaining 
%                              convergence of the light time solution. 
%                              The "CN" correction computation also  
%                              requires a significantly greater number 
%                              of CPU cycles than does the  
%                              one-iteration light time correction. 
%  
%                   'CN+S'     Converged Newtonian light time 
%                              and stellar aberration corrections. 
%  
%  
%                The following values of 'abcorr' apply to the 
%                "transmission" case in which photons *depart* from 
%                the observer's location at 'et' and arrive at the 
%                target's location at the light-time corrected epoch 
%                et+lt: 
%  
%                   'XLT'      "Transmission" case:  correct for 
%                              one-way light time using a Newtonian 
%                              formulation. This correction yields the 
%                              state of the target at the moment it 
%                              receives photons emitted from the 
%                              observer's location at 'et'. 
%  
%                   'XLT+S'    "Transmission" case:  correct for 
%                              one-way light time and stellar 
%                              aberration using a Newtonian 
%                              formulation  This option modifies the 
%                              state obtained with the "XLT" option to 
%                              account for the observer's velocity 
%                              relative to the solar system 
%                              barycenter. The position component of 
%                              the computed target state indicates the 
%                              direction that photons emitted from the 
%                              observer's location must be "aimed" to 
%                              hit the target. 
%  
%                   'XCN'      "Transmission" case:  converged  
%                              Newtonian light time correction. 
%  
%                   'XCN+S'    "Transmission" case:  converged  
%                              Newtonian light time and stellar  
%                              aberration corrections. 
%  
%  
%                Neither special nor general relativistic effects are 
%                accounted for in the aberration corrections applied 
%                by this routine. 
%  
%                Neither letter case or embedded blanks are significant 
%                in the 'abcorr' string. 
%
%      obs       the scalar string name of a observing body. 
%                Optionally, you may supply the integer ID code 
%                for the object as an integer string, i.e. both 
%                'MOON' and '301' are legitimate strings that 
%                indicate the Moon is the observing body.
%
%   the call:
%
%      starg = mice_spkezr(targ, et, ref, abcorr, obs)
%
%   returns:
%
%      starg   the scalar or 1xN  array of structures, each structure
%              consisting of two fields:
%
%                  state    a double precision 6x1 array containing the state
%                           of the target body in kilometers and 
%                           kilometers-per-second relative the specified
%                           observer (the first three components to of 'state' 
%                           represent the x-, y- and z-components of the 
%                           target's position, the last three components 
%                           form the corresponding velocity vector)
%
%                  lt       the double precision scalar one-way light time
%                           the observer and target in seconds; if the target 
%                           position is corrected for aberrations, then 'lt' 
%                           is the one-way light time between the observer 
%                           and the light time corrected target location
%
%              'starg' returns with the same vectorization
%              measure (N) as 'et'.
%
%      mice_spkezr also accepts the string form of the integer NAIF IDs 
%      as inputs to 'targ' and 'obs', e.g.
%
%         targ = 'Mars'
%         obs  = 'Earth'
%
%      or (remember, string representation of the integers)
%
%         targ = '499'
%         obs  = '399'
%
%      Note, If needed the user can extract the field data from vectorized
%      'starg' structures into separate arrays.
%
%      Extract the N 'state' field data to a 6XN array 'posvel':
%
%         posvel = reshape( [starg(:).state], 6, [] )
%
%      Extract the N 'lt' field data to a 1XN array 'lighttime':
%
%         lighttime = reshape( [starg(:).lt], 1, [] )
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      % 
%      %  Load a set of kernels: an SPK file, a PCK
%      %  file and a leapseconds file. Use a meta
%      %  kernel for convenience.
%      % 
%      cspice_furnsh( 'standard.tm' )
%   
%      % 
%      %  Define parameters for a state lookup:
%      % 
%      %  Return the state vector of Mars (499) as seen from
%      %  Earth (399) in the J2000 frame
%      %  using aberration correction LT+S (light time plus
%      %  stellar aberration) at the epoch 
%      %  July 4, 2003 11:00 AM PST.
%      % 
%      target   = 'Mars';
%      epoch    = 'July 4, 2003 11:00 AM PST';
%      frame    = 'J2000';
%      abcorr   = 'LT+S';
%      observer = 'Earth';
%   
%      % 
%      %  Convert the epoch to ephemeris time.
%      % 
%      et = cspice_str2et( epoch );
%      
%      % 
%      %  Look-up the state for the defined parameters.
%      % 
%      starg = mice_spkezr( target, et, frame, abcorr, observer);
%
%      % 
%      %  Output...
%      % 
%      txt = sprintf( 'The position of    : %s', target);
%      disp( txt )
%      
%      txt = sprintf( 'As observed from   : %s', observer );
%      disp( txt )
%
%      txt = sprintf( 'In reference frame : %s', frame );
%      disp( txt )
%      disp( ' ' )
%
%      % 
%      %  The first three entries of state contain the
%      %  X, Y, Z position components. The final three contain
%      %  the Vx, Vy, Vz velocity components.
%      % 
%      txt = sprintf( 'Scalar' );
%      disp( txt )
%
%      utc_epoch = cspice_et2utc( et, 'C', 3 );
%
%      txt = sprintf( 'At epoch           : %s', epoch );
%      disp( txt )
%
%      txt = sprintf( '                   : i.e. %s', utc_epoch );
%      disp( txt )
%
%      txt = sprintf( ['R (kilometers)     : ' ...
%                        '%12.4f %12.4f %12.4f'], starg.state(1:3) );
%      disp( txt )
%
%      txt = sprintf( ['V (kilometers/sec) : ' ...
%                        '%12.7f %12.7f %12.7f'], starg.state(4:6) );
%      disp( txt )
%
%      txt = sprintf( 'Light time (secs)  : %12.7f', starg.lt );
%      disp( txt )
%
%      disp(' between observer' )
%      disp(' and target' )
%      disp( ' ' )
%
%      %
%      % Create a vector of et's, starting at 'epoch'
%      % in steps of 100000 ephemeris seconds.
%      %
%      vec_et = [0:4]*100000. + et;
%
%      disp( 'Vector' )
%      vec_epoch = cspice_et2utc( vec_et, 'C', 3 );
%
%      %
%      % Look up the state vectors and light time values
%      % corresponding to the vector of input
%      % ephemeris time 'vec_et'.
%      %
%      starg = mice_spkezr( target, vec_et, frame, abcorr, observer );
%
%      for i=1:5
%
%         txt = sprintf( 'At epoch (UTC)     : %s', vec_epoch(i,:) );
%         disp( txt )
%
%         txt = sprintf( ['R (kilometers)     : ' ...
%                        '%12.4f %12.4f %12.4f'], starg(i).state(1:3) );
%         disp( txt )
%
%         txt = sprintf( ['V (kilometers/sec) : ' ...
%                        '%12.7f %12.7f %12.7f'], starg(i).state(4:6) );
%         disp( txt )
%
%         txt = sprintf( 'Light time (secs)  : %12.7f', starg(i).lt );
%         disp( txt )
%
%         disp(' between observer' )
%         disp(' and target' )
%         disp( ' ' )
%
%      end
%
%      % 
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      % 
%      cspice_kclear
%
%   MATLAB outputs:
%
%      The position of    : Mars
%      As observed from   : Earth
%      In reference frame : J2000
%       
%      Scalar
%      At epoch           : July 4, 2003 11:00 AM PST
%                         : i.e. 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3105 -27127918.9985 -18741306.3015
%      V (kilometers/sec) :   -6.8090923    7.5131823    3.0009890
%      Light time (secs)  :  269.6898814
%       between observer
%       and target
%       
%      Vector
%      At epoch (UTC)     : 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3105 -27127918.9985 -18741306.3015
%      V (kilometers/sec) :   -6.8090923    7.5131823    3.0009890
%      Light time (secs)  :  269.6898814
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 05 22:46:40.000
%      R (kilometers)     : 73140185.4144 -26390524.7797 -18446763.0348
%      V (kilometers/sec) :   -6.8317855    7.2333512    2.8893940
%      Light time (secs)  :  266.5640394
%       between observer
%       and target
%       
%      At epoch (UTC)     : 2003 JUL 07 02:33:20.000
%      R (kilometers)     : 72456239.6608 -25681031.0146 -18163339.1448
%      V (kilometers/sec) :   -6.8470343    6.9552228    2.7786326
%      Light time (secs)  :  263.4803533
%       between observer
%       and target
%       
%      At epoch (UTC)     : 2003 JUL 08 06:20:00.000
%      R (kilometers)     : 71771127.0087 -24999259.4606 -17890946.6362
%      V (kilometers/sec) :   -6.8551544    6.6789442    2.6687919
%      Light time (secs)  :  260.4395234
%       between observer
%       and target
%       
%      At epoch (UTC)     : 2003 JUL 09 10:06:40.000
%      R (kilometers)     : 71085543.8280 -24345021.1811 -17629490.7100
%      V (kilometers/sec) :   -6.8564772    6.4045794    2.5599191
%      Light time (secs)  :  257.4422002
%       between observer
%       and target
%       
%-Particulars
%
%   A sister version of this routine exists named cspice_spkezr that returns 
%   the structure field data as separate arguments.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkezr_c.
%
%   MICE.REQ
%   SPK.REQ
%   NAIF_IDS.REQ 
%   FRAMES.REQ 
%   TIME.REQ 
%
%-Version
%
%   -Mice Version 1.0.1, 22-DEC-2008, EDW (JPL)
%   
%      Header edits performed to improve argument descriptions.
%      These descriptions should now closely match the descriptions
%      in the corresponding CSPICE routine.
%
%   -Mice Version 1.0.0, 16-DEC-2005, EDW (JPL)
%
%-Index_Entries
% 
%   using body names get target state relative to an observer 
%   get state relative to observer corrected for aberrations 
%   read ephemeris data 
%   read trajectory data 
% 
%-&

function [starg] = mice_spkezr(targ, et, ref, abcorr, obs)

   switch nargin
      case 5

         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);

      otherwise

         error ( ['Usage: [_starg_] = ' ...
                  'mice_spkezr( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type 
   % return argument.
   %
   try
      [starg] = mice('spkezr_s', targ, et, ref, abcorr, obs);
   catch
      rethrow(lasterror)
   end






