%-Abstract
%
%   CSPICE_FURNSH loads SPICE kernel files into MATLAB.
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
%      file   the string scalar or NXM character array of SPICE
%             kernel file names (the kernel file may be either binary 
%             or text).
%   
%   the call:
%   
%      cspice_furnsh( file )
%   
%   loads 'file' into the SPICE kernel system. Once loaded, 
%   the subsystem processes the file based on kernel type, 
%   loading the kernels with the appropriate loader. 
%   
%   If 'file' is a SPICE meta-kernel containing initialization 
%   instructions (through use of the correct kernel pool variables), 
%   the kernel subsystem will load the listed files.
%   
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapseconds kernel.
%      %
%      cspice_furnsh( 'naif0008.tls' )
%
%   or
%   
%      %
%      % Load planetary ephemeris SPK file de405s.bsp.
%      %
%      cspice_furnsh( 'de405s.bsp' )
%
%   or
%
%      %
%      % Load a meta kernel that lists leapseconds, SPK,
%      % and PCK kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%-Particulars
%
%   Text kernels input to this routine need not have native line
%   terminators for the platform. Lower level CSPICE routines can
%   read and process non-native text files. This functionality does
%   not exist in the FORTRAN SPICELIB.
%
%   Kernel pool variable names are restricted to a length of 32 
%   characters or less.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine furnsh_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 10-FEB-2010, EDW (JPL)
%
%      Added mention of the length restriction on the kernel pool variable
%      names.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
% 
%   Load SPICE data from a list of items 
% 
%-&

function cspice_furnsh(file)

   switch nargin
      case 1

         file = zzmice_str(file);

      otherwise

         error ( 'Usage: cspice_furnsh(_`file`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('furnsh_c',file);
   catch
      rethrow(lasterror)
   end


