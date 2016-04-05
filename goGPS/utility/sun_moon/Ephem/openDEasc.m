function [s,error] = openDEasc( s, denum, fn, dontSave )
% open an ASCII JPL planetary ephemeris file and sets initial values
% short int ephem_open (char *ephem_name,
%                       double *jd_begin, double *jd_end,
%                       short int *de_number)
% 
% ------------------------------------------------------------------------
% 
%    PURPOSE:
%       This function opens a JPL planetary ephemeris file and
%       sets initial values.  This function must be called
%       prior to calls to the other JPL ephemeris functions.
% 
%    REFERENCES:
%       Standish, E.M. and Newhall, X X (1988). "The JPL Export
%          Planetary Ephemeris"; JPL document dated 17 June 1988.
% 
%    INPUT
%    ARGUMENTS:
%       *ephem_name (char)
%          Name of the direct-access ephemeris file.
% 
%    OUTPUT
%    ARGUMENTS:
%       *jd_begin (double)
%          Beginning Julian date of the ephemeris file.
%       *jd_end (double)
%          Ending Julian date of the ephemeris file.
%       *de_number (short int)
%          DE number of the ephemeris file opened.
% 
%    RETURNED
%    VALUE:
%       (short int)
%           0   ...file exists and is opened correctly.
%           1   ...file does not exist/not found.
%           2-10...error reading from file header.
%           11  ...unable to set record length; ephemeris (DE number)
%                  not in look-up table.
% 
%    GLOBALS
%    USED:
%       SS                eph_manager.h
%       JPLAU             eph_manager.h
%       PC                eph_manager.h
%       VC                eph_manager.h
%       TWOT              eph_manager.h
%       EM_RATIO          eph_manager.h
%       BUFFER            eph_manager.h
%       IPT               eph_manager.h
%       LPT               eph_manager.h
%       NRL               eph_manager.h
%       KM                eph_manager.h
%       NP                eph_manager.h
%       NV                eph_manager.h
%       RECORD_LENGTH     eph_manager.h
%       EPHFILE           eph_manager.h
% 
%    FUNCTIONS
%    CALLED:
%       fclose            stdio.h
%       free              stdlib.h
%       fopen             stdio.h
%       fread             stdio.h
%       calloc            stdlib.h
% 
%    VER./DATE/
%    PROGRAMMER:
%       V1.0/06-90/JAB (USNO/NA)
%       V1.1/06-92/JAB (USNO/AA): Restructure and add initializations.
%       V1.2/07-98/WTH (USNO/AA): Modified to open files for different
%                                 ephemeris types. (200,403,404,405,406)
%       V1.3/11-07/WKP (USNO/AA): Updated prolog.
%       V1.4/09-10/WKP (USNO/AA): Changed ncon and denum variables and
%                                 sizeof ipt array to type 'int' for
%                                 64-bit system compatibility.
%       V1.5/09-10/WTH (USNO/AA): Added support for DE421, default case
%                                 for switch, close file on error.
%       V1.6/10-10/WKP (USNO/AA): Renamed function to lowercase to
%                                 comply with coding standards.
% 
%    NOTES:
%       KM...flag defining physical units of the output states.
%          = 1, km and km/sec
%          = 0, AU and AU/day
%       Default value is 0 (KM determines time unit for nutations.
%                           Angle unit is always radians.)
% 
% ------------------------------------------------------------------------
  if nargin < 4
    dontSave = false;
  end
  if nargin < 2
    denum = s.de_number;
  end
  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  s.de_number = 0;
  s.jd_begin = 0.0;
  s.jd_end = 0.0;
  [ s.scan, error ] = Ephem.read_asc( denum, fn, dontSave );
  if isempty(s.scan)
    return;
  end
  switch denum
    case '102'
      offset1 = 3;
    case '200'
      offset1 = 4;
    case '202'
      offset1 = 4;
    case '403'
      offset1 = 4;
    case '404'
      offset1 = 3;
    case '405'
      offset1 = 4;
    case '406'
      offset1 = 3;
    case '410'
      offset1 = 4;
    case '413'
      offset1 = 4;
    case '414'
      offset1 = 4;
    case '418'
      offset1 = 4;
    case '421'
      offset1 = 4;
    case '422'
      offset1 = 4;
    case '423'
      offset1 = 4;
    case '424'
      offset1 = 4;
    case '430'
      offset1 = 4;
    case '430t'
      offset1 = 2;
    case '431'
      offset1 = 2;
    otherwise
      %   An unknown DE file was opened. Close the file and return an error
      %   code.
      if s.output > 0
        fprintf(s.output,'openDEasc: Unknown DE%d\n',denum);
      end
      error = 11;
      return;
  end
  s.indexInc = offset1+s.RECORD_LENGTH/8;
  s.indexOffset = 1;
  len = length(s.scan);
  offset = 3;
  
  s.numRecs = len / s.indexInc;
  if mod(s.numRecs,1.0) ~= 0
    if s.output > 0
      fprintf(s.output,'openDEasc: DE%s Bad indexInc (%g) numRecs (%g) should be integer\n',...
        denum,s.indexInc,s.numRecs);
    end
    error = 12;
    return;
  end
  
  s.jd_begin = s.scan(offset); % first record's start
  s.jd_end = s.scan((s.numRecs-1)*s.indexInc+offset); % last record's start
  s.jd_inc = (s.jd_end-s.jd_begin)/(s.numRecs-1);
  s.SS(1) = s.jd_begin;
  s.SS(2) = s.jd_end;
  s.SS(3) = s.jd_inc;
  s.de_number = denum;

  error = 0;
end
