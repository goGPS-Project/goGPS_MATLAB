function [s,error] = openDEeph( s, denum, ineph, dontSave )
% open a binary JPL planetary ephemeris file and sets initial values
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
  persistent m;
  if isempty(m)
    m = containers.Map('KeyType','char','ValueType','any');
  end
  if nargin < 4
    dontSave = false;
  end
  s.de_number = 0;
  s.jd_begin = 0.0;
  s.jd_end = 0.0;
  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  error = 0;
  if isempty(ineph)
    error = 1;
    return;
  elseif exist(ineph, 'file') ~= 2
    inasc = ineph;
    inasc(1:3) = 'asc';
    Ephem.asc2eph(denum,inasc,ineph);
  end
  if exist(ineph, 'file') == 2
    % found locally
    if m.isKey(ineph)
      st = m(ineph);
    else
      EPHFILE = fopen (ineph,'rb');
      if EPHFILE ~= -1
        for i=1:3
          [st.ttl{i},count] = fread (EPHFILE, Ephem.TTLsize, '*char');
          if count ~= Ephem.TTLsize
            error = 1;
            fprintf(2,'openDEeph: Error TTL retrieved %d from "%s", expected %d\n',...
              count,ineph,Ephem.TTLsize);
          end
          st.ttl{i} = strcat(transpose(st.ttl{i}));
        end
        % pos 252
        for i=1:Ephem.OLDCONMAX
          [st.CNAM{i},count] = fread (EPHFILE, Ephem.CNAMsize, '*char');
          if count ~= Ephem.CNAMsize
            error = 1;
            fprintf(2,'openDEeph: Error CNAM retrieved %d from "%s", expected %d\n',...
              count,ineph,Ephem.CNAMsize);
          end
          st.CNAM{i} = strcat(transpose(st.CNAM{i}));
        end
        % pos 2652
        [st.SS,count] = fread (EPHFILE, 3, 'double');
        if count ~= 3
          error = 1;
          fprintf(2,'openDEeph: Error SS retrieved %d from "%s", expected %d\n',...
            count,ineph,3);
        end
        [st.NCON,count] = fread (EPHFILE, 1, 'int32');
        if count ~= 1
          error = 1;
          fprintf(2,'openDEeph: Error NCON retrieved %d from "%s", expected %d\n',...
            count,ineph,1);
        end
        [st.JPLAU,count] = fread (EPHFILE, 1, 'double');
        if count ~= 1
          error = 1;
          fprintf(2,'openDEeph: Error JPLAU retrieved %d from "%s", expected %d\n',...
            count,ineph,1);
        end
        [st.EM_RATIO,count] = fread (EPHFILE, 1, 'double');
        if count ~= 1
          error = 1;
          fprintf(2,'openDEeph: Error EM_RATIO retrieved %d from "%s", expected %d\n',...
            count,ineph,1);
        end
        % pos 2652+8*3+4+8+8
        for i=1:12
          for j=1:3
            [st.IPT(j,i),count] = fread (EPHFILE, 1, 'int32');
            if count ~= 1
              error = 1;
              fprintf(2,'openDEeph: Error IPT retrieved %d from "%s", expected %d\n',...
                count,ineph,1);
            end
          end
        end
        % pos 2652+8*3+4+8+8+12*3*4
        [st.de_number,count] = fread (EPHFILE, 1, 'int32');
        st.de_number = sprintf('%d',st.de_number);
        if count ~= 1
          error = 1;
          fprintf(2,'openDEeph: Error de_number retrieved %d from "%s", expected %d\n',...
            count,ineph,1);
        end
        % pos 2652+8*3+4+8+8+12*3*4+4
        if strcmp(st.de_number,denum(1:min(length(st.de_number),length(denum)))) == false
          if s.output > 0
            fprintf(s.output,'openDEeph: Warning DE# retrieved %s from "%s", expected %s\n',...
              st.de_number,ineph,denum);
          end
        end
        [LPT,count] = fread (EPHFILE, 3, 'int32');
        if count ~= 3
          error = 1;
          fprintf(2,'openDEeph: Error LPT retrieved %d from "%s", expected %d\n',...
            count,ineph,3);
        end
        % pos 2652+8*3+4+8+8+12*3*4+4+3*4
        for i=Ephem.OLDCONMAX+1:st.NCON
          [st.CNAM{i},count] = fread (EPHFILE, Ephem.CNAMsize, '*char');
          if count ~= Ephem.CNAMsize
            error = 1;
            fprintf(2,'openDEeph: Error CNAM retrieved %d from "%s", expected %d\n',...
              count,ineph,Ephem.CNAMsize);
          end
          st.CNAM{i} = strcat(transpose(st.CNAM{i}));
        end
        [RPT,count] = fread (EPHFILE, 3, 'int32');
        if count ~= 3
          error = 1;
          fprintf(2,'openDEeph: Error RPT retrieved %d from "%s", expected %d\n',...
            count,ineph,3);
        end
        [TPT,count] = fread (EPHFILE, 3, 'int32');
        if count ~= 3
          error = 1;
          fprintf(2,'openDEeph: Error TPT retrieved %d from "%s", expected %d\n',...
            count,ineph,3);
        end
        for i=1:3
          st.IPT(i,13) = LPT(i);
          st.IPT(i,14) = RPT(i);
          st.IPT(i,15) = TPT(i);
        end
        st.NCOEFF = Ephem.numCoeff( st.IPT );
        st.RECORD_LENGTH = 8*st.NCOEFF;
        status = fseek(EPHFILE,(2-1)*st.RECORD_LENGTH,'bof');
        if status == -1
          fprintf(2,'openDEeph seek error at %d\n',2);
        end
        [st.CVAL,count] = fread (EPHFILE, st.NCON, 'double');
        if count ~= st.NCON
          error = 1;
          fprintf(2,'openDEeph: Error CVAL retrieved %d from "%s", expected %d\n',...
            count,ineph,st.NCON);
        end
        for i=1:st.NCON
          cnm = deblank(st.CNAM{i});
          indx = strfind(cnm,'.');
          for j=1:length(indx)
            cnm(indx(j)) = '_';
          end
          if ~isempty(cnm)
            st.(cnm) = st.CVAL(i);
          else
            st.('blank') = st.CVAL(i);
          end
        end
        st.masses(Ephem.Mercury) = st.GM1;
        st.masses(Ephem.Venus) = st.GM2;
        st.masses(Ephem.Mars) = st.GM4;
        st.masses(Ephem.Jupiter) = st.GM5;
        st.masses(Ephem.Saturn) = st.GM6;
        st.masses(Ephem.Uranus) = st.GM7;
        st.masses(Ephem.Neptune) = st.GM8;
        st.masses(Ephem.Pluto) = st.GM9;
        st.masses(Ephem.Sun) = st.GMS;
        st.masses(Ephem.EarthMoonBarycenter) = st.GMB;
        st.masses(Ephem.Moon) = st.GMB/(1+s.EM_RATIO);
        st.masses(Ephem.Earth) = st.GMB*s.EM_RATIO/(1+s.EM_RATIO);
        mm = 0;
        for i=Ephem.Mercury:Ephem.Sun
          mm = mm + st.masses(i);
        end
        st.masses(Ephem.SolarSystemBarycenter) = mm;

        pos = 1;
        status = fseek(EPHFILE,(3-1)*st.RECORD_LENGTH,'bof');
        if status == -1
          if s.output > 0
            fprintf(s.output,'openDEeph: Bad fseek of "%s"\n',ineph);
          end
          fclose (EPHFILE);
          error = 9;
          return;
        end
        st.scan = zeros(1,0);
        for NROUT = 3:281474976710655
          [sc,count] = fread (EPHFILE,st.NCOEFF,'double');
          if isempty(sc) || count == 0 % expected end
            break;
          end
          if st.NCOEFF ~= count
            fprintf(2,'openDEeph: Error COEFF retrieved %d from "%s", expected %d\n',...
              count,ineph,st.NCOEFF);
            %break;
          else
            for K=1:st.NCOEFF
              st.scan(pos) = sc(K);
              pos = pos + 1;
            end
          end
        end
        fclose (EPHFILE);
        if isempty(st.scan)
          if s.output > 0
            fprintf(s.output,'openDEeph: Bad textscan of "%s"\n',ineph);
          end
          error = 9;
          return;
        end
        if dontSave == false
          m(ineph) = st;
        end
      else
        if s.output > 0
          fprintf(s.output,'openDEeph: Could not open "%s"\n',ineph);
        end
        error = 10;
        return;
      end
    end
  else
    if s.output > 0
      fprintf(2,'openDEeph: Bad file access for "%s"\n',ineph);
    end
    error = 10;
    return;
  end

  s.scan = st.scan;
  s.SS = st.SS;
  s.de_number = denum;
  s.RECORD_LENGTH = st.RECORD_LENGTH;
  s.IPT = st.IPT;
  s.JPLAU = st.JPLAU;
  s.EM_RATIO = st.EM_RATIO;
  for i=Ephem.Mercury:Ephem.EarthMoonBarycenter
    s.masses(i) = st.masses(i);
  end
  
  s.jd_begin = s.SS(1);
  s.jd_end = s.SS(2);
  s.jd_inc = s.SS(3);
  s.numRecs =  1 + (s.jd_end-s.jd_begin)/s.jd_inc;
  if mod(s.numRecs,1.0) ~= 0
    if s.output > 0
      fprintf(s.output,'openDEeph: DE%s Bad jd_inc (%g) numRecs (%g) should be integer\n',...
        denum,s.indexInc,s.numRecs);
    end
    s.jd_begin = 0.0;
    s.jd_end = 0.0;
    s.jd_inc = 0.0;
    error = 12;
    return;
  end
  s.indexInc = s.RECORD_LENGTH / 8;
  s.indexOffset = -1;
end
