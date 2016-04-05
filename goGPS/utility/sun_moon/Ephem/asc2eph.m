function [ error ] = asc2eph( denum, INFILE, OUTFILE )
%ASC2EPH creates a binary format JPL Planetary Ephemeris file from one or more ascii text files
%
%      ASC2EPH creates a binary format JPL Planetary Ephemeris file from
%      one or more ascii text files.
%
%$ Disclaimer
%
%     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%
%      This program, 'asc2eph', requires (via standard input) an ascii
%      header file ('header.XXX'), followed by one or more ascii ephemeris 
%      data files ('ascSYYYY.XXX').  All files must have the same ephemeris
%      number, XXX.  Further, the data files must be consecutive in time
%      with no gaps between them. 
%
%      By default, the output ephemeris will span the same interval as the input
%      text file(s).  if you are interested in only a portion of data, set the
%      below T1 and T2 to the begin and end times of the span you desire.  T1
%      and T2 must be specified in  Julian Ephemeris Days (ET).
%
%      A sample sequence of files might be:
%
%        header.405  asc+1920.405 asc+1940.405 asc+1960.405 asc+1980.405  
%
%      This program is written in standard Fortran-77.  
%
% **********************************************************************************
%
%                                    *** NOTE ***
%
%      However, the units in which the length of a direct access record is specified 
%      are PROCESSOR DEPENDENT.  The parameter NRECL, the number of units per word, 
%      controls the length of a record in the direct access ephemeris.
%      The user MUST select the correct value of NRECL by editing one of the 
%      comemnted lines defining PARAMETER (NRECL) below.
%
% **********************************************************************************
%
%     Updated 02 March  2013 to accommodate more than 400 dynamical parameters
%     Updated 02 August 2013 to accommodate reading of TT-TDB
%     Updated 15 August 2013 to accommodate negative Julian dates
%
% **********************************************************************************

  % *****  The user must choose one of the following statements  *****
  %            ( Usually NRECL = 4 is used on Unix platforms)

  NRECL = 4;
  %NRECL = 1;

  % **********************************************************************************

  SS = zeros(1,3);

  FIRST = true;
  %error = 0;
  debugOut = 1;
  %debugOut = 0;

  % ***********************************************************************
  %
  %     By default, the output ephemeris will span the same interval as the
  %     input ascii data file(s).  The user may reset these to other JED's.
  %
  DB2Z = -Inf;
  T1   = -Inf;
  T2   =  Inf;

  if(NRECL ~= 1 && NRECL ~= 4)
    fprintf(2,'*** ERROR: User did not set NRECL ***\n');
    error = 1;
    return;
  end

  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  % Write a fingerprint to the screen.

  if debugOut > 0
    fprintf(debugOut,' JPL ASCII-TO-DIRECT-I/O program for DE%s. Last modified 15-Aug-2013.\n', ...
      denum);
  end

  [ st, error ] = readHeader( denum );
  if error
    fprintf(2,'*** ERROR: bad header for DE%s ***\n',denum);
    return;
  end
  IRECSZ = NRECL * st.KSIZE;
  % Open direct-access output file ('JPLEPH')
  UNIT = fopen ( OUTFILE, 'wb' );
  if UNIT == -1
    fprintf(2,'*** ERROR: bad fopen %s ***\n',OUTFILE);
    error = 1;
    return;
  end
  % write over first two records since fseek does not work until there is a
  % place to seek.  These will be filled out at the end
  for I=1:st.KSIZE
    count = fwrite(UNIT,0.0,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite 0.0 returned %d, expected 1 ***\n',count);
      fclose(UNIT);
      return;
    end
  end
  NROUT  = 0;
  DB = zeros(1,st.NCOEFF);
  if ~iscell(INFILE)
    INFILE = {INFILE};
  end
  for yr = 1:length(INFILE)
    [ scan, error ] = read_asc( denum, INFILE{yr} );
    len = length(scan);
    if len == 0 || error
      continue;
    end
    IN     = 0;
    pos = 1;
    NRW = scan(pos);
    pos = pos + 1;
    while NRW == 0.0
      NRW = scan(pos);
      pos = pos + 1;
    end
    NCOEFF = scan(pos);
    pos = pos + 1;

    for K=1:NCOEFF
      DB(K) = scan(pos);
      pos = pos + 1;
    end

    while ( ( IN    == 0 ) && ...
            ( DB(2) < T2) )

      if ( NCOEFF ~= st.NCOEFF )
        error = 1;
        fprintf (2,'ERROR #%8d  %s\n',  NCOEFF, ' 2*NCOEFF not equal to KSIZE');
        break;
      end

      % Skip this data block if the end of the interval is less
      % than the specified start time or if the it does not begin
      % where the previous block ended.

      if  ( (DB(2) >= T1) && (DB(1) >= DB2Z) )

        if ( FIRST )

          % Don't worry about the intervals overlapping
          % or abutting if this is the first applicable
          % interval.

          DB2Z  = DB(1);
          FIRST = false;
        end

        if (DB(1) ~= DB2Z )

          % Beginning of current interval is past the end
          % of the previous one.

          error = 1;
          fprintf (2,'ERROR #%8d  %s\n',  NRW, 'Records do not overlap or abut');
          break;
        end

        DB2Z  = DB(2);
        NROUT = NROUT + 1;
%         STATUS = fseek(UNIT,(NROUT+2-1)*IRECSZ,'bof');
%         if STATUS == -1
%           fprintf(2,'asc2eph seek error at %d\n',NROUT+2);
%         end
        count = fwrite (UNIT,DB,'double');

        if ( count ~= NCOEFF )
          error = 1;
          fprintf (2,'ERROR #%8d th record not written because of error %d\n', ...
            NROUT,count);
        end

        % Save this block's starting date, its interval span, and its end
        % date.

        if (NROUT == 1)
          SS(1) = DB(1);
          SS(3) = DB(2) - DB(1);
        end
        SS(2) = DB(2);

        % Update the user as to our progress every 100th block.

        if ( mod(NROUT,100) == 1 )
          if debugOut > 0
            if ( DB(1) >= T1 )
              fprintf (debugOut,'%6d EPHEMERIS RECORDS WRITTEN.  LAST JED = %12.2f\n', NROUT, DB(2));
            else
              fprintf (debugOut,' Searching for first requested record...\n');
            end
          end
        end

      end
      if pos >= len
        break;
      end
      NRW = scan(pos);
      pos = pos + 1;
      if pos >= len
        break;
      end
      while NRW == 0.0
        NRW = scan(pos);
        pos = pos + 1;
        if pos >= len
          break;
        end
      end
      if pos >= len
        break;
      end
      NCOEFF = scan(pos);
      pos = pos + 1;
      if pos >= len
        break;
      end
      for K=1:NCOEFF
        DB(K) = scan(pos);
        pos = pos + 1;
      end
    end

    if debugOut > 0
      fprintf (debugOut,'%6d EPHEMERIS RECORDS WRITTEN.  LAST JED = %12.2f\n', ...
        NROUT, DB(2));
    end
  end
  % Write header records onto output file.

  NROUT = 1;

  status = fseek(UNIT,(NROUT-1)*IRECSZ,'bof');
  if status == -1
    error = 1;
    fprintf(2,'asc2eph seek error at %d\n',NROUT);
  end
  TTL = cell(1,Ephem.TTLsize);
  for j=1:3
    for i=1:min(Ephem.TTLsize,length(st.TTL{j}))
      TTL{i} = st.TTL{j}(i);
    end
    for i=min(Ephem.TTLsize,length(st.TTL{j}))+1:Ephem.TTLsize
      TTL{i} = ' ';
    end
    count = fwrite(UNIT,char(TTL),'char');
    if count ~= Ephem.TTLsize
      error = 1;
      fprintf(2,'*** ERROR: fwrite TTL returned %d, expected %d ***\n',count,Ephem.TTLsize);
      break;
    end
  end
  CNAM = cell(1,Ephem.CNAMsize);
  if(st.NCON <= Ephem.OLDCONMAX)
    for J=1:st.NCON
      for i=1:min(Ephem.CNAMsize,length(st.CNAM{J}))
        CNAM{i} = st.CNAM{J}(i);
      end
      for i=min(Ephem.CNAMsize,length(st.CNAM{J}))+1:Ephem.CNAMsize
        CNAM{i} = ' ';
      end
      count = fwrite(UNIT,char(CNAM),'char');
      if count ~= Ephem.CNAMsize
        error = 1;
        fprintf(2,'*** ERROR: fwrite CNAM returned %d, expected %d ***\n',count,Ephem.CNAMsize);
        break;
      end
    end
    for I=st.NCON+1:Ephem.OLDCONMAX
      count = fwrite(UNIT,'      ');
      if count ~= Ephem.CNAMsize
        error = 1;
        fprintf(2,'*** ERROR: fwrite " " returned %d, expected %d ***\n',count,Ephem.CNAMsize);
        break;
      end
    end
    count = fwrite(UNIT,SS,'double');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite SS returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.NCON,'int32');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite NCON returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.AU,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite AU returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.EMRAT,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite EMRAT returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    for i=1:12
      for j=1:3
        count = fwrite (UNIT, st.IPT(j,i), 'int32');
        if count ~= 1
          error = 1;
          fprintf(2,'*** ERROR: fwrite IPT returned %d, expected %d ***\n',count,1);
          fclose(UNIT);
          return;
        end
      end
    end
    count = fwrite(UNIT,st.DENUM,'int32');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite NUMDE returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.LPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite LPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.RPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite RPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.TPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite TPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
  else
    K = Ephem.OLDCONMAX+1;
    for J=1:Ephem.OLDCONMAX
      for i=1:min(Ephem.CNAMsize,length(st.CNAM{J}))
        CNAM{i} = st.CNAM{J}(i);
      end
      for i=min(Ephem.CNAMsize,length(st.CNAM{J}))+1:Ephem.CNAMsize
        CNAM{i} = ' ';
      end
      count = fwrite(UNIT,char(CNAM),'char');
      if count ~= Ephem.CNAMsize
        error = 1;
        fprintf(2,'*** ERROR: fwrite CNAM returned %d, expected %d ***\n',count,Ephem.CNAMsize);
        fclose(UNIT);
        return;
      end
    end
    count = fwrite(UNIT,SS,'double');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite SS returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.NCON,'int32');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite NCON returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.AU,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite AU returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.EMRAT,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite EMRAT returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    for i=1:12
      for j=1:3
        count = fwrite (UNIT, st.IPT(j,i), 'int32');
        if count ~= 1
          error = 1;
          fprintf(2,'*** ERROR: fwrite IPT returned %d, expected %d ***\n',count,1);
          fclose(UNIT);
          return;
        end
      end
    end
    count = fwrite(UNIT,st.DENUM,'int32');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite DENUM returned %d, expected %d ***\n',count,1);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.LPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite LPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    for J=K:st.NCON
      for i=1:min(Ephem.CNAMsize,length(st.CNAM{J}))
        CNAM{i} = st.CNAM{J}(i);
      end
      for i=min(Ephem.CNAMsize,length(st.CNAM{J}))+1:Ephem.CNAMsize
        CNAM{i} = ' ';
      end
      count = fwrite(UNIT,char(CNAM),'char');
      if count ~= Ephem.CNAMsize
        error = 1;
        fprintf(2,'*** ERROR: fwrite CNAM returned %d, expected %d ***\n',count,Ephem.CNAMsize);
        fclose(UNIT);
        return;
      end
    end
    count = fwrite(UNIT,st.RPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite RPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
    count = fwrite(UNIT,st.TPT,'int32');
    if count ~= 3
      error = 1;
      fprintf(2,'*** ERROR: fwrite TPT returned %d, expected %d ***\n',count,3);
      fclose(UNIT);
      return;
    end
  end

  NROUT = 2;

  status = fseek(UNIT,(NROUT-1)*IRECSZ,'bof');
  if status == -1
    error = 1;
    fprintf(2,'asc2eph seek error at %d\n',NROUT);
  end
  count = fwrite(UNIT,st.CVAL,'double');
  if count ~= st.NCON
    error = 1;
    fprintf(2,'*** ERROR: fwrite CVAL returned %d, expected %d ***\n',count,st.NCON);
    fclose(UNIT);
    return;
  end
  for I=st.NCON+1:Ephem.OLDCONMAX
    count = fwrite(UNIT,0.0,'double');
    if count ~= 1
      error = 1;
      fprintf(2,'*** ERROR: fwrite 0.0 returned %d, expected %d ***\n',count,1);
      break;
    end
  end
  
  % We're through.  Wrap it up.

  fclose (UNIT);

end

