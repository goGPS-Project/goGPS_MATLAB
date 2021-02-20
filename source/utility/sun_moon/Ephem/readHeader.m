function [ st, error ] = readHeader( denum )
%READHEADER reads the header file for JPL Ephemerides from hard-drive or Internet

  NIPT = 12;

  st.TTL = cell(1,3);

  st.SS = zeros(1,3);

  st.IPT = zeros(3,NIPT,'int32'); % Pointers to number of coefficients for bodies
  st.LPT = zeros(1,3,'int32'); % Pointer to number of coefficients for lunar librations
  st.RPT = zeros(1,3,'int32'); % Pointer to number of coefficients for lunar euler angle rates
  st.TPT = zeros(1,3,'int32'); % Pointer to number of coefficients for TT-TDB

  error = 0;
  debugOut = 0;

  % ***********************************************************************
  %
  %     By default, the output ephemeris will span the same interval as the
  %     input ascii data file(s).  The user may reset these to other JED's.
  %

  %      Write a fingerprint to the screen.

  if debugOut > 0
    fprintf(debugOut,' JPL ASCII-TO-DIRECT-I/O program. Last modified 15-Aug-2013.');
  end

  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  %      Read the size and number of main ephemeris records.
  %        (from header.xxx)
  if strcmp(denum,'430') || strcmp(denum,'431')
    header = sprintf('header.%s_572',denum);
  else
    header = sprintf('header.%s',denum);
  end
  hd = fopen(header,'r');
  if hd == -1
     % retrieve from NASA
    url = sprintf('%s/de%s/%s',Ephem.JPLftpSite,denum,header);
    [EPHFILE,status] = urlread(url);
    if status == 1
      lines = textscan(EPHFILE,'%s','Delimiter','\n');
      lines = lines{1};
    else
      st = cell.empty;
      error = 1;
      return;
    end
  else
    lines = textscan(hd,'%s','Delimiter','\n');
    lines = lines{1};
    fclose(hd);
  end
  pos = 1;
  KSIZE = textscan(lines{pos},'KSIZE=%d NCOEFF=  %d');
  st.NUMDE = denum;
  st.NCOEFF = KSIZE{2}(1);
  st.KSIZE = KSIZE{1}(1);
  st.TTL = cell(1,3);
  if debugOut > 0
    fprintf(debugOut,'KSIZE =%6d\n',st.KSIZE);
  end

  %      Now for the alphanumeric heading records (GROUP 1010)

  pos = NXTGRP ( lines, pos );

  if strcmp(lines{pos-1}(1:12),'GROUP   1010') == false
    ERRPRT ( 1010, 'NOT HEADER' );
    return;
  end

  st.TTL{1} = lines{pos+1};
  st.TTL{2} = lines{pos+2};
  st.TTL{3} = lines{pos+3};
  if debugOut > 0
    fprintf(debugOut,'%s\n%s\n%s\n', st.TTL{1}, st.TTL{2}, st.TTL{3});
  end

  %      Read start, end and record span  (GROUP 1030)

  pos = NXTGRP ( lines, pos );

  if strcmp(lines{pos-1}(1:12), 'GROUP   1030' ) == false
    ERRPRT ( 1030, 'NOT HEADER' );
    return;
  end

  sc = textscan (lines{pos+1},'%f %f %f');
  st.SS(1) = sc{1}(1);
  st.SS(2) = sc{2}(1);
  st.SS(3) = sc{3}(1);

  %      Read number of constants and names of constants (GROUP 1040/4).

  pos = NXTGRP ( lines, pos );

  if strcmp(lines{pos-1}(1:12), 'GROUP   1040') == false
    ERRPRT ( 1040, 'NOT HEADER' );
    return;
  end

  N = textscan (lines{pos+1},'%6d');
  st.NCON = N{1};
  st.CNAM = cell(1,st.NCON);
  pos = pos + 2;
  I = 1;
  while I<=st.NCON
    sc = textscan (lines{pos},'%s');
    for J=1:length(sc{1})
      st.CNAM{I} = sc{1}{J};
      I = I + 1;
    end
    pos = pos + 1;
  end

  %      Read number of values and values (GROUP 1041/4)

  pos = NXTGRP ( lines, pos );

  if strcmp(lines{pos-1}(1:12), 'GROUP   1041' ) == false
    ERRPRT ( 1041, 'NOT HEADER' );
    return;
  end

  N = textscan (lines{pos+1},'%d');
  N = N{1};
  if N ~= st.NCON
    ERRPRT ( 1041, 'Num Constants do not match' );
  end
  st.CVAL = zeros(1,st.NCON);
  pos = pos + 2;
  I=1;
  while I<=st.NCON
    str = textscan(lines{pos},'%f');
    for J=1:length(str{1})
      if I > st.NCON
        break;
      end
      st.CVAL(I) = str{1}(J);
      I = I + 1;
    end
    pos = pos + 1;
  end

  for  I = 1:st.NCON
    ind = strfind(st.CNAM{I},'.');
    nm = st.CNAM{I};
    if ~isempty(ind)
      for J=1:length(ind)
        nm(ind(J)) = '_';
      end
    end
    st.(nm) = st.CVAL(I);
  end

  if debugOut > 0
    for I = 1:st.NCON
      fprintf (debugOut,'%s = %24.16f\n', st.CNAM{I},st.CVAL(I));
    end
  end

  %      Zero out pointer arrays

  for I = 1:3
    for J = 1:NIPT
      st.IPT(I,J) = 0;
    end
    st.LPT(I) = 0;
    st.RPT(I) = 0;
    st.TPT(I) = 0;
  end

  %      Read pointers needed by INTERP (GROUP 1050)

  pos = NXTGRP ( lines, pos );

  if strcmp( lines{pos-1}(1:12), 'GROUP   1050' ) == false
    ERRPRT ( 1050, 'NOT HEADER' );
    return;
  end

  if debugOut > 0
    fprintf(debugOut,'\n');
  end

  line1 = lines{pos+1};
  line2 = lines{pos+2};
  line3 = lines{pos+3};
  sc{1} = textscan(line1,'%d');
  sc{2} = textscan(line2,'%d');
  sc{3} = textscan(line3,'%d');
  num = length(sc{1}{1});
  for I=1:3
    for J=1:NIPT
      st.IPT(I,J) = sc{I}{1}(J);
    end
    if num > NIPT
      st.LPT(I) = sc{I}{1}(NIPT+1);
      st.IPT(I,NIPT+1) = st.LPT(I);
    end
    if num > NIPT+1
      st.RPT(I) = sc{I}{1}(NIPT+2);
      st.IPT(I,NIPT+2) = st.RPT(I);
    end
    if num > NIPT+2
      st.TPT(I) = sc{I}{1}(NIPT+3);
      st.IPT(I,NIPT+3) = st.TPT(I);
    end
    if debugOut > 0
      for J=1:NIPT
        fprintf (debugOut,'%6d ',st.IPT(I,J));
      end
      fprintf (debugOut,'%6d %6d %6d\n',st.LPT(I),st.RPT(I),st.TPT(I));
    end
  end

  %     Read and write the ephemeris data records (GROUP 1070).

  pos = NXTGRP ( lines, pos );

  if strcmp( lines{pos}(1:12), 'GROUP   1070' ) == false
    ERRPRT(1070,'NOT HEADER');
  end

end
function ERRPRT (I, MSG)

  fprintf (2,'ERROR #%8d  %s\n',  I, MSG);

end


function pos = NXTGRP ( lines, pos )

  %     Start with nothing.

  %     The next non-blank line we encounter is a header record.
  %     The group header and data are seperated by a blank line.

  while length(lines{pos})<6 || ...
        strcmp(lines{pos}(1:5), 'GROUP') == false
    pos = pos + 1;
  end

  %     Found the header.  Read the blank line so we can get at the data.

  if strcmp(lines{pos}(1:12),'GROUP   1070') == false
    pos = pos + 1;
  end

end

