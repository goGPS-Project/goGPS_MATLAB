function [ scan, error ] = read_asc( denum, fn, dontSave )
%READ_ASC read an ascii JPL ephemeris file from hard-drive or internet
  persistent m;
  if isempty(m)
    m = containers.Map('KeyType','char','ValueType','any');
  end
  if nargin < 3
    dontSave = false;
  end
  error = 0;
  if exist(fn, 'file') == 2
    % found locally
    if m.isKey(fn)
      scan = m(fn);
    else
      EPHFILE = fopen (fn,'r');
      if EPHFILE ~= -1
        scan = textscan(EPHFILE,'%f');
        fclose (EPHFILE);
        scan = scan{1};
        if length(scan) < 4
          fprintf(2,'read_asc: Bad textscan of %s\n',fn);
          error = 9;
        else
          if dontSave == false
            m(fn) = scan;
          end
        end
      else
        scan = zeros(0,0);
        fprintf(2,'read_asc: Could not open %s\n',fn);
        error = 10;
      end
    end
  else
    if isnumeric(denum)
      denum = sprintf('%d',denum);
    end
    % retrieve from NASA
    url = sprintf('%s/de%s/%s',Ephem.JPLftpSite,denum,fn);
    if m.isKey(url)
      scan = m(url);
    else
      [EPHFILE,status] = urlread(url);
      if status == 1 && ~isempty(EPHFILE)
        scan = textscan(EPHFILE,'%f');
        scan = scan{1};
        if isempty(scan)
          fprintf(2,'read_asc: Bad textscan of %s\n',url);
          error = 9;
        else
          if dontSave == false
            m(url) = scan;
          end
        end
      else
        scan = zeros(0,0);
        fprintf(2,'read_asc: Could not download %s\n',url);
        error = 10;
      end
    end
  end
end

