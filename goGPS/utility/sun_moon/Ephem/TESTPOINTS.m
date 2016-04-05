function [ rtn ] = TESTPOINTS( denum, testCase )
%TEST test a wide range of known values from JPL's ephemerides
  output = 1;
  rtn = 0;
  
  if nargin < 1
    denum = 405;
  end
  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  if nargin < 2
    testCase = -1;
  end
  fprintf(output,'Testing binary ephemeris DE%s test points\n',denum);
  testfile = sprintf('testpo.%s',denum);
  tfd = fopen(testfile,'r');
  if tfd == -1
    % retrieve from NASA
    url = sprintf('%s/de%s/%s',Ephem.JPLftpSite,denum,testfile);
    [tfd,status] = urlread(url);
    if status == 1
      lines = textscan(tfd,'%s','Delimiter','\n');
      lines = lines{1};
    else
      str = sprintf('TEST: Error, could not retrieve %s\n',url);
      if nargin > 1
        testCase.verifyEqual(status,1,str);
      else
        fprintf(output,str);
      end
      rtn = 1;
      return;
    end
  else
    lines = textscan(tfd,'%s','Delimiter','\n');
    lines = lines{1};
    fclose(tfd);
  end
  % find start of test points
  for i=1:100
    if length(lines{i}) >= 3 && strcmp(lines{i}(1:3),'EOT')
      break;
    end
  end
  s = Ephem(denum);
  s.output = output;
  s.KM = 0;
  for i=i+1:length(lines)
    %430  1550.01.01 2287195.5 13 12  2        0.80392166046084840000
    year   = sscanf(lines{i}(4:9),'%d');
    jed    = sscanf(lines{i}(16:25),'%f');
    target = sscanf(lines{i}(26:28),'%d');
    center = sscanf(lines{i}(29:31),'%d');
    x      = sscanf(lines{i}(32:34),'%d');
    value  = sscanf(lines{i}(35:end),'%f');
    fn = Ephem.eph_name(denum,year);
    [s,error] = s.openDEeph(denum, fn);
    if error
      continue;
    end

    if target == center
      [ target_pos, target_vel, target_acc, error ] = s.state( jed, target );
    else
      [ target_pos, target_vel, target_acc, error ] = s.planet_ephemeris( jed, target, center );
    end
    if error
      continue;
    end
    CY = jed/365.25;
    if (target == Ephem.TTmTDB && x < 2) || ...
       (target == Ephem.Nutations && x < 3) || ...
       (target ~= Ephem.Nutations && target ~= Ephem.TTmTDB && x < 4)
      perr = 2.0e-7;
      tv = target_pos(x);
    else
      perr = 8.0e-6;
      if target == Ephem.TTmTDB
        tv = target_vel(x-1);
      elseif target == Ephem.Nutations
        tv = target_vel(x-2);
      else
        tv = target_vel(x-3);
      end
    end
    diff = (tv-value)/value;
    if abs(diff) > perr
      str = sprintf('%d %10.1f %d %d X(%d) %g ~= %g (%g) CY %g\n',i,jed,target,center,x, ...
        tv,value,diff,CY);
      if nargin > 1
        testCase.verifyEqual(tv,value,'RelTol',perr,str);
      else
        fprintf(output,str);
      end
      rtn = rtn + 1;
    end
  end
end
