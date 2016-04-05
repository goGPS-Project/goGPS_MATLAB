function [ rtn ] = TESTZERODAY( denum, testCase )
%TEST test the known zero-day values from JPL's ephemerides
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
  [ st, error ] = Ephem.readHeader( denum );
  if error == 0 && isfield(st,'JDEPOC') % contains zero data
    year = floor(st.JDEPOC/365.2425 - 4712.019); % possibly off by +1 year at New Years Day
    s = Ephem(denum);
    s.output = output;
    s.KM = 0;
    pos_err = 5.0e-8;
    vel_err = 7.0e-6;
    for i=1:2
      if i == 1
        fprintf(output,'Testing ascii ephemeris DE%s zero day\n',denum);
        fn = Ephem.asc_name(denum,year);
        [s,error] = s.openDEasc(denum, fn);
      elseif i == 2
        fprintf(output,'Testing binary ephemeris DE%s zero day\n',denum);
        fn = Ephem.eph_name(denum,year);
        [s,error] = s.openDEeph(denum, fn);
      end
      if error
        %if nargin > 1
        %  testCase.verifyEqual(error,0,'bad openDE');
        %end
        break;
      end
      target = Ephem.Mercury;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X1')
        check(st,'1',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Venus;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X2')
        check(st,'2',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Earth;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'XB')
        check(st,'B',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Mars;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X4')
        check(st,'4',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Jupiter;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X5')
        check(st,'5',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Saturn;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X6')
        check(st,'6',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Uranus;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X7')
        check(st,'7',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Neptune;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X8')
        check(st,'8',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Pluto;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'X9')
        check(st,'9',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      target = Ephem.Moon;
      [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
      if ~error && isfield(st,'XM')
        check(st,'M',target_pos,target_vel,pos_err,vel_err,testCase);
      end
      if isfield(st,'XS') && st.XS ~= 0 && st.YS ~= 0 && st.ZS ~= 0 % has Sun data
        target = Ephem.Sun;
        [ target_pos, target_vel,~, error ] = s.state( st.JDEPOC, target );
        check(st,'S',target_pos,target_vel,pos_err,vel_err,testCase);
      end
    end
  end
end
function rtn = check(st,n,pos,vel,pos_err,vel_err,testCase)
  rtn = true;
  output = 1;
  real_pos = [st.(sprintf('X%s',n)),st.(sprintf('Y%s',n)),st.(sprintf('Z%s',n))];
  real_vel = [st.(sprintf('XD%s',n)),st.(sprintf('YD%s',n)),st.(sprintf('ZD%s',n))];
  diff = (pos(1) - real_pos(1))/real_pos(1);
  if abs(diff) >= pos_err
    str = sprintf('X%c %g ~= %g (%g)\n',n,pos(1),real_pos(1),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyEqual(pos(1),real_pos(1),'RelTol',pos_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
  diff = (pos(2) - real_pos(2))/real_pos(2);
  if abs(diff) >= pos_err
    str = sprintf('Y%c %g ~= %g (%g)\n',n,pos(2),real_pos(2),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyLessThan(abs(diff),pos_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
  diff = (pos(3) - real_pos(3))/real_pos(3);
  if abs(diff) >= pos_err
    str = sprintf('Z%c %g ~= %g (%g)\n',n,pos(3),real_pos(3),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyLessThan(abs(diff),pos_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
  diff = (vel(1) - real_vel(1))/real_vel(1);
  if abs(diff) >= vel_err
    str = sprintf('XD%c %g ~= %g (%g)\n',n,vel(1),real_vel(1),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyLessThan(abs(diff),vel_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
  diff = (vel(2) - real_vel(2))/real_vel(2);
  if abs(diff) >= vel_err
    str = sprintf('YD%c %g ~= %g (%g)\n',n,vel(2),real_vel(2),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyLessThan(abs(diff),vel_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
  diff = (vel(3) - real_vel(3))/real_vel(3);
  if abs(diff) >= vel_err
    str = sprintf('ZD%c %g ~= %g (%g)\n',n,vel(3),real_vel(3),diff);
    if nargin > 6 && testCase ~= -1
      testCase.verifyLessThan(abs(diff),vel_err,str);
    else
      fprintf(output,str);
    end
    rtn = false;
  end
end
