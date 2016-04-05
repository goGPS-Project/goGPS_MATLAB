function [ fn ] = asc_name( denum, year )
%ASC_NAME creates a naming convention following the JPL FTP site for the ascii ephemerides

  if nargin < 2
    year = datevec(now());
    year = year(1);
  end
  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  index = {};
  for i=1:length(Ephem.de_numbers)
    if strcmp(Ephem.de_numbers{i},denum)
      index = i;
      break;
    end
  end
  if isempty(index)
    fprintf(2,'asc_name: Bad DE number %s\n',denum);
    fn = '';
    return;
  end
  firstyear = Ephem.minyears(index);
  incyear = Ephem.incyears(index);
  lastyear = Ephem.maxyears(index) + incyear;
  year = year - mod(year-firstyear,incyear);
  if year == firstyear - incyear
    year = firstyear;
  elseif year == lastyear
    year = year - incyear;
  end
  if year < firstyear || year > lastyear
    fprintf(2,'asc_name: Bad year %d (%d to %d)\n',year,firstyear,lastyear);
    fn = '';
    return;
  end
  if year < 0
    prefix = 'ascm';
  else
    prefix = 'ascp';
  end
  if strcmp(denum,'430t') || strcmp(denum,'431')
    fn = sprintf('%s%05d.%s',prefix,abs(year),denum);
  else
    fn = sprintf('%s%04d.%s',prefix,abs(year),denum);
  end

end

