function [date1, leapsec] = utc2gps(date0)
%UTC2GPS Convert UTC(GMT) time tags to GPS time accounting for leap seconds
%   UTC2GPS(date) corrects an array of UTC dates(in any matlab format) for
%   leap seconds and returns an array of GPS datenums where:
%   GPS = UTC + steptime
%   Currently step times are through Jan 1 2009, but need to be added below
%   as they are instuted. All input dates must be later than the start of
%   GPS time on Jan 6 1980 00:00:00
%
%	See also GPS2UTC.

%   Copyright 2008 Ian M. Howat, ihowat@gmail.com
%   $Version: 1.0 $  $Date: 23-Aug-2008 13:56:44 $
%   Adapted by Eugenio Realini, 2013 (added 'Jul 1 2012' leap date and leap
%   seconds output)

%% ADD NEW LEAP DATES HERE:
stepdates = [...
    'Jan 6 1980'
    'Jul 1 1981'
    'Jul 1 1982'
    'Jul 1 1983'
    'Jul 1 1985'
    'Jan 1 1988'
    'Jan 1 1990'
    'Jan 1 1991'
    'Jul 1 1992'
    'Jul 1 1993'
    'Jul 1 1994'
    'Jan 1 1996'
    'Jul 1 1997'
    'Jan 1 1999'
    'Jan 1 2006'
    'Jan 1 2009'
    'Jul 1 2012'
    'Jul 1 2015'];

%% Convert Steps to datenums and make step offsets
stepdates = datenum(stepdates)'; %step date coversion
steptime = (0:length(stepdates)-1)'./86400;  %corresponding step time (sec)

%% Arg Checking
if ~isnumeric(date0) %make sure date0 are datenums, if not try converting
    date0 = datenum(date0); %will error if not a proper format
end

if ~isempty(find(date0 < stepdates(1)))%date0 must all be after GPS start date
    error('Input dates must be after 00:00:00 on Jan 6th 1980') 
end

%% Array Sizing
sz = size(date0);
date0 = date0(:);

date0 = repmat(date0,[1 size(stepdates,2)]);
stepdates = repmat(stepdates,[size(date0,1) 1]);

%% Conversion
delta = steptime(sum((date0 - stepdates) >= 0,2));
date1 = date0(:,1) + delta;

%% Reshape Output Array
date1 = reshape(date1,sz);

%% Leap seconds
leapsec = datevec(delta);
leapsec = leapsec(:,end);


