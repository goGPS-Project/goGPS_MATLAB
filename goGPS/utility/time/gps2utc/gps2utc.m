function [date1, leapsec] = gps2utc(date0)
%GPS2UTC Convert GPS time tags to UTC(GMT) time accounting for leap seconds
%   GPS2UTC(date) corrects an array of GPS dates(in any matlab format) for
%   leap seconds and returns an array of UTC datenums where:
%   UTC = GPS - steptime
%   Currently step times are through Jan 1 2009, but need to be added below
%   as they are instuted. All input dates must be later than the start of
%   GPS time on Jan 6 1980 00:00:00
%
%	See also UTC2GPS.

%   Copyright 2008 Ian M. Howat, ihowat@gmail.com
%   $Version: 1.0 $  $Date: 23-Aug-2008 13:56:44 $
%   Adapted by Eugenio Realini, 2013 (added 'Jul 1 2012' leap date and leap
%   seconds output; bug-fixed 'stepdates' array)

%% ADD NEW LEAP DATES HERE:
% stepdates = [...
%     'Jan 6 1980 0:00:00'
%     'Jul 1 1981 0:00:01'
%     'Jul 1 1982 0:00:02'
%     'Jul 1 1983 0:00:03'
%     'Jul 1 1985 0:00:04'
%     'Jan 1 1988 0:00:05'
%     'Jan 1 1990 0:00:06'
%     'Jan 1 1991 0:00:07'
%     'Jul 1 1992 0:00:08'
%     'Jul 1 1993 0:00:09'
%     'Jul 1 1994 0:00:10'
%     'Jan 1 1996 0:00:11'
%     'Jul 1 1997 0:00:12'
%     'Jan 1 1999 0:00:13'
%     'Jan 1 2006 0:00:14'
%     'Jan 1 2009 0:00:15'
%     'Jul 1 2012 0:00:16'
%     'Jul 1 2015 0:00:17'];

stepdates = [7.231860000000000, ...  % 'Jan 6 1980 0:00:00'
             7.237280000115741, ...  % 'Jul 1 1981 0:00:01'
             7.240930000231481, ...  % 'Jul 1 1982 0:00:02'
             7.244580000347223, ...  % 'Jul 1 1983 0:00:03'
             7.251890000462963, ...  % 'Jul 1 1985 0:00:04'
             7.261030000578704, ...  % 'Jan 1 1988 0:00:05'
             7.268340000694444, ...  % 'Jan 1 1990 0:00:06'
             7.271990000810185, ...  % 'Jan 1 1991 0:00:07'
             7.277460000925926, ...  % 'Jul 1 1992 0:00:08'
             7.281110001041667, ...  % 'Jul 1 1993 0:00:09'
             7.284760001157407, ...  % 'Jul 1 1994 0:00:10'
             7.290250001273148, ...  % 'Jan 1 1996 0:00:11'
             7.295720001388888, ...  % 'Jul 1 1997 0:00:12'
             7.301210001504630, ...  % 'Jan 1 1999 0:00:13'
             7.326780001620371, ...  % 'Jan 1 2006 0:00:14'
             7.337740001736111, ...  % 'Jan 1 2009 0:00:15'
             7.350510001851852, ...  % 'Jul 1 2012 0:00:16'
             7.361460001967592, ...  % 'Jul 1 2015 0:00:17'
             7.366960002083334]*1e5; % 'Jan 1 2017 0:00:18'

%% Convert Steps to datenums and make step offsets
% stepdates = datenum(stepdates)'; %step date coversion
steptime = (0:length(stepdates)-1)'./86400; %corresponding step time (sec)

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
date1 = date0(:,1) - delta;

%% Reshape Output Array
date1 = reshape(date1,sz);

%% Leap seconds
leapsec = datevec(delta);
leapsec = leapsec(:,end);

