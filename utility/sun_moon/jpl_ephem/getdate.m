function [m, d, y] = getdate

% interactive request and input of calendar date

% output

%  m = calendar month
%  d = calendar day
%  y = calendar year

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for itry = 1:1:5
    fprintf('\nplease input the calendar date');

    fprintf('\n(1 <= month <= 12, 1 <= day <= 31, year = all digits!)\n');

    cdstr = input('? ', 's');

    tl = size(cdstr);

    ci = findstr(cdstr, ',');

    % extract month, day and year

    m = str2double(cdstr(1:ci(1)-1));

    d = str2double(cdstr(ci(1)+1:ci(2)-1));

    y = str2double(cdstr(ci(2)+1:tl(2)));

    % check for valid inputs

    if (m >= 1 && m <= 12 && d >= 1 && d <= 31)
        break;
    end
end
