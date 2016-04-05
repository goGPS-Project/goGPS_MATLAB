function idss = idss_novas (name)

% this function returns the id number of a solar system body
% for the version of solsys (or solsys-auxpos combination) in use.

%     name   = name of body whose id number is desired, e.g.,
%              'sun', 'moon, 'mercury', etc., expressed as all
%              upper-case letters (in)
%     idss   = id number of body, for use in calls to solsys

% note 1: in this version, only the first three letters of the
% body's name are used for identification.  alternative versions
% might use more letters.

% note 2: if name is 'jd', idss returns idss = 2 if solsys processes
% split julian dates (in successive calls), idss = 1 otherwise

% note 3: all versions of idss must return idss = -9999 for objects
% that it cannot identify or are unsupported by solsys.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

names = ['sun'; 'moo'; 'ear'; 'mer'; 'ven'; 'mar'; 'jup'; ...
         'sat'; 'ura'; 'nep'; 'plu'];

tname = cellstr(names);
     
ids = [10;    11;     3;     1;     2;     4;     5; ...
        6;     7;     8;     9];

num = 11;

idss = -9999.0;

namein = name(1:3);

% look through list of body names to find match

for i = 1:num
    
    if (strcmp(namein, tname(i)) == 1)
        
        idss = ids(i);
        
    end
    
end

% if no match, check for inquiry about split julian dates

if (strcmp(namein, 'jd ') == 1)
    
    % in this case, set idss = 2 if solsys processes split
    % julian dates (in successive calls), idss = 1 otherwise
    
    idss = 2;
    
end

