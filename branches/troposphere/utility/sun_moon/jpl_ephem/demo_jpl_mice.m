% demo_jpl_mice.m    December 18, 2012

% this script demonstrates how to interact
% with the jpleph_mice ephemeris function

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global iephem km ephname au

% read leap seconds data file

readleap;

% de421 value for astronomical unit (kilometers)

au = 149597870.699626200;

% define target body name vector

pdata = ['Mercury                '; ...
         'Venus                  '; ...
         'Earth                  '; ...
         'Mars                   '; ...
         'Jupiter                '; ...
         'Saturn                 '; ...
         'Uranus                 '; ...
         'Neptune                '; ...
         'Pluto                  '; ...
         'Moon                   '; ...
         'Sun                    '];
   
tname = cellstr(pdata);

% define central body name vector

pdata = ['Mercury                '; ...
         'Venus                  '; ...
         'Earth                  '; ...
         'Mars                   '; ...
         'Jupiter                '; ...
         'Saturn                 '; ...
         'Uranus                 '; ...
         'Neptune                '; ...
         'Pluto                  '; ...
         'Moon                   '; ...
         'Sun                    '];
      
cname = cellstr(pdata);

% define name of ephemeris binary file

ephname = 'de421.bsp';

% initialize ephemeris function

iephem = 1;

% request output in kilometers and kilometers/second

km = 1;

% request calendar date

clc; home;
   
fprintf('\ndemo_jpl_mice - demonstrates how to use the jpleph_mice.m function\n');
   
fprintf('\nplease input a UTC calendar date\n');
   
[month, day, year] = getdate;

jdutc = julian(month, day, year);
   
% compute tdb julian date

jdtdb = utc2tdb(jdutc);

% select target body

while(1) 
    
   fprintf('\n    target body menu\n');
       
   fprintf('\n  <1> Mercury');
   fprintf('\n  <2> Venus');
   fprintf('\n  <3> Earth');
   fprintf('\n  <4> Mars');
   fprintf('\n  <5> Jupiter');
   fprintf('\n  <6> Saturn');
   fprintf('\n  <7> Uranus');
   fprintf('\n  <8> Neptune');
   fprintf('\n  <9> Pluto');
   fprintf('\n  <10> Moon');
   fprintf('\n  <11> Sun');
   
   fprintf('\n\nplease select the target body\n');
          
   itarg = input('? ');
            
   if (itarg >= 1 && itarg <= 11)
       
      break;
      
   end
   
end

% select central body

while(1)
    
   fprintf('\n    central body menu\n');

   fprintf('\n  <1> Mercury');
   fprintf('\n  <2> Venus');
   fprintf('\n  <3> Earth');
   fprintf('\n  <4> Mars');
   fprintf('\n  <5> Jupiter');
   fprintf('\n  <6> Saturn');
   fprintf('\n  <7> Uranus');
   fprintf('\n  <8> Neptune');
   fprintf('\n  <9> Pluto');
   fprintf('\n  <10> Moon');
   fprintf('\n  <11> Sun');

   fprintf('\n\nplease select the central body\n');
          
   icent = input('? ');
            
   if (icent >= 0 && icent <= 11)
       
      break;
      
   end
   
end

% evalute ephemeris
      
rrd = jpleph_mice(jdtdb, itarg, icent);

% extract position vector

r = rrd(1:3);

% extract velocity vector

v = rrd(4:6);

% print results

fprintf('\n\nprogram demo_jpl_mice\n');

fprintf('\n\nephemeris file     ');

disp(ephname);

fprintf('\ntarget body    ');

disp(tname(itarg));

fprintf('central body   ');

disp(cname(icent));

[cdstr, utstr] = jd2str(jdtdb);

fprintf('\nUTC calendar date       ');

disp(cdstr);

fprintf('\nUTC Julian date          %14.8f \n', jdutc);

fprintf('\nTDB Julian date          %14.8f \n\n', jdtdb);

fprintf('\nstate vector\n');
  
svprint(r, v);

% unload ephemeris

cspice_unload('de421.bsp');

