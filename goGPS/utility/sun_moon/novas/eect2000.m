function eect2k = eect2000 (date1, date2)

%  equation of the equinoxes complementary terms, consistent with
%  iau 2000 resolutions.

%  annexe to iers conventions 2000, chapter 5

%  capitaine, n., wallace, p.t., & mccarthy, d.d. (2003). astron. &
%  astrophys. 406, pp. 1135-1149, table 3.
%  iers conventions (2010), chapter 5, p. 60, table 5.2e.
%  (table 5.2e presented in the printed publication is a truncated
%  series. the full series, which is used in novas, is available on
%  the iers conventions center website in file tab5.2e.txt.)
%  ftp://tai.bipm.org/iers/conv2010/chapter5/

%  input

%   date1, date2 = tt date (jd = date1 + date2)

%  output

%   eect2k = complementary terms (radians)

%  this revision:  2002 november 13
%                  references updated 2010 november 26

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% 2 pi

d2pi = 6.283185307179586476925287d0;

% arc seconds to radians

das2r = 4.848136811095359935899141d-6;

% reference epoch (j2000), jd

dj0 = 2451545.0d0;

% days per julian century

djc = 36525.0d0;

%  -----------------------------------------
%  the series for the ee complementary terms
%  -----------------------------------------

%  number of terms in the series

ne0 = 33;

ne1 = 1;

%  argument coefficients for t^0

ke0 = [0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0;
    0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1;
    0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0;
    1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;];

% transpose

ke0 = ke0';

%  argument coefficients for t^1

ke1 = [0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0];

% transpose

ke1 = ke1';

%  sine and cosine coefficients for t^0

se0 = [+2640.96d-6,          -0.39d-6;
    +63.52d-6,          -0.02d-6;
    +11.75d-6,          +0.01d-6;
    +11.21d-6,          +0.01d-6;
    -4.55d-6,          +0.00d-6;
    +2.02d-6,          +0.00d-6;
    +1.98d-6,          +0.00d-6;
    -1.72d-6,          +0.00d-6;
    -1.41d-6,          -0.01d-6;
    -1.26d-6,          -0.01d-6;
    -0.63d-6,          +0.00d-6;
    -0.63d-6,          +0.00d-6;
    +0.46d-6,          +0.00d-6;
    +0.45d-6,          +0.00d-6;
    +0.36d-6,          +0.00d-6;
    -0.24d-6,          -0.12d-6;
    +0.32d-6,          +0.00d-6;
    +0.28d-6,          +0.00d-6;
    +0.27d-6,          +0.00d-6;
    +0.26d-6,          +0.00d-6;
    -0.21d-6,          +0.00d-6;
    +0.19d-6,          +0.00d-6;
    +0.18d-6,          +0.00d-6;
    -0.10d-6,          +0.05d-6;
    +0.15d-6,          +0.00d-6;
    -0.14d-6,          +0.00d-6;
    +0.14d-6,          +0.00d-6;
    -0.14d-6,          +0.00d-6;
    +0.14d-6,          +0.00d-6;
    +0.13d-6,          +0.00d-6;
    -0.11d-6,          +0.00d-6;
    +0.11d-6,          +0.00d-6;
    +0.11d-6,          +0.00d-6;];

se0 = se0';

%  sine and cosine coefficients for t^1

se1 = [-0.87d-6, +0.00d-6];

% transpose

se1 = se1';

%  interval between fundamental epoch j2000.0 and current date (jc)

t = ((date1 - dj0) + date2) / djc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fundamental arguments (from iers conventions 2000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  mean anomaly of the moon

fa(1) = anmp ((485868.249036d0 + (715923.2178d0 + (31.8792d0 ...
    + (0.051635d0 + (-0.00024470d0) * t) * t) * t) * t) * das2r ...
    + mod (1325.0d0 * t, 1.0d0) * d2pi);

%  mean anomaly of the sun

fa(2) = anmp ((1287104.793048d0 + (1292581.0481d0 + (-0.5532d0 ...
    + (+0.000136d0 + (-0.00001149d0) * t) * t) * t) * t) * das2r ...
    + mod (99.0d0 * t, 1.0d0) * d2pi);

%  mean longitude of the moon minus mean longitude of the ascending
%  node of the moon

fa(3) = anmp (( 335779.526232d0 + (295262.8478d0 + (-12.7512d0 ...
    + (-0.001037d0 + (0.00000417d0) * t) * t) * t) * t) * das2r ...
    + mod (1342d0 * t, 1d0) * d2pi);

%  mean elongation of the moon from the sun

fa(4) = anmp ((1072260.703692d0 + (1105601.2090d0 + (-6.3706d0 ...
    + (0.006593d0 + (-0.00003169d0) * t) * t) * t) * t) * das2r ...
    + mod (1236d0 * t, 1d0) * d2pi);

%  mean longitude of the ascending node of the moon

fa(5) = anmp ((450160.398036d0 + (-482890.5431d0 + (7.4722d0 ...
    + (0.007702d0 + (-0.00005939d0) * t) * t) * t) * t) * das2r ...
    + mod (-5.0d0 * t, 1.0d0) * d2pi);

fa(6) = anmp (4.402608842d0 + 2608.7903141574d0 * t);
fa(7) = anmp (3.176146697d0 + 1021.3285546211d0 * t);
fa(8) = anmp (1.753470314d0 +  628.3075849991d0 * t);
fa(9) = anmp (6.203480913d0 +  334.0612426700d0 * t);
fa(10) = anmp (0.599546497d0 +   52.9690962641d0 * t);
fa(11) = anmp (0.874016757d0 +   21.3299104960d0 * t);
fa(12) = anmp (5.481293872d0 +    7.4781598567d0 * t);
fa(13) = anmp (5.311886287d0 +    3.8133035638d0 * t);
fa(14) =      (0.024381750d0 +    0.00000538691d0 * t) * t;

%  evaluate the ee complementary terms

s0 = 0.0d0;

s1 = 0.0d0;

for i = ne0: -1: 1
    
    a = 0.0d0;
    
    for j = 1:14
        
        a = a + ke0(j, i) * fa(j);
        
    end
    
    s0 = s0 + (se0(1, i) * sin(a) + se0(2, i) * cos(a));
    
end

for i = ne1: -1: 1
    
    a = 0.0d0;
    
    for j = 1:14
        
        a = a + ke1(j, i) * fa(j);
        
    end
    
    s1 = s1 + (se1(1, i) * sin(a) + se1(2, i) * cos(a));
    
end

eect2k = (s0 + s1 * t) * das2r;

