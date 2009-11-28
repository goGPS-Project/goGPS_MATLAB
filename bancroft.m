function [pos] = bancroft(B_pass)

% SYNTAX:
%   [pos] = bancroft(B_pass);
%
% INPUT:
%   B_pass = Bancroft matrix
%
% OUTPUT:
%   pos = approximated ground position (X,Y,Z coordinates)
%
% DESCRIPTION:
%   Bancroft algorithm for the computation of ground coordinates
%   having at least 4 visible satellites.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) Kai Borre 
% Kai Borre 04-30-95, improved by C.C. Goad 11-24-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

global v_light %made global for usage in goGPS

pos = zeros(4,1);

for iter = 1:2
   B = B_pass;
   [m,n] = size(B);
   for i = 1:m
      x = B(i,1);
      y = B(i,2);
      if iter == 1
         traveltime = 0.072;
      else
         z = B(i,3);
         rho = (x-pos(1))^2+(y-pos(2))^2+(z-pos(3))^2;
         traveltime = sqrt(rho)/v_light;
      end
      angle = traveltime*7.292115147e-5;
      cosa = cos(angle);
      sina = sin(angle);
      B(i,1) =	cosa*x + sina*y;
      B(i,2) = -sina*x + cosa*y;
   end; % i-loop
   
   if m > 4
      BBB = inv(B'*B)*B';
   else
      BBB = inv(B);
   end
   e = ones(m,1);
   alpha = zeros(m,1);
   for i = 1:m
      alpha(i) = lorentz(B(i,:)',B(i,:)')/2; 
   end
   BBBe = BBB*e;
   BBBalpha = BBB*alpha;
   a = lorentz(BBBe,BBBe);
   b = lorentz(BBBe,BBBalpha)-1;
   c = lorentz(BBBalpha,BBBalpha);
   root = sqrt(b*b-a*c);
   r(1) = (-b-root)/a;
   r(2) = (-b+root)/a;
   possible_pos = zeros(4,2);
   for i = 1:2
      possible_pos(:,i) = r(i)*BBBe+BBBalpha;
      possible_pos(4,i) = -possible_pos(4,i);
   end
   for j =1:m
      for i = 1:2
         c_dt = possible_pos(4,i);
         calc = norm(B(j,1:3)' -possible_pos(1:3,i))+c_dt;
         omc = B(j,4)-calc;
         abs_omc(i) = abs(omc);
      end
   end; % j-loop
   
   % discrimination between roots
   if abs_omc(1) > abs_omc(2)
      pos = possible_pos(:,2);
   else 
      pos = possible_pos(:,1);
   end
end; % iter loop
%%%%%%%%%%%%  end bancroft.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%