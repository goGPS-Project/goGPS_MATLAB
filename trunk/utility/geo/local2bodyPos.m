function [C] = local2bodyPos(Lm, Lr1, Lr2, Lr3, Lr4)
%
% SYNTAX:
%   [B] = local2bodyPos(L, R, P, Y);
%
% INPUT:
%   L = local level position vector(s) [m]
%   R = Roll [deg]
%   P = Pitch [deg]
%   Y = Yaw [deg]
%
% OUTPUT:
%   C = body position vector(s) in b1, b2, b3 [m]
%
% DESCRIPTION:
% -> Local level uses sequence n,e,d to body system b1 = fuselage,
%    b2 = right, b3 = down.
% -> Maximum antennae 5, minimum antennae 2
%----------------------------------------------------------------------------------------------
% Hendy F. Suhandri
% Institute of Navigation
% 2013
%----------------------------------------------------------------------------------------------
%
% Note: This frame puts the b1-axis from very right side of the vehicle,
% where other receivers are located on the left side of the b1-axis (in this
% case the Y-axes from other receivers, who are not aligned with b1-axis,  
% always give negative values). Note that for other system one should be 
% aware with the signs of coordinates P2.
%
%--------------------------------------------------------------------------
%
%%find angle alpha between baseline b_mr1 and b_mr2
b_mr1  = norm(Lr1 - Lm)
b_mr2  = norm(Lr2 - Lm);
b_r1r2 = norm(Lr2 - Lr1);
alpha  = acosd((b_r1r2^2 - b_mr1^2 - b_mr2^2)/(-2*b_mr1*b_mr2))

%%body baseline coordinates for three main antennae
M = [0;0;0];
P1= [b_mr1;0;0];
P2= [b_mr2*cosd(alpha);-b_mr2*sind(alpha);0];

%%define rotation angle parameters from three previous coordinates
%heading and elevation (H,E)
H = atan2d(Lr1(2),Lr1(1));
E = atan2d(-Lr1(3),sqrt(Lr1(1)^2+Lr1(2)^2));

%bank (B)
R3 = [cosd(H) sind(H) 0;
    -sind(H) cosd(H) 0;
    0 0 1];

R2 = [cosd(E) 0 -sind(E);
    0 1 0;
    sind(E) 0 cosd(E)];

Lr2u = R3*Lr2;
Lr2uu= R2*Lr2u;

B = atand(Lr2uu(3)/Lr2uu(2));

R321=[cosd(H)*cosd(E) sind(H)*cosd(E) -sind(E);
      cosd(H)*sind(E)*sind(B)-sind(H)*cosd(B) sind(H)*sind(E)*sind(B)+cosd(H)*cosd(B) cosd(E)*sind(B);
      cosd(H)*sind(E)*cosd(B)+sind(H)*sind(B) sind(H)*sind(E)*cosd(B)-cosd(H)*sind(B) cosd(E)*cosd(B)];

P3 = R321*Lr3;

P4 = R321*Lr4;

switch nargin
    case 2
        C = [M, P1];
    case 3
        C = [M, P1, P2];
    case 4
        C = [M, P1, P2, P3];
    case 5
        C = [M, P1, P2, P3, P4];
end

