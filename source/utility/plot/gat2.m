function [cmap] = gat(nColors, cZero)
if nargin < 1
    nColors = 255;
end

clim = caxis();
N = nColors;

L = ones(nColors,1)*100;
a = ones(nColors,1)*0;
b = ones(nColors,1)*0;

% RED     100  128  128
% VIOLET  100  128 -128
% GREEN     x -128  128
% BLUE      x -128 -128

nVal = floor(N/5);
Lab1 = [96 -17 11];
Lab2 = [86 11 80];
step = 1/(nVal);
L(1:1+nVal) = interp1([0 1]', [Lab1(1) Lab2(1)], [0:step:1]','linear');
a(1:1+nVal) = interp1([0 1]', [Lab1(2) Lab2(2)], [0:step:1]','linear');
b(1:1+nVal) = interp1([0 1]', [Lab1(3) Lab2(3)], [0:step:1]','linear');
s = 1+nVal;

nVal = floor(N/5);
Lab1 = Lab2;
Lab2 = [72 -29 31];
step = 1/(nVal+1);
L(s:s+nVal) = interp1([0 1]', [Lab1(1) Lab2(1)], [step:step:1]','pchip');
a(s:s+nVal) = interp1([0 1]', [Lab1(2) Lab2(2)], [step:step:1]','linear');
b(s:s+nVal) = interp1([0 1]', [Lab1(3) Lab2(3)], [step:step:1]','linear');
s = s+nVal;

nVal = floor(N/5);
Lab1 = Lab2;
Lab2 = [58 -33 -21];
step = 1/(nVal+1);
L(s:s+nVal) = interp1([0 1]', [Lab1(1) Lab2(1)], [step:step:1]','pchip');
a(s:s+nVal) = interp1([0 1]', [Lab1(2) Lab2(2)], [step:step:1]','linear');
b(s:s+nVal) = interp1([0 1]', [Lab1(3) Lab2(3)], [step:step:1]','linear');
s = s+nVal;

nVal = floor(N/5);
Lab1 = Lab2;
Lab2 = [45 2 -60];
step = 1/(nVal+1);
L(s:s+nVal) = interp1([0 1]', [Lab1(1) Lab2(1)], [step:step:1]','pchip');
a(s:s+nVal) = interp1([0 1]', [Lab1(2) Lab2(2)], [step:step:1]','linear');
b(s:s+nVal) = interp1([0 1]', [Lab1(3) Lab2(3)], [step:step:1]','linear');
s = s+nVal;

nVal = length(L)-s;
Lab1 = Lab2;
Lab2 = [10 27 -80];
step = 1/(nColors-s+1);
L(s:end) = interp1([0 1]', [Lab1(1) Lab2(1)], [step:step:1]','linear');
a(s:end) = interp1([0 1]', [Lab1(2) Lab2(2)], [step:step:1]','linear');
b(s:end) = interp1([0 1]', [Lab1(3) Lab2(3)], [step:step:1]','linear');

%L = [[0:100/(nColors/2-1):100]' ; flipud([0:100/(nColors/2-1):100]')];
%L = flipud([25:(100-25)/(nColors-1):100]');
cmap = flipud(Lab2RGB(L, a, b )./255);

end

function [R, G, B] = Lab2RGB(L, a, b)

    var_Y = ( L + 16 ) ./ 116;
    var_X = a ./ 500 + var_Y;
    var_Z = var_Y - b ./ 200;
    
    pos = (var_Y.^3 > 0.008856);
    var_Y(pos) = var_Y(pos).^3;
    var_Y(~pos) = (var_Y(~pos) - 16 / 116 ) ./ 7.787;
    
    pos = (var_X.^3 > 0.008856);
    var_X(pos) = var_X(pos).^3;
    var_X(~pos) = (var_X(~pos) - 16 / 116 ) ./ 7.787;
    
    pos = (var_Z.^3 > 0.008856);
    var_Z(pos) = var_Z(pos).^3;
    var_Z(~pos) = (var_Z(~pos) - 16 / 116 ) ./ 7.787;
    
    % Observer= 2?, Illuminant= D65
    ref_X =  95.047;
    ref_Y = 100.000;
    ref_Z = 108.883;
    
    X = ref_X .* var_X;
    Y = ref_Y .* var_Y;
    Z = ref_Z .* var_Z;
    
    var_X = X ./ 100;        % X from 0 to  95.047      (Observer = 2?, Illuminant = D65)
    var_Y = Y ./ 100;        % Y from 0 to 100.000
    var_Z = Z ./ 100;        % Z from 0 to 108.883
    
    var_R = var_X .*  3.2406 + var_Y .* -1.5372 + var_Z .* -0.4986;
    var_G = var_X .* -0.9689 + var_Y .*  1.8758 + var_Z .*  0.0415;
    var_B = var_X .*  0.0557 + var_Y .* -0.2040 + var_Z .*  1.0570;
    
    pos = ( var_R > 0.0031308 );
    var_R(pos) = 1.055 .* ( var_R(pos) .^ ( 1 / 2.4 ) ) - 0.055;
    var_R(~pos) = 12.92 .* var_R(~pos);
    
    pos = ( var_G > 0.0031308 );
    var_G(pos) = 1.055 .* ( var_G(pos) .^ ( 1 / 2.4 ) ) - 0.055;
    var_G(~pos) = 12.92 .* var_G(~pos);
    
    pos = ( var_B > 0.0031308 );
    var_B(pos) = 1.055 .* ( var_B(pos) .^ ( 1 / 2.4 ) ) - 0.055;
    var_B(~pos) = 12.92 .* var_B(~pos);
    
    R = max(0,min(var_R .* 255,255));
    G = max(0,min(var_G .* 255,255));
    B = max(0,min(var_B .* 255,255));
    
    if ((nargout == 1) || (nargout == 0))
        R = [R,G,B];
    end

end
