function w = anmp (a)

% normalize angle into the range -pi <= a < +pi.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

w = mod(a, 2.0 * pi);

if (abs(w) >= pi)
    
    w = w - sign(a) * (2.0 * pi);
    
end

