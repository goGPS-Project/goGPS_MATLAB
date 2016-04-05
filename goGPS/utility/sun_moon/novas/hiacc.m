function hiacc

global imode lmode

lmode = imode;

if (mod (imode, 2) == 1)
    
    imode = imode - 1;
end
