%%========================================
%%     Toolbox for attitude determination
%%     Zhen Dai
%%     dai@zess.uni-siegen.de
%%     ZESS, University of Siegen, Germany
%%     Last Modified  : 1.Sep.2008
%%========================================
%% Functions:
%%      Convert a angle from arc to degree
%% Input parameters:     
%%      angle_arc -> Angle in arcs.     
%% Output:
%%      angle_deg -> Angle in degrees from -180 to 180

function angle_deg=GetAngleDeg(angle_arc)

ang=angle_arc*180/pi;
 
while (ang>180) || (ang<=-180),
    if ang>180,
        ang=ang-360;
    elseif ang<-180,
        ang=ang+360;
    end
end
angle_deg=ang;
    
