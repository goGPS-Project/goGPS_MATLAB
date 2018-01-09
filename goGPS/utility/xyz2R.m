function R = xyz2R(xyz)
%DESCRIPTION: given a point on earth retrun the rotation matrix with:
%           - fist lien directed to loacal east
%           - second line directed  as north
%           - third line directed as local up
[phi, lam] = cart2geod(xyz(1), xyz(2), xyz(3));

    %rotation matrix from global to local reference system
    R = [-sin(lam) cos(lam) 0;
         -sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi);
         +cos(phi)*cos(lam) +cos(phi)*sin(lam) sin(phi)];
end
