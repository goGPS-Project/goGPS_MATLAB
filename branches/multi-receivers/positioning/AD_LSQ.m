%%========================================
%%     Toolbox for attitude determination
%%     Zhen Dai
%%     dai@zess.uni-siegen.de
%%     ZESS, University of Siegen, Germany
%%     Last Modified  : 1.Sep.2008
%%========================================
%% Functions:
%%      Direct attitude determination
%% Input parameters:     
%%      mAntXYZ_local
%%          -> Antenna loca level coordinates in form of [antenna, parameter]
%%      mAntXYZ_body
%%          -> Antenna body frame coordinates in form of [antenna,parameter]
%%      yaw_init_deg,roll_init_deg,pitch_init_deg (optional)
%%          -> Starting point of LSQ adjustment
%% Output:
%%      yaw_deg,roll_deg,pitch_deg:
%%          -> Attitude parameters expressed in Euler angles
%% Remarks:
%%       1. Actually the coordinates of master antenna are not required,
%%       since they are [0 0 0] in both of local level and antenna body
%%       frame
%%       2. Starting point (initial yaw,roll,pitch) can be simply obtained
%%       from the direct attitude computation.
%%       3. Covariance matrix for each coordinate parameters are set all 1.
%%       This can be improved by introducing different weights derived from
%%       the differential positioning.
%%  Reference:
%%      Lu, Gang (1994) Development of a GPS Multi-Antenna System 
%%      for Attitude Determination. Ph.D. Thesis, Published as Report No. 20073, 
%%      Department of Geomatics Engineering, University of Calgary.
%% Remarks:
%%      In the input parameter, the yaw,roll,pitch should be in DEGREEs.

function [yaw_deg,roll_deg,pitch_deg]=AD_LSQ(mAntXYZ_local,mAntXYZ_body,yaw_init_deg,roll_init_deg,pitch_init_deg);
%% Check the input and set the innitial guess of attitude parameters
if nargin<5, 
    yaw_arc=0;
    roll_arc=0;
    pitch_arc=0;
elseif nargin==5,
    yaw_arc=yaw_init_deg*pi/180;
    roll_arc=roll_init_deg*pi/180;
    pitch_arc=pitch_init_deg*pi/180;
else 
    error('Insufficient input parameters')
end

totalantenna=size(mAntXYZ_local,1); 
 
%% mCl and mCb are the covariance matrices for local level coordinates and for antenna body
%% frame coordinates, respectively. Here we use equal weight, not only for
%% each parameter but also for each antenna. Also we neglect the covariance
%% between antennas. More precisely, we can also assign different weights
%% depending on the quality of the origional measurements
mCl=eye(3);  
mCb=eye(3); 
iter=0;
while(iter<20),
    %% Abbrevation for sin and cos
    cr=cos(roll_arc);
    sr=sin(roll_arc);
    cp=cos(pitch_arc);
    sp=sin(pitch_arc);
    cy=cos(yaw_arc);
    sy=sin(yaw_arc);

    mSum1=zeros(3,3);
    vSum2=zeros(3,1);

    for antid=1:1:totalantenna,
        clear vW mAr mR0
        %% Abbrevation for body and local frame
        xb=mAntXYZ_body(antid,1);
        yb=mAntXYZ_body(antid,2);
        zb=mAntXYZ_body(antid,3);
        xl =mAntXYZ_local(antid,1);
        yl =mAntXYZ_local(antid,2);
        zl =mAntXYZ_local(antid,3);
        %% Derivative w.r.t Euler angles
        a11=(-cr*xl-sr*sp*yl)*sy+(-sr*sp*xl+cr*yl)*cy;
        a12=(sr*zl)*sp+(-sr*sy*xl+sr*cy*yl)*cp;
        a13=(-cy*xl-sy*yl)*sr+(-sp*sy*xl+sp*cy*yl-cp*zl)*cr;
        a21=(-cp*yl)*sy-(cp*xl)*cy;
        a22=(sy*xl-cy*yl)*sp+zl*cp;
        a23=0;
        a31=(-sr*xl+cr*sp*yl)*sy+(cr*sp*xl+sr*yl)*cy;
        a32=(-cr*zl)*sp+(cr*sy*xl-cr*cy*yl)*cp;
        a33=(-sp*sy*xl+sp*cy*yl-cp*zl)*sr+(sy*yl+cy*xl)*cr;
        %% Combined rotation matrix
        mR0=[cr*cy-sr*sp*sy    cr*sy+sr*sp*cy    -sr*cp;
            -cp*sy                       cp*cy                sp;
            sr*cy+cr*sp*sy    sr*sy-cr*sp*cy     cr*cp]; 

        %% Matrix Ai
        mAi=[a11 a12 a13;a21 a22 a23; a31 a32 a33];
        %% Residual
        vW=mR0*[xl; yl; zl]-[xb; yb; zb]; 
        %% Add to matrix sum1 and sum2
        mSum1=mSum1+mAi'*inv(mR0'*mCl*mR0+mCb)*mAi;
        vSum2=vSum2+mAi'*inv(mR0'*mCl*mR0+mCb)*vW;
    end

    iter= iter+1;
    %% Correction values
    dt=-inv(mSum1)*vSum2;
    %% Update
    yaw_arc= yaw_arc+dt(1);
    pitch_arc=pitch_arc+dt(2);
    roll_arc=roll_arc+dt(3);
    %% Exit of the adjustment
    if  norm(dt)<1e-7;
        %% Format the result in degree
        yaw_deg=GetAngleDeg(yaw_arc);
        pitch_deg=GetAngleDeg(pitch_arc);
        roll_deg=GetAngleDeg(roll_arc)  ;         
        return
    end
end
%% If the update is not accomplished within 20 steps, we consider it as a
%% error. Most probably, the body frame has incorrectly fixed.
error('Least squares adjustment for attitude determination can not converge!!!')