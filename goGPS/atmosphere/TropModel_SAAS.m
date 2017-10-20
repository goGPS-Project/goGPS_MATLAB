%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "TropModel_SAAS" is a subroutine that Computes the wet,dry*
%             and Total Tropospheric Delays Using  Saastamoinen Model          *
%USAGE:
%      [tropTOTAL,tropDRY,tropWET] = TropModel_SAAS(ReceiverPos,SatPos,... *
%                                         Temperature,Pressure,RelHumidity)*
%***The function Calls the following Subroutine(s):                        *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%                                                            RefEllipsoid);*
%4.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
 
%INPUT:                                                                    *
%      Generally, the function accepts five(5) sets of inputs:
%1.    ReceiverPos--> Receiver position in either Latitude,Longitude &     *
%                     height or ECEF(XYZ) Coordinates                      *
%2.    SatPos  -----> Satellite position(s) in ECEF(XYZ) Coordinates       *
%3.    Temperature--> Atmospheric Temperature in Degree Celsius(C)         *
%4.    Pressure  ---> Atmospheric Pressure in millibars(mbar /hPa)         *
%5.    RelHumidity--> Relative Humidity in(%) eg: 50
%Other Considerations:
%*******For four(4) inputs,default Relative Humidity of 50% is used        *
%*******For three(3)inputs,the function checks the number of columns in    *
%       Temperature and if the number of columns are three(3),then all
%       meteolorogical parameters are available or provided.In that case,  *  
%       Temperature(Temp) will be in column 1, Pressure in column 2 &      *
%       RelHumidity in column 3. i.e.[Temperature  Pressure  RelHumidity]. *
%       Any other missing column,default parameters/values are assigned. 
%       Default values are [1013.25(mbar) 291.15(K) 50(%)].
%*******For two(2)inputs,the function,examines the columns in the first two*
%       input arguements( ReceiverPos and SatPos),that is if nargin==2,and * 
%       if the columns are one(1) each in ReceiverPos and SatPos,then the  *
%       function assumes elevation angle and height respectively as inputs * 
%       where: satEL = ReceiverPos i.e. satellite elevation (satEL)        *
%                 ht = SatPos i.e. Ellipsoidal Height,                     *
%       NOTE:
%            Elevation Angles should be decimal degrees(eg:26.981.78.102)
%            HeightS are also given in meters
%       otherwise, entry will be seen as Receiver & satellite positions    *
%*******For only one(1) put, error message is given and the program        *
%       Terminates / stop working given empty matrices([]) or zeros as     * 
%       outputs                                                            *

%REMARKS:
%        All these considerations are in away, other input formats one can *
%        consider when using this subroutine.              
%NOTE:
%1.   a)Elevation Angles should be decimal degrees(eg:26.981.78.102)
%     b)HeightS are also given in meters
%2.   Receiver positions are in the ff formats:
%     a)For Geographic Coordinates (n x m):[Latitude Longitude height(h)]
%      eg:DMS:[6 40 21 -1 33 17 187.76] / DM:[6 40.35 -1 33.2833 187.76] /
%         Dec.deg:[6.6725 -1.5547 187.76]
%     b)For ECEF(XYZ) Coordinates (3 x n) / (n x 3)matrix : This is also
%       applicable to Satellite Positions(SatPos).
%     I.e.: 3xn=|[X eg:[6332942.597  6332976.932  6332957.890  6332977.582|*
%               | Y    -172955.641  -172805.004  -172878.972   -172804.786|*
%               | Z]    737935.003   737647.856   737824.057   737648.519]|*
%                ----------------------------------------------------------*
%           nx3=|[ X Y Z] eg:[6332942.597  -172955.641  737935.003 |       *
%               |             6332976.932  -172805.004  737647.856 |       *
%               |             6332957.890  -172878.972  737824.057 |       *
%               |             6332977.582  -172804.786  737648.519]|       *
%                --------------------------------------------------                                      *
%OUTPUT:                                                                   *
%1.     tropTOTAL => Total Tropospheric Delay/Error Correction in meters   *
%2.     tropDRY => Hydrostaic Tropospheric Delay in meters                 *
%3.     tropWET => slant Wet Tropospheric Delay  in meters                 *
% =========================================================================
%REFERENCE:                                                                *
%1.        Bosy J.,Borkowski J.:Troposphere Modeling in Local GPS Network  * 
%2.        J. Sanz Subirana,et al., GNSS Data Processing, Vol. I:          *
%          Fundamentals and Algorithms(ESA TM-23/1, May 2013)              *
%3.        "GPS Theory and application",edited by B.Parkinson,J.Spilker,   *
%                                                         P.Enge, AIAA,1996*
%4.        GPS Theory and Practice,Hofmann-Wellenhof, Lichtenegger, and... *
%                            Collins,Fifth, revised edition  pages 106-115.* 
%5         "Global Positioning System, Mishra & Enge", pg 172              *
%6.        Modeling of Tropospheric Delays Using ANFIS,Wayan Suparta;...   *
%                                                Kemal Maulana Alhasa(2016)*
%7.        GPS Theory,Algorithms and Applications 2ed,Guochang Xu,June 2007*
%                                                    pg 66 eqns(5.95-5.108)*
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

function [tropTOTAL,tropDRY,tropWET] = TropModel_SAAS(ReceiverPos,SatPos,Temperature,...
                                                      Pressure,RelHumidity)
switch nargin
    
    case {5,4,3,2} %When all inputs are provided
        
        if (any(nargin==[5,4,3,2]))
            
          if nargin ==5
            %Assignments
            Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h)
            Satpos=SatPos;%Satellite position in XYZ
            T=Temperature;%Atmospheric Temperature in Degree Celcius
            P=Pressure;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
            RH=RelHumidity;%Relative Humidity(%)
            
          elseif   nargin ==4
                %Assignments
                Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h)
                Satpos=SatPos;%Satellite position in XYZ
                T=Temperature;%Atmospheric Temperature in Degree Celcius
                P=Pressure;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
                RH=[];%Assign empty matrix([]) to Relative Humidity
                
          elseif  nargin==3
                %Assignments
                Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h)
                Satpos=SatPos;%Satellite position in XYZ
                
                %CHECK THE NUMBER OF COLUMNs IN Temperature
                Tcol=size(Temperature,2);
                
                switch Tcol
                    case 1 %if num of col is 1
                        T=Temperature;%Atmospheric Temperature in Degree Celcius
                        P=[];%Assign empty matrix([]) to Atmospheric Pressure
                        RH=[];%Assign empty matrix([]) to Relative Humidity
                    case 2 %if num of col is 2
                        T=Temperature(:,1);%Assign 1st Column of Temperature
                                            %Temperature in Degree Celcius
                        P=Temperature(:,2);%Assign 2nd Column of Temperature
                                           %        to Atmospheric Pressure
                        RH=[];%Assign empty matrix([]) to Relative Humidity
                    case 3 %if num of col is 3
                        T=Temperature(:,1);%Assign 1st Column of Temperature
                                            %Temperature in Degree Celcius
                        P=Temperature(:,2);%Assign 2nd Column of Temperature
                                           %        to Atmospheric Pressure
                        RH=Temperature(:,3);%Assign 3rd Column of Temperature
                                            %to Relative Humidity
                end %switch Tcol
                
          elseif nargin==2
                %Assignments
                T=[];%Assign empty matrix([]) to Temperature in Degree Celcius
                P=[];%Assign empty matrix([]) to Atmospheric Pressure
                RH=[];%Assign empty matrix([]) to Relative Humidity
              
              %CHECK IF SATELLITE ELEVATION ANGLE & HEIGHT OF RECEIVER ARE
              %                                             ENTERED INSTEAD
              %***ROUND TO THE NEAREST THOUSAND
              t1=roundmod(ReceiverPos,1000);
              t2=roundmod(SatPos,1000);
               if size(ReceiverPos,2)==1 && size(SatPos,2)==1
                  if (all(t1(:,1))==0 && all(t2(:,1))==0)
                     satEL_deg=ReceiverPos;%Assigning Rpos satEL(Satellite Elevation)
                     h=SatPos;%Assigning SatPos to ht(Elipsoidal Height)
                  elseif (all(t1(:,1))==0 && all(t2(:,1))~=0)
                        satEL_deg=ReceiverPos;%Assigning Rpos satEL(Satellite Elevation)
                        h=SatPos;%Assigning SatPos to ht(Elipsoidal Height)
                  elseif (all(t1(:,1))~=0 && all(t2(:,1))==0)
                        satEL_deg=SatPos;%Assigning SatPos satEL(Satellite Elevation)
                        h=ReceiverPos;%Assigning ReceiverPos to ht(Elipsoidal Height)
                  end %if (all(t1(:,1))==0 & all(t2(:,1))==0)
                  
               elseif size(ReceiverPos,2)==3 && size(SatPos,2)==1
                     satEL_deg=SatPos;%Assigning SatPos satEL(Satellite Elevation)
                     Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h)
               else
                   Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h)
                   Satpos=SatPos;%Satellite position in XYZ
                   
               end %if size(ReceiverPos,2)==1 && size(SatPos,2)==1
          end %if nargin ==5
        end %if (any(nargin==[5,4,3,2]))
         
        %*******IF USER INPUTs are Coordinates,.....
        if exist('Rpos','var')
        %CHECK RECEIVER POSITION TYPE(XYZ/LAT LONG h)
        [Rrow,Rcol]=size(Rpos);%Get number of rows & columns of Rpos
         %***ROUND TO THE NEAREST THOUSAND
          t3=roundmod(ReceiverPos,1000);
          
          if (Rcol>7 & Rrow==3)
             X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
             Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
             Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
             
          elseif (Rcol==7 &(Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS
                lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS
                h=Rpos(:,end);%Assigning end (7th) column to heights
                
          elseif (Rcol==7 & Rrow==3)
                if all(t3(:,1:6))==0 %if columns 1-6 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS
                  lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS
                  h=Rpos(:,end);%Assigning end (7th) column to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:6))~=0 | ...
                                                        any(t3(:,1:6))~=0)
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                      
                end %if all(t3(:,1:6))==0
                                     
          elseif (Rcol==6 && (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS
                lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                          
          elseif (Rcol==6 & Rrow==3)
                if all(t3(:,1:6))==0 %if columns 1-6 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS
                  lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                  
                elseif (all(t3(:,1:end))~=0 | any(t3(:,1:end))~=0)
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                end %if all(t3(:,1:6))==0
                
          elseif (Rcol==5 && (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM
                lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM
                 h=Rpos(:,end);%Assigning end (5th) column to  heights
                 
          elseif  (Rcol==5 & Rrow==3)
                if all(t3(:,1:4))==0 %if columns 1-4 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM
                  lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM
                   h=Rpos(:,end);%Assigning end (5th) column to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:4))~=0 | ...
                                                        any(t3(:,1:4))~=0)
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                end %if all(t3(:,1:4))==0
                
          elseif (Rcol==4 & (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM
                lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                 
          elseif (Rcol==4 & Rrow==3)
                if all(t3(:,1:4))==0 %if columns 1-4 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM
                  lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:4))~=0 | ...
                                                        any(t3(:,1:4))~=0)
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                end %if all(t3(:,1:4))==0
                
          elseif  Rcol==3
              
                if all(t3(:,1:2))==0 %if columns 1-2 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D
                  lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D
                   h=Rpos(:,end);%Assigning end (3rd) column to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:2))~=0 | ...
                                                        any(t3(:,1:2))~=0)
                      if (Rrow > 3 | Rrow < 3 )
                        X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates
                        Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates
                        Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates
                        
                      elseif (Rrow == 3 )
                            if (all(Rpos(:,1)>0) & all(Rpos(:,2)<0) &...
                                                          all(Rpos(:,3)>0))
                              X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates
                              Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates
                              Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates
                              
                            elseif (all(Rpos(1,:)>0) & all(Rpos(2,:)<0) &...
                                                          all(Rpos(3,:)>0))
                                  X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                                  Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                                  Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                                  
                            elseif  (all(Rpos(1,:)>0) & any(Rpos(2,:)<0) &...
                                                          all(Rpos(3,:)>0))
                                  X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                                  Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                                  Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                     
                            elseif (all(Rpos(:,:)>0))
                                %IDENTIFYING DATA ARRANGEMENT IN Rpos if
                                %                    ALL XYZ INPUTs ARE +VE
                                mrow=mean(Rpos,2);%Find mean of Rpos along the rows
                                mcol=mean(Rpos,1);%Find mean of Rpos along the columns
                                %ROUND EACH MEAN VALUEs to 2 sifnificant
                                %                                   Figures
                                rmrow_2=round(mrow,2,'significant');
                                rmcol_2=round(mcol,2,'significant');
                                
                                if (strcmp(num2str(rmcol_2(1)),num2str(rmcol_2(2)))...
                                        | strcmp(num2str(rmcol_2(2)),num2str(rmcol_2(3)))...
                                        | strcmp(num2str(rmcol_2(1)),num2str(rmcol_2(3))))
                                    
                                   X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                                   Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                                   Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                                   
                                elseif (all(abs(mean(rmcol_2)- rmcol_2)==0) |...
                                               (rmcol_2(1)==rmcol_2(2) | ...
                                                rmcol_2(2)==rmcol_2(3) |... 
                                                rmcol_2(1)==rmcol_2(3)))
                                            
                                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                                      
                                elseif (~strcmp(num2str(rmrow_2(1)),num2str(rmrow_2(2)))...
                                        | ~strcmp(num2str(rmrow_2(2)),num2str(rmrow_2(3)))...
                                        | ~strcmp(num2str(rmrow_2(1)),num2str(rmrow_2(3))))
                                    
                                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates

                                else
                                    X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates
                                    Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates
                                    Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates
                                    
                                end %if (strcmp(num2str(rmcol_2(1)),...
                 
                            end %if (all(Rpos(:,1)>0) || all(Rpos(:,2)<0)

                      end %if (Rrow > 3 || Rrow < 3 )
                             
                end %if all(t3(:,1:4))==0
                
          elseif  (Rcol==2 &(Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D
                lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                 
          elseif   (Rcol==2 & Rrow ==3)
                if all(t3(:,1:2))==0 %if columns 1-2 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D
                  lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                else
                    X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                    Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                    Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                end %if all(t3(:,1:2))==0
                
          elseif ((Rcol==1 & Rrow ==3) & all(t3(:,1))~=0)
                X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates
                Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates
                Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                
          end %if (Rcol>7 && Rrow==3)
          
        end %if exist('Rpos','var')
            
        %*******CONVERT Lat,Long,h IF ANY TO XYZ
        if exist('lat','var') || exist('lon','var')%Check if lat/lon exist
          %******CONVERT LAT,LONG h TO XYZ
          [X,Y,Z] = geo2xyz(lat,lon,h);
          
        elseif (exist('X','var')|| exist('Y','var') || exist('Z','var'))
              %**CONVERT USER POSITION(XYZ) TO lat,long,h
              [latRAD,longRAD,h] = xyz2LLH(X,Y,Z); %#ok<*ASGLU>
                  
        end %exist('lat','var') || exist('lon','var')
          
        if ~exist('satEL_deg','var')%If Satellite Elevation is not given,then
                                %Compute Satellite Elevation with the given
                                %Receiver & Satellite positions
         %COMPUTE SATELLITE ELEVATION(<) & AZIMUTH
         %Call the function "satAzEl"
         [satEL_deg]=satAzEl([X Y Z],Satpos);
         
        end %if ~exist('satEL','var'
        
                                                               
    otherwise
             %ISSUE ERROR MESSAGE INPUT IS ONE
              beep%Give a beep sound
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check file / Data format & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal')
              return
end %switch nargin

%************TROPOHERIC DELAY MODELING/CORRECTION USING SAASTAMOINEN MODEL
%INITIALIZING OUTPUT VARIABLEs
ZHD=zeros(size(satEL_deg,1),size(satEL_deg,2));%Assign zeros of nxm to SHD
[ZTD,ZWD]=deal(ZHD);%copy the contents of SHD to all the requested outputs
                                                  
zhd=zeros(size(satEL_deg,1),1);%Assign zeros of nx1 to ZHD
[ztd,zwd,MFh,MFw]=deal(zhd);%copy the contents of ZHD to all the requested outputs
                                                
%GET METEOROLOGICAL PARAMETERS
if all(isempty([T P RH]))
  %Use Standard atmosphere @ MSL(h=0)- Berg, 1948 (Bernese)
	Ps = 1013.25;%pressure [mbar] @ SEALEVEL
	Ts = 291.15;%temperature [K] @ SEALEVEL
	RHs = 50.0;%Relative Humidity(%)@ SEALEVEL
    
    %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
    P = Ps .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ altitude(h)
	T = Ts - 0.0065.*h; %temperature [K] @ altitude(h)
	RH = RHs * exp(-0.0006396.*h);%Relative Humidity(%)@ altitude(h)
    
elseif (~isempty(T) && all(isempty([P RH])))
      T=T+273.16; %Convert temperture degree Celcius(C) to kelvin(K)
      Ps = 1013.25;%pressure [mbar] @ SEALEVEL
	  RHs = 50.0;%Relative Humidity(RH)[%] @ SEALEVEL

      %***COMPUTE PRESSURE & RH
       P = Ps .* (1-0.0000226.*h).^5.225; %pressure [mbar] @ altitude(h)
	  RH = RHs * exp(-0.0006396.*h); %Relative Humidity(%)@ altitude(h)
     
elseif (all(~isempty([T P])) && ~isempty(RH))
      T=T+273.16; %Convert temperture degree Celcius(C) to kelvin(K)
      RHs = 50.0;%Relative Humidity(RH)[%] @ SEALEVEL
      
      %***COMPUTE PRESSURE & RH
	  RH = RHs * exp(-0.0006396.*h); %Relative Humidity(%)@ altitude(h)
    
end %if all(isempty([T P RH]))

%COMPUTE PARTIAL PRESSURE OF WATER VAPOR(e)
e = 0.01 * RH .* exp(-37.2465 + 0.213166.*T - 0.000256908.*T.^2);

%*******COMPUTE DRY,WET  TOTAL TROPO DELAYS

try
   %FIND NUMBER OF ROWS & COLUMNs IN satEL_deg
   [nrow,ncol]=size(satEL_deg);
   for i=1:ncol %Loop over the Number of Receiver positions
    
      for j=1:nrow %Loop over the Number of Satellite Elevations
    
         %1.COMPUTE ZENITH ANGLEs(Z) FROM ELEVATION ANGLEs & CONVERT TO RADIAN
         %  ------------------------------------------------------------------
         %First, Check if the elevation is less than 0, set it to .1 deg
         %REMARKs:
         %       Valid elevations are between -pi/2 and pi/2.Elevations below .1.
         %       will have the a delay mapped to .1 degree.No mapping below zero
         %       degree elevation is modeled.
						
         EL_zero = find(satEL_deg(j,i) < .0017);
         if ~isempty(EL_zero)
           satEL_deg(EL_zero) = ones(size(satEL_deg(EL_zero)))*.0017;
         end   % if ~isempty(EL_zero)

         Z =(90-satEL_deg(j,i)).* pi / 180;%Zenith Angle(s) converted to radian
   
         %1.COMPUTING DRY TROPO DELAYS(Hofmann-Wellenhof et al,2001)
         zhd(j,1)=(0.002277./cos(Z)) .* P ;
   
         %2.COMPUTING WET TROPO DELAYS(Hofmann-Wellenhof et al,2001)
         zwd(j,1)=(0.002277./cos(Z)) .* (((1255 ./ T) + 0.05).* e - tan(Z).^2 );
   
        %3.COMPUTING TOTAL TROPO DELAYS(Hofmann-Wellenhof et al,2001,Eqn 6.115)
        try
           ztd(j,1)= (0.002277./cos(Z)) .* (P + ((1255 ./ T) + 0.05).* e - tan(Z).^2 );
        catch
             ztd(j,1) =  zhd(j,1) + zwd(j,1);
        end
      end  %for j=1:nrow
    %***ARRANGE DATA FOR EACH USER POSITION/STATION
    ZHD(:,i)=zhd;  %Dry Hydrostatic tropospheric Delays
    ZWD(:,i)=zwd; %Wet tropospheric Delays
    ZTD(:,i)=ztd ; %Total tropospheric Delays
   end  %for i=1:ncol
   
%******ASSIGNMENTs
tropDRY=ZHD;%Assigning ZHD to tropDRY(dry Tropospheric delay)
tropWET= ZWD;%Assigning ZWD to tropWET(wet Tropospheric  delay)
tropTOTAL=ZTD;%Assigning ZTD to tropTOTAL(total Tropospheric delay)
    
catch
     %1.DRY TROPO DELAYS
     try
        %COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY
        %                                ACCELERATION
        f_lat_h=1+0.0026.*cos(2.*latRAD)+0.00028e-3.*h;
        tropDRY=0.002277.*f_lat_h.*P  ;
     catch
          %COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY
          %                                ACCELERATION
          f_lat_h=1-0.00266.*cos(2.*latRAD)-0.00028e-3.*h;
          tropDRY=0.0022768.*(P./f_lat_h);
     
     end %try(line 485)

    %2.WET TROPO DELAYS
    try
       tropWET=0.002277.*((1255/T)+0.05).*e;
    
    catch
         %COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY
         %                                ACCELERATION
         f_lat_h=1-0.00266.*cos(2.*latRAD)-0.00028e-3.*h;
         tropWET=(0.0022768.*((1255./T)+0.05).*e)./f_lat_h;
    end %try (line 499)
   
    %COMPUTE TOTAL TROPO DELAY(STD)
    try
       %COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY
       %                                ACCELERATION
       f_lat_h=1-0.00266.*cos(2.*latRAD)-0.00028e-3.*h;
       tropTOTAL=0.002277.*(P+(((1255./T)+0.05)).*e )./f_lat_h;
    catch
         tropTOTAL=tropDRY+tropWET;
    end  %try(line 510)
    
end %try(line 442)

%%%%%%%%%%%%%%%%%%%%%%%%END OF TropModel_Hopfield.m %%%%%%%%%%%%%%%%%%%%%%%


%1.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
function [X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid)

% *************************************************************************
% *************************************************************************
%DESCRIPTION:                                                              *
%           geo2xyz  calculate Cartesian coordinates(XYZ) given geodetic
%           coordinates Latitude(degree/dms,Longitude(degree/dms,and ...   * 
%           height(m)above reference ellipsoid(RefEllipsoid) given by the  * 
%                                                                     user.* 

%The function Calls the following function(s):                             *
%1. [a,finv] = Elipsoidpara(RefEllipsoid);                                 *

%INPUT:                                                                    *
%1      latitude -> latitude of point in decimal degrees or DMS            * 
%2.     longitude -> longitude of point in decimal degrees or DMS          *
%3.     height -> ellipsoidal height in meters(m)                          * 
%4.     RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84',...    *
%                                                       'GRS80','Pz-90')   *
%NOTE:                                                                     * 
%  if there is no height put zero(0)                                       *
% =========================================================================
%OUTPUT:                                                                   *
%1.    X --> 3D Cartesian X coordinates in meters                          *
%2.    Y --> 3D Cartesian Y coordinates in meters                          *
%3.    Z-->  3D Cartesian Z coordinates in meters                          *
% =========================================================================
%REFERENCE:*
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   *
%          and practice. 2001. Fifth revised edition. Springer,Wien,New ... 
%                                                           York.p.280-282 *
%2        GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS     *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************
%**************************************************************************

%****FORMAT DATA INPUT
switch nargin
    
    case 1
         RefEllipsoid='WGS84';
         lat=latitude;%Assigning latitude to lat
         
         %***GET THE NUMBER OF COLUMNS IN lat
         clat=size(lat,2);%No. of columns of lat matrix

         if clat==7
           latDeg  = dms2degrees(lat(:,1:3));%Latitude in DMS to decimal degrees
           longDeg = dms2degrees(lat(:,4:6));%Longitude in DMS to decimal degrees
           h       = lat(:,7);%Ellipsoidal heights in the 7TH Column
         elseif  clat==6
               latDeg  = dms2degrees(lat(:,1:3));%Latitude in DMS to decimal degrees
               longDeg = dms2degrees(lat(:,4:6));%Longitude in DMS to decimal degrees
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h)
         elseif  clat==5
               latDeg  = dm2degrees(lat(:,1:2));%Latitude in DM to decimal degrees
               longDeg = dm2degrees(lat(:,3:4));%Longitude in DM to decimal degrees
               h       = lat(:,5);%Ellipsoidal heights in the 5TH Column
         elseif  clat==4
               latDeg  = dm2degrees(lat(:,1:2));%Latitude in DM to decimal degrees
               longDeg = dm2degrees(lat(:,3:4));%Longitude in DM to decimal degrees
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h)
         elseif  clat==3
               latDeg  = lat(:,1);%Latitude in decimal degrees in the 1ST Column
               longDeg = lat(:,2);%Longitude in decimal degrees in the 2ND Column
               h       = lat(:,3);%Ellipsoidal heights in the 3RD Column
         elseif  clat==2
               latDeg  = lat(:,1);%Latitude in decimal degrees in the 1ST Column
               longDeg = lat(:,2);%Longitude in decimal degrees in the 2ND Column
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h)
         end    %if clat==7

    case {4,3,2} 
         lat  = latitude;%Assigning latitude to latDeg
         long = longitude;%Assigning longitude to lonitude
         
         %%***GET THE NUMBER OF COLUMNS IN lat & long
         clat = size(lat,2);%No. of columns of latDeg matrix
         clong= size(long,2);%No. of columns of longDeg matrix
         
         if clat==3 %if No. of columns is 2 then data is in DMS
           latDeg=dms2degrees(lat);%Latitude in DMS to decimal degrees
         elseif   clat==2 %if No. of columns is 2 then data is in DM
               latDeg=dm2degrees(lat);%Latitude in DM to decimal degrees
         else   %if No. of columns is 1 then data is in decimal degrees
             latDeg=lat;
         end     %if clat==3
    
           if clong==3 %if No. of columns is 2 then data is in DMS
             longDeg=dms2degrees(long);%Latitude in DMS to decimal degrees
    
           elseif   clong==2 %if No. of columns is 2 then data is in DM
                 longDeg=dm2degrees(long);%Latitude in DM to decimal degrees

           else   %if No. of columns is 1 then data is in decimal degrees
               longDeg=long;
           end   %if clong==3
           
           if (any(nargin==[4,3,2]))
            
          if nargin ==4
            h = height;%Assigning height to h
          elseif nargin ==3
                RefEllipsoid='WGS84';
                h = height;%Assigning height to h
          elseif nargin ==2
                RefEllipsoid='WGS84';
                h = zeros(size(latitude,1),1);%Assigning zeros to ...
                                              %       Ellipsoidal height(h)
          end %if nargin ==4
          
           end %if (any(nargin==[4,3,2]))
           
           %CHECK IF SIZE OF h CONFORMS WITH SIZE OF lat/long
           if (size(h,1)==1 && (size(lat,1) && size(long,1))>1)
            if h==0
              h=repmat(h,size(lat,1),1);%Replicating height to conform with
                                        %              size of lat or long
            else
                h=[h;repmat(zeros(size(lat,1)-1,1),1)];%Replicating height to
                                                      %conform with size of
                                                      %         lat or long
            end %if h==0
            
           elseif (size(h,1)>1 && (size(lat,1) && size(long,1))>1)
               
                 if size(h,1)<size(lat,1)
                   h=[h;repmat(zeros(size(lat,1)-size(h,1),1),1)];%Replicating height
                                                                 %to conform
                                                                 %with size of
                                                                 %lat or long
                 end %if size(h,1)<size(lat,1)
                                                                                                                                                         
           end   %if (size(h,1)==1 && (size(lat,1)&& size(long,1))>1)
       
end   % switch nargin

%ISSUE ERROR MESSAGE IF SIZE OF lat & long ARE UNEQUAL

if size(lat,1)>size(long,1) %Check sizes
   beep%Give a beep sound
   errmsg{1}=sprintf('Size of Latitude Coordinates  >  Size of Longitude Coordinates   %d > %d .',size(lat,1),size(long,1));
   errmsg{2}='';
   errmsg{3}='Please ensure that size of Latitude & Longitude Coordinates are the same .';
   errordlg(errmsg,'Coordinate(s) Input Error','modal')
   return
   
elseif size(lat,1)<size(long,1) %Check sizes
      beep%Give a beep sound
      errmsg{1}=sprintf('Size of Latitude Coordinates  <  Size of Longitude Coordinates  %d > %d .',size(lat,1),size(long,1));
      errmsg{2}='';
      errmsg{3}='Please ensure that size of Latitude & Longitude Coordinates are the same .';
      errordlg(errmsg,'Coordinate(s) Input Error','modal')
      return
      
end %if size(lat,1)>size(long,1) %Check sizes

%******CALCULATION
%
%INITIALIZE OUTPUT
X=zeros(size(latDeg,1),1);
Y=deal(X);
Z=deal(X);

%****GET ELLIPSOID PARAMETERS
%
[a,finv] = Ellipsoidpara(RefEllipsoid);

f=1/finv;%flattening

%***compute square of eccentricity
e2 = (2.*f)-((f).^2);% first eccentricity

for i=1:size(latDeg,1)
    
%***CONVERT LatDeg & LongDeg to Radians
latRad=deg2rad(latDeg(i));%converting latitude values in degrees to radian
longRad=deg2rad(longDeg(i));%converting longitude values in degrees to radian

%***COMPUTE RADIUS OF PRIME VERTICAL(N)
N = a./sqrt(1-(e2.*(sin(latRad)).^2));%prime vertical radius

%***COMPUTE DISTANCE FROM Z-AXIS(P)
P = (N + h(i)).*cos(latRad);

%COMPUTE XYZ COORDINATES
Z(i,1) = roundmod((((1-e2).*N + h(i)) .* sin(latRad)),1e-3); %cartesian Z coordinates in meters
X(i,1) = roundmod(P.*cos(longRad),1e-3); %cartesian X coordinates in meters
Y(i,1) = roundmod(P.*sin(longRad),1e-3); %cartesian Y coordinates in meters

fprintf('\n X =%12.3f\n Y =%12.3f\n Z =%12.3f\n',X(i),Y(i),Z(i));

end
%%%%%%%%%%%%%%%%%%%%%%%%END OF geo2xyz.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
%2.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
function [latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid)
% *************************************************************************+
%DESCRIPTION:                                                              +
%           xyz2LLH Converts Cartesian(XYZ)Coordinates to Geographic/...
%                                                    Ellipsoidal coords    +
%           (lat,lon,alt) on WGS-84 according to a non iterative method...
%                                                         (Bowring method) +

%INPUT:                                                                    +
%1       X                                                                 + 
%2.      Y  > vectors of Cartesian(ECEF) Coordinates  (m)                  +
%3.      Z                                                                 + 
%4.      RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84','GRS80+
%                                                                ','Pz-90')
% =========================================================================+
%OUTPUT:                                                                   +
%1.    latRAD  --> latitude in radians                                     +
%2.    longRAD --> longitude in radians                                    +
%3.          h --> ellipsoidal height in units like the input              +
%4.    latDEC  --> latitude in decimal degrees                             +
%5.    longDEC --> longitude in decimal degrees                            +
% =========================================================================+
%REFERENCE:
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   +
%          and practice. 2001. Fifth revised edition. Springer,Wien,...    +
%                                                       New York.p.280-282 +
%2         GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS
% =========================================================================+
%  Written by OSAH SAMUEL, Msc Geomatic Engineering ,2016
%       Email: osahsamuel@yahoo.ca
%         Tel:+233 (0)246137410/+233 (0)509438484
%==========================================================================+
%**************************************************************************+
%**************************************************************************+
%****FORMAT DATA INPUT

switch nargin
    case 1
          RefEllipsoid='WGS84';
          XYZ=X  ;
         [nrow,ncol]=size(XYZ);

         if (ncol==3 && (nrow>ncol || nrow<ncol ))
           X=XYZ(:,1);
           Y=XYZ(:,2);
           Z=XYZ(:,3) ;
         elseif  ((nrow==3) && (ncol>nrow || ncol<nrow ))
               X=(XYZ(1,:))';
               Y=(XYZ(2,:))';
               Z=(XYZ(3,:))';
         elseif  ((nrow==3) && (ncol==3 ))
              %***Round X to 2 digit numbers
              X_vpa=vpa(roundmod(XYZ,5),100);%Variable precision arithmetic(vpa).
              Xm=mean(X_vpa);%Mean
              if (find(X_vpa(:,1)-Xm(1)==0)| find(X_vpa(:,2)-Xm(2)==0)| find(X_vpa(:,3)-Xm(3)==0))
                X=XYZ(:,1);
                Y=XYZ(:,2);
                Z=XYZ(:,3);
              elseif (find(X_vpa(:,1)-Xm(1)~=0)| find(X_vpa(:,2)-Xm(2)~=0)| find(X_vpa(:,3)-Xm(3)~=0))
                    X=(XYZ(1,:))';
                    Y=(XYZ(2,:))';
                    Z=(XYZ(3,:))';
              end
         end
 case 2
RefEllipsoid=Y;

%***FORMAT DATA
 XYZ=X  ;
[nrow,ncol]=size(XYZ);
if (ncol==3 && (nrow>ncol || nrow<ncol ))
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3) ;
elseif ((nrow==3) && (ncol>nrow || ncol<nrow ))
X=(XYZ(1,:))';
Y=(XYZ(2,:))';
Z=(XYZ(3,:))';
elseif ((nrow==3) && (ncol==3 ))
    
%***Round X to 2 digit numbers
X_vpa=vpa(roundmod(XYZ,5),100);%Variable precision arithmetic(vpa).
Xm=mean(X_vpa);%Mean
if (find(X_vpa(:,1)-Xm(1)==0)| find(X_vpa(:,2)-Xm(2)==0)| find(X_vpa(:,3)-Xm(3)==0))
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
elseif (find(X_vpa(:,1)-Xm(1)~=0)| find(X_vpa(:,2)-Xm(2)~=0)| find(X_vpa(:,3)-Xm(3)~=0))
X=(XYZ(1,:))';
Y=(XYZ(2,:))';
Z=(XYZ(3,:))';
end
end
 case 3
     
RefEllipsoid='WGS84';
 
end
        
%******CALCULATION
%
%INITIALIZE OUTPUT
latRAD=zeros(size(X,1),1);
longRAD=deal(latRAD);
latDEC=deal(latRAD);
longDEC=deal(latRAD);
h=deal(latRAD);

%****GET ELLIPSOID PARAMETERS
%
[a,finv] = Ellipsoidpara(RefEllipsoid);

f=1/finv;% flattening
b=a*(1-f);%radus of semi-minor axis
e2=(2*f)-f.^2; %eccentrcity squared
se2=e2./(1-e2) ; % second eccentricity squared

for i=1:size(X,1)
%***COMPUTE DISTANCE IN XY PLANE FROM GEOCENTER
p=sqrt((X(i)).^2+(Y(i)).^2);
 
try
theta=atan((a.*Z(i))./(b.*p));
catch
    theta = atan2(Z(i).*a,p.*b);
end

%**COMPUTE LONGITUDE OF POINT IN RADIANS
longRAD(i,1)=atan2(Y(i),X(i));%longitude of points
% Format the longitude value
if (X(i) < 0) && (Y(i) > 0)
    longRAD(i,1)= pi + longRAD(i,1);
elseif (X(i) < 0) && (Y(i) < 0)
    longRAD(i,1) = longRAD(i,1) - pi;
end

%**COMPUTE LATITUDE OF POINT IN RADIANS
latRAD(i,1)=atan((Z(i)+(se2.*b.*((sin(theta)).^3)))./(p-(e2*a*((cos(theta)).^3))));%latitude of points

%****COMPUTE ELLIPSOIDAL HEIGHT(h)
V = a./sqrt(1-(e2.*(sin(latRAD(i,1))).^2)); % prime vertical radius of curvature(v)
h(i,1)=(p./cos(latRAD(i,1)))-V;% ellipsoidal height

%****CONVERT LATITUDE & LONGITUDE IN RADIAN TO DECIMAL DEGREES
latDEC(i,1) = rad2deg(latRAD(i,1));
longDEC(i,1)= rad2deg(longRAD(i,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%END OF xyz2LLH.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

function [a,finv,RefEllipsoid] = Ellipsoidpara(RefEllipsoid)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            Elipsoidpara Returns the radius of the Semi Major Axis(a) and * 
%            inverse flattening of a given Ellipsoid                       * 
%INPUT:                                                                    *
%      RefEllipsoid -The name of the Ellipsoid in single quote(eg:'WGS84') *

%OUTPUT:                                                                   *
%       a    = Semi Major Axis of the given ellipsoid in meters            *
%       finv = inverse flattening of the given ellipsoid                   *

%REFERENCE:                                                                *
%         Department of Defence World Geodetic System 1984                 *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *                                       
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

if ischar(RefEllipsoid)
 
  refE=RefEllipsoid;%Refernce Ellipsoid

 %***GPS REFERENCE FRAME
if (strcmpi(refE,'WGS84')|| (strcmpi(refE,'WGS1984')||strcmpi(refE,'GPS')))
RefEllipsoid='WGS84';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257223563;%inverse flattening on the WGS 84 ellipsoid

elseif (strcmpi(refE,'WGS72')|| strcmpi(refE,'WGS1972') )
RefEllipsoid='WGS72';
a=6378135.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.26;%inverse flattening on the WGS 84 ellipsoid

%***Geodetic Reference System 1980
elseif (strcmpi(refE,'GRS80')|| strcmpi(refE,'GRS1980') )
RefEllipsoid='GRS80';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid

%Geodetic Reference System 1967
elseif (strcmpi(refE,'GRS67') || strcmpi(refE,'GRS1967'))
RefEllipsoid='GRS67';
a=6378160.0;%radius of major axis
finv=298.247167427;%inverse flattening
 
%GLONASS REFERENC FRAME
elseif (strcmpi(refE,'Pz-90')||(strcmpi(refE,'Pz-90.02')||strcmpi(refE,'GLONASS')))
RefEllipsoid='Pz-90';
a=6378136.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257839303;%inverse flattening on the WGS 84 ellipsoid

%***BEIDOU REFERENC FRAME
elseif (strcmpi(refE,'CGCS2000')|| strcmpi(refE,'BEIDOU') )
RefEllipsoid='CGCS2000';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid

%***GHANA
%OLD
elseif (strcmpi(refE,'WAROFFICE') || strcmpi(refE,'WO'))
RefEllipsoid='WAROFFICE';
a=6378299.996;%radius of major axis
finv=296;%inverse flattening

elseif (strcmpi(refE,'ACCRADATUM') || strcmpi(refE,'AD'))
RefEllipsoid='WAROFFICE';
a=6378299.996;%radius of major axis
finv=296;%inverse flattening

elseif (strcmpi(refE,'LEIGONDATUM') || strcmpi(refE,'LD'))
RefEllipsoid='Clarke1880';
a=20926202 ;%radius of major axis
finv=293.465;%inverse flattening

%NEW DATUM(Ghana Geodetic Datum)
elseif (strcmpi(refE,'GGD') || strcmpi(refE,'Ghana Geodetic Datum'))
RefEllipsoid='GRS80'; %Geodetic Reference System 1980
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid
    
elseif (strcmpi(refE,'Clarke1866') || strcmpi(refE,'Clarke66'))
RefEllipsoid='Clarke66';
a=6378206.4;%radius of major axis
finv=294.9786982; %inverse flattening

elseif (strcmpi(refE,'Clarke1880') || strcmpi(refE,'Clarke80'))
RefEllipsoid='Clarke80';
a=6378249.145;%radius of major axis
finv=293.465;%inverse flattening

elseif (strcmpi(refE,'Airy1830') || strcmpi(refE,'Airy30'))
RefEllipsoid='Airy30';
a=6377563.396;%radius of major axis
finv=299.3249646;%inverse flattening

elseif (strcmpi(refE,'Everest30') || strcmpi(refE,'Everest1830'))
RefEllipsoid='Everest30';
a=6377276.0;%radius of major axis
finv=300.8;%inverse flattening
  
elseif (strcmpi(refE,'Bessel41') || strcmpi(refE,'Bessel1841'))
RefEllipsoid='Bessel41';
a=6377397.0;%radius of major axis
finv=299.15;%inverse flattening
end

else
RefEllipsoid='WGS84';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257223563;%inverse flattening on the WGS 84 ellipsoid
end
%%%%%%%%%%%%%%%%%%%%%%%%END OF Elipsoidpara.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4.******SUBROUTINE TO COMPUTE SATELLITE AZIMUTH & ELEVATION ANGLE

function [satEL_deg,satEL_rad,satAZ_deg,satAZ_rad,Range]=satAzEl(UserXYZ,SatXYZ,...
                                                              RefEllipsoid)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "satAzEl" Computes Satellites Azimuth(AZ),Elevation Angle(EL) *
%             and Range(D),relative to observer position/Ground Station /..*
%             from User to Satellite given User and Satellite Position(s)  *   
  
%The function Calls the following Subroutine(s):                             *
%1. [a,finv] = Elipsoidpara(RefEllipsoid);                                 *
%2. [latRAD,longRAD,h] = xyz2LLH(X,Y,Z,RefEllipsoid);                      *

%INPUT:                                                                    *
%      UserXYZ - Observer position/Ground station coordinates(ECEF)in(m) ..*
%                                                  (3 x n) / (n x 3)matrix *
%      SatXYZ  - ECEF satellite coordinates (m), (3 x n) / (n x 3) matrix  *

%     I.e.: 3xn=|[X eg:[6332942.597  6332976.932  6332957.890  6332977.582|*
%               | Y    -172955.641  -172805.004  -172878.972   -172804.786|*
%               | Z]    737935.003   737647.856   737824.057   737648.519]|*
%                ----------------------------------------------------------*
%           nx3=|[ X Y Z] eg:[6332942.597  -172955.641  737935.003 |       *
%               |             6332976.932  -172805.004  737647.856 |       *
%               |             6332957.890  -172878.972  737824.057 |       *
%               |             6332977.582  -172804.786  737648.519]|       *
%                --------------------------------------------------        *
%      RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84','GRS80',*
%                                                                   'Pz-90)*
%NOTE:
%      If RefEllipsoid is not indicated,system uses 'WGS84' as default     *

%OUTPUT:                                                                   *
%       satAZ_rad = Azimuth (radians ClockWise from North)                 *
%		satEL_rad = Elevation Angle (radians)                              *
%       satAZ_deg = Azimuth (degrees ClockWise from North)                 *
%		satEL_deg = Elevation Angle (degrees)                              *
%		    Range = Receiver-Satellite Distance in units like the input    *

% REMARKS:                                                                 *
%         All the catch routins in this m-file are alternative solutions to*
%         their corresponding try(s).all give the same results.                               *

% =========================================================================
%REFERENCE:                                                                *
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   *
%          and practice. 2001. Fifth revised edition. Springer, Wien, New..*
%                                                            York.p.280-282*
%2         GPS Theory and application",edited by B.Parkinson,J.Spilker     *
%3.        Alfred Leick, GPS Satellite Surveying, 2nd ed.,...              *
%          Wiley-Interscience,John Wiley & Sons, New York, 1995.           *                     * 
%4.        GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS    *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *                                      
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

%****FORMAT DATA INPUT

switch nargin
    
    case 2
               
RefEllipsoid='WGS84';
          
end

[nrow,ncol]=size(UserXYZ);

if(ncol==3 && (nrow>ncol || nrow<ncol ))

UserXYZ=UserXYZ';

end
[nrow,ncol]=size(SatXYZ);

if(ncol==3 && (nrow>ncol || nrow<ncol ))

SatXYZ=SatXYZ';

end

%INITIALIZING OUTPUT VARIABLEs
satAZ_rad=zeros(size(SatXYZ,2),size(UserXYZ,2));%Assign zeros of nxm to satAZ_rad
[satEL_rad,satAZ_deg,satEL_deg,Range]=deal(satAZ_rad);%copy the contents of satAZ_rad
                                                  %to all the requested outputs
Az_rad=zeros(size(SatXYZ,2),1);%Assign zeros of nx1 to Az_rad
[El_rad,Az_deg,El_deg,hdist,range]=deal(Az_rad);%copy the contents of Az_rad
                                                %to all the requested outputs
for j=1:size(UserXYZ,2)%LOOP THROUGH ALL USER POSITIONs
    
for i=1:size(SatXYZ,2)%LOOP THROUGH ALL SATELLITE POSITIONs

%**CONVERT USER POSITION(XYZ) TO latitude and longitude
[lat,long] = xyz2LLH(UserXYZ(:,j),RefEllipsoid);

%FINDING Baseline differences
dxyz = SatXYZ(:,i)-UserXYZ(:,j);

try

%*****COMPUTING Rotation matrix
%REMARKS:
%      R is direction cosine matrix to transform position vector components
%      from geocentric equatorial frame into the topocentric horizon frame.

try
R = [-sin(long)  -sin(lat)*cos(long)  cos(lat)*cos(long);
      cos(long)  -sin(lat)*sin(long)  cos(lat)*sin(long);
       0              cos(lat)            sin(lat)     ];
 %**COMPUTE UNIT VECTORS
  vENU = R'*dxyz;
   
catch
   
 R = [   -sin(long)                        cos(long)         0     ;
              -sin(lat)*cos(long)  -sin(lat)*sin(long)  cos(lat);
              cos(lat)*cos(long)   cos(lat)*sin(long)  sin(lat)];
          
%**COMPUTE UNIT VECTORS
vENU = R*dxyz ;
        
end

%***EXTRACTING INDIVIDUAL LOCAL COORDINATES
E = vENU(1);% 'East'-coordinate relative to local origin (meters)
N = vENU(2);% 'North'-coordinate relative to local origin (meters)
U = vENU(3);% Up-coordinate relative to local origin (meters)

%***COMPUTING AZIMUTH(Az) AND ELEVATION ANGLE(El)FROM NEU
hdist(i,1) = sqrt(E.^2+N.^2);%Horizontal Distance
if hdist(i,1) < 1.e-20
   Az_rad(i,1) = 0; % Radians
   El_rad(i,1) = pi/2;% Radians
else
    
 %***COMPUTE AZIMUTH(Az);
try
   Az_rad(i,1) = atan2(E,N);% Radians
catch
   Az_rad(i,1)=rad2deg(atan2(E/norm(vENU),N/norm(vENU)));% Radians
end
if Az_rad(i,1) < 0
   Az_rad(i,1) = Az_rad(i,1)+ 2*pi;% Radians
end
%****COMPUTE ELEVATION ANGLE(El)
try
   El_rad(i,1) = atan2(U,hdist(i,1));% Radians
catch
   El_rad(i,1)=asin(U/norm(vENU));% Radians
end

end

catch
%REMARK:
%       This an alternative way to compute satellite Azimuth and Elevation
%       Angle.Both(the try and catch) are suppose to produce same results.

%**COMPUTE UNIT VECTOR FROM OBSERVATION STATION(USER) TO SATELLITE POSITION*
r=sqrt((dxyz(1)).^2+(dxyz(2)).^2+(dxyz(3)).^2);
Dxyz=dxyz./r;
dx = Dxyz(1);
dy = Dxyz(2);
dz = Dxyz(3);

%****COMPUTE ROTATED UNIT VECTORS IN ENU FORM OBSERVATION STATION TO ...   *
%                                                        SATELLITE POSITION*
try
%*****COMPUTING Rotation matrix
%REMARKS:
%       R is direction cosine matrix to transform position vector components
%       from geocentric equatorial frame into the topocentric horizon frame.    

R = [ -sin(long)                        cos(long)         0     ;
      -sin(lat)*cos(long)  -sin(lat)*sin(long)  cos(lat);
      cos(lat)*cos(long)   cos(lat)*sin(long)  sin(lat)];

%**COMPUTE UNIT VECTORS
vENU = R*Dxyz;

%***EXTRACTING INDIVIDUAL LOCAL COORDINATES
E = vENU(1);% 'East'-coordinate relative to local origin (meters)
N = vENU(2);% 'North'-coordinate relative to local origin (meters)
U = vENU(3);% Up-coordinate relative to local origin (meters)
   
catch
    
N =  sin(lat).*(-dx.*cos(long) - dy.*sin(long))+dz.*cos(lat);
E = - dx.* sin(long)+ dy.*cos(long) ;
U = cos(lat).*(dx.*cos(long)+dy.*sin(long)) + dz.*sin(lat);
               
end

%***COMPUTING AZIMUTH(Az) AND ELEVATION ANGLE(El)FROM NEU
        
%****COMPUTE ELEVATION ANGLE(El)
%ElRAD=asin(U/norm(vENU));% Radians
try
	El_rad(i,1)= (pi / 2 - acos(U));% Radians
    
catch
   LE = sqrt(E.^2+N.^2);
   El_rad(i,1) =(pi/2) - atan2(LE,U);% Radians
end
   El_deg(i,1)=rad2deg(El_rad(i,1));%Degrees
    
%***COMPUTE AZIMUTH(Az);
try
Az_rad(i,1) = atan2(E,N);% Radians
catch
Az_rad(i,1)=atan2(E/norm(vENU),N/norm(vENU));% Radians
end
%***CHECK FOR NEGATIVE ANGLES
if (Az_rad(i,1) < 0);
Az_rad(i,1) = Az_rad(i,1) + 2 * pi;
end
Az_deg(i,1) = rad2deg(Az_rad(i,1));	% degrees
end

%***COMPUTE RANGE
range(i,1) = sqrt(dxyz(1)^2+dxyz(2)^2+dxyz(3)^2);

%***CONVERT ANGLES IN RADIANS TO DEGREES
Az_deg(i,1) = rad2deg(Az_rad(i,1));%Degrees
El_deg(i,1) = rad2deg(El_rad(i,1));%Degrees
    
end

%***ARRANGE DATA FOR EACH USER POSITION/STATION
satAZ_rad(:,j)=Az_rad;  %Satellite Azimuth in radian
satEL_rad(:,j)=El_rad; %Satellite Elevation in radian
satAZ_deg(:,j)=Az_deg ; %Satellite Azimuth in decimal degrees
satEL_deg(:,j)=El_deg; %Satellite Elevation in decimal degrees
Range(:,j)=range;% Range(Distance)b/n User Position & Satellite

end %for j=1:size(UserXYZ,2)

%%%%%%%%%%%%%%%%%%%%%%%%END OF satAzimuthElevation.m %%%%%%%%%%%%%%%%%%%%%