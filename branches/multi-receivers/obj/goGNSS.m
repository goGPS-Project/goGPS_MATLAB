%   OBJECT goGNSS
% =========================================================================
%
% DESCRIPTION
%   Object to manage a useful functions / standard parameters of goGPS
%   Under development no methods are used up to now...
%
% EXAMPLE
%   go = goGNSS();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Stefano Caldera, Andrea Gatti, Lisa Pertusini
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------------------

classdef goGNSS < handle
    
    % Constant values
    % => to discriminate them from function (in autocompletition) they are
    % written in capital letters
    properties (Constant)
        % GENERIC CONSTANTS -----------------------------------------------
        
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
        
        ELL_A = 6378137;                      % Ellipsoid semi-major axis [m]
        ELL_F = 1/298.257222101;              % Ellipsoid flattening
        ELL_E = sqrt(1-(1-goGNSS.ELL_F)^2);   % Eccentricity
        
        OMEGAE_DOT = 7.2921151467e-5;         % Angular velocity of the Earth rotation [rad/s]
        GM = 3.986004418e14;                  % Gravitational constant * (mass of Earth) [m^3/s^2]
        
        MAX_SAT      = 32                     % Maximum number of active satellites in a constellation
        
        % CONSTELLATION SPECIFIC ------------------------------------------
        
        FG = [1575.420 1227.600]*1e6;         % GPS carriers frequencies [Hz]
        LAMBDAG = goGNSS.V_LIGHT ./ goGNSS.FG;% GPS carriers wavelengths [m]
        
        % CONSTELLATIONS IDs ----------------------------------------------
        
        ID_GPS     = 1 % Id of GPS constellation for goGPS internal use
        ID_GLONASS = 2 % Id of GLONASS constellation for goGPS internal use
        ID_GALILEO = 3 % Id of Galileo constellation for goGPS internal use
        ID_BEIDOU  = 4 % Id of BeiDou constellation for goGPS internal use
        ID_QZSS    = 5 % Id of QZSS constellation for goGPS internal use
        ID_SBAS    = 6 % Id of SBAS constellation for goGPS internal use
        
        % goGPS MODES -----------------------------------------------------
        
        MODE_RT_NAV          = 24;  % Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)
        MODE_RT_R_MON        = 21;  % Real Time Rover Monitor
        MODE_RT_M_MON        = 22;  % Real Time Master Monitor
        MODE_RT_RM_MON       = 23;  % Real Time Master + Rover Monitor
        
        MODE_PP_LS_C_SA      = 1;   % Post Proc Least Squares on Code Stand Alone
        MODE_PP_LS_CP_SA     = 3;   % Post Proc Least Squares on Code Stand Alone (BASE FOR FUTURE PPP IMPLEMENTATION)
        MODE_PP_LS_C_DD      = 11;  % Post Proc Least Squares on Code Double Differences
        MODE_PP_LS_CP_DD_L   = 13;  % Post Proc Least Squares on Code Double Differences with Lambda
        MODE_PP_LS_CP_VEL    = 3.1; % Post Proc Least Squares on Code and Phase for Velocity estimation
        
        MODE_PP_KF_C_SA      = 2;   % Post Proc Kalman Filter on Code Stand Alone
        MODE_PP_KF_C_DD      = 12;  % Post Proc Kalman Filter on Code Double Differencies
        MODE_PP_KF_CP_SA     = 4;   % Post Proc Kalman Filter on Code and Phase Stand Alone
        MODE_PP_KF_CP_DD     = 14;  % Post Proc Kalman Filter on Code and Phase Double Differences
        MODE_PP_KF_CP_DD_MR  = 15;  % Post Proc Kalman Filter on Code and Phase Double Differences Multiple Receivers
        
        GMODE_PP = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of post processing modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_C_DD ...
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_C_DD ...
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
        GMODE_RT = [ goGNSS.MODE_RT_NAV ...         % Group of real time modes
            goGNSS.MODE_RT_R_MON ...
            goGNSS.MODE_RT_M_MON ...
            goGNSS.MODE_RT_RM_MON];
        
        GMODE_SA = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of stand alone modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_CP_SA];
        
        GMODE_DD = [ goGNSS.MODE_PP_LS_C_DD ...     % Group of double differences modes
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_KF_C_DD ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
    end
    
    % Creator (empty)
    methods
        % Object to manage a useful functions / standard parameters of goGPS
        function obj = goGNSS()
        end
    end

    %   COMPATIBILITY FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to keep compatibility with the past
    methods (Static, Access = 'public')
        function f1 = F1()
        % GPS carriers frequencies F1
            f1 = goGNSS.FG(1);
        end        
        function f1 = F2()
        % GPS carriers frequencies F2
            f1 = goGNSS.FG(2);
        end
        
        function lambda1 = LAMBDA1()
        % GPS carriers frequency 1 [Hz]
            lambda1 = goGNSS.LAMBDAG(1);
        end
        
        function lambda2 = LAMBDA2()
        % GPS carriers frequency 2 [Hz]
            lambda2 = goGNSS.LAMBDAG(2);        

        end        
    end
    
    %   MODE FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to detect a certain kind of processing
    methods (Static, Access = 'public') 
        function isPostProcessing = isPP(mode)
        % return wheather or not the mode given in use is a Post Processing mode
            isPostProcessing = sum(intersect(mode, goGNSS.GMODE_PP));
        end
        
        function isRealTime = isRT(mode)
        % return wheather or not the mode given in use is a Real Time mode
            isRealTime = sum(intersect(mode, goGNSS.GMODE_RT));
        end
        
        function isDoubleDifferences = isDD(mode)
        % return wheather or not the mode given in use is a Double Difference mode
            isDoubleDifferences = sum(intersect(mode, goGNSS.GMODE_DD));
        end
        
        function isStandAlone = isSA(mode)
        % return wheather or not the mode given in use is a Stand Alone mode
            isStandAlone = sum(intersect(mode, goGNSS.GMODE_SA));
        end
    end
    
    %   USEFUL FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    methods (Static, Access = 'public')
        
        function [XR, dtR] = getBancroftPos(XS, dtS, prR)
            % get a bancroft solution from one receiver

            matB = [XS(:,:), prR(:) + goGNSS.V_LIGHT * dtS(:)]; % Bancroft matrix
            b = goGNSS.bancroft(matB);
            XR = b(1:3);
            dtR = b(4)/goGNSS.V_LIGHT;
        end
    
        
        function chiamalacomevuoi(goObs, goIni)

            ID_GNSS=1; % <- must be taken from the object!
            
            nRec=goObs.getNumRec();
            
            time_GPS=goObs.getTime_Ref();
            
            nFreq=goObs.getGNSSnFreq(ID_GNSS);
            Eph=goObs.getGNSSeph(ID_GNSS);
            iono=goObs.getIono();
            SP3_time=goObs.getGNSS_SP3time();
            SP3_coor=goObs.getGNSS_SP3coordinates();
            SP3_clck=goObs.getGNSS_SP3clock();
            
            % master receiver: preprocessing
            % ------------------------------
            time_M=goObs.getTime_M();
            [pos_M flag_M]= goObs.getPos_M(0);
            
            pr1_M=goObs.getGNSSpr_M(ID_GNSS, 0, 0, 1);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
            ph1_M=goObs.getGNSSph_M(ID_GNSS, 0, 0, 1);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
            snr_M=goObs.getGNSSsnr_M(ID_GNSS, 0, 0, 1);  %snr = getGNSSsnr_M(obj, idGNSS, idSat, idObs, nFreq)
            dop1_M=goObs.getGNSSdop_M(ID_GNSS, 0, 0, 1); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
            
            
            pr2_M=zeros(size(pr1_M));
            ph2_M=zeros(size(pr1_M));
            dop2_M=zeros(size(pr1_M));
            if nFreq==2
                pr2_M=goObs.getGNSSpr_M(ID_GNSS, 0, 0, 2);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
                ph2_M=goObs.getGNSSph_M(ID_GNSS, 0, 0, 2);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
                dop2_M=goObs.getGNSSdop_M(ID_GNSS, 0, 0, 2); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
            end
            
            fprintf('Master station ');
            [pr1_M, ph1_M, pr2_M, ph2_M, dtM, dtMdot] = pre_processing_clock(time_GPS, time_M, pos_M(:,1), pr1_M, ph1_M, ...
                pr2_M, ph2_M, snr_M, dop1_M, dop2_M, Eph, SP3_time, SP3_coor, SP3_clck, iono);
            fprintf('\n');
            % ------------------------------
            
            time_R=zeros(length(time_M),1,nRec);
            pr1_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            ph1_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            snr_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            dop1_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            pr2_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            ph2_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            dop2_R=zeros(goGNSS.MAX_SAT,length(time_M),nRec);
            
            
            dtR=NaN(length(time_M),1,nRec);
            dtRdot=NaN(length(time_M),1,nRec);
            
            
            for i=1:nRec
                time_R(:,1,i)=goObs.getTime_R(i); %time = getTime_R(obj, idRec)
                pos_R =[];
                pr1_R(:,:,i)=goObs.getGNSSpr_R(ID_GNSS,0,i,0,1);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                ph1_R(:,:,i)=goObs.getGNSSph_R(ID_GNSS,0,i,0,1);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                snr_R(:,:,i)=goObs.getGNSSsnr_R(ID_GNSS,0,i,0,1);      %snr = getGNSSsnr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                dop1_R(:,:,i)=goObs.getGNSSdop_R(ID_GNSS, 0, i, 0, 1); %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                
                if nFreq==2
                    pr2_R(:,:,i)=goObs.getGNSSpr_R(ID_GNSS,0,i,0,2);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    ph2_R(:,:,i)=goObs.getGNSSph_R(ID_GNSS,0,i,0,2);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    dop2_R(:,:,i)=goObs.getGNSSdop_R(ID_GNSS, 0, i, 0, 2); %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                end
                
                fprintf('Rover #%d ',i);
                [pr1_R(:,:,i), ph1_R(:,:,i), pr2_R(:,:,i), ph2_R(:,:,i), dtR(:,:,i), dtRdot(:,:,i)] = pre_processing_clock(time_GPS, time_R(:,1,i), pos_R, pr1_R(:,:,i), ph1_R(:,:,i), ...
                    pr2_R(:,:,i), ph2_R(:,:,i), snr_R(:,:,i), dop1_R(:,:,i), dop2_R(:,:,i), Eph, SP3_time, SP3_coor, SP3_clck, iono);
                fprintf('\n');
                
                %--> come modifico il contenuto globale di %%goObs.getGNSSpr_R(ID_GNSS,0,i,0,1) ????
                %--> perchè non mi tengo già le coordinate dei rover che escono da qui come valori a priori,
                %    invece di farle con bancroft ancora dopo?
                
                
            end
            
            
            
            % vanno tagliate le epoche all'inizio e alla fine, ci sono zeri!!!!!!
            index_epoch_common=find(time_R>0);
            first_epoch=index_epoch_common(1);
            
            
            
            
            %  QUESTO VA FATTO ADESSO? PENSO DI SI'
            %if (~flag_SP3)  <-- sistemare il flag_SP3 prendenolo dall'obj
                %remove satellites without ephemerides (GPS)
                delsat = setdiff(1:32,unique(Eph(1,:)));
                pr1_R(delsat,:,:) = 0;
                pr1_M(delsat,:,:) = 0;
                pr2_R(delsat,:,:) = 0;
                pr2_M(delsat,:,:) = 0;
                ph1_R(delsat,:,:) = 0;
                ph1_M(delsat,:,:) = 0;
                ph2_R(delsat,:,:) = 0;
                ph2_M(delsat,:,:) = 0;
                dop1_R(delsat,:,:) = 0;
                dop1_M(delsat,:,:) = 0;
                dop2_R(delsat,:,:) = 0;
                dop2_M(delsat,:,:) = 0;
                snr_R(delsat,:,:) = 0;
                snr_M(delsat,:,:) = 0;
            %end
            
            
            
            % Processing
            % ----------
            %  - for each rover receiver:
            %       - enhance coordinates with code and phase DD in single
            %         epoch (lambda)
            %  - estimation of apriori attitude
            %  - enhance solution with constrained least squares (DD in
            %         single epoch (lambda))
            
            
            
            
            % instrumental RS coordinates
            % ---------------------------
            % get geometry
            [geometry ev_point]=goIni.getGeometry();
            
            % barycenter definition
            %xb=mean(geometry,2);
            % barycentric instrumental RS coordinates
            %xR=geometry-repmat(xb,1,nRec);
            
            % non barycentric! the origin is the first receiver
            xR=geometry;
            
            
            % for each rover receiver: enhance coordinates with code and phase DD in single epoch (lambda)
            % --------------------------------------------------------------------------------------------
            
            % cosa sono???
            check_on = 0;
            check_off = 0;
            check_pivot = 0;
            check_cs = 0;
            plot_t = 1;
            %
            
            

            

            global Xhat_t_t  % forse conviene aggiungere output alla funzione goGPS_LS_DD_code_phase
            
            
            global azR azM elR elM distR distM;
            satCoord = struct('az',zeros(goGNSS.MAX_SAT,nRec+1),'el',zeros(goGNSS.MAX_SAT,nRec+1),'dist',zeros(goGNSS.MAX_SAT,nRec+1)); %first column: MASTER, further columns: ROVERS
            err_iono = NaN(goGNSS.MAX_SAT, nRec+1);
            err_tropo = NaN(goGNSS.MAX_SAT, nRec+1);
            
            
  
            
keyboard
            
            for t = first_epoch : first_epoch
                XR_DD=NaN(3,1,nRec);
                %cartesian to geodetic conversion of ROVER coordinates
                [phiM, lamM, hM] = cart2geod(pos_M(1,t), pos_M(2,t), pos_M(3,t));
                
                %radians to degrees
                phiM = phiM * 180 / pi;
                lamM = lamM * 180 / pi;
                
                
                %                 Eph_t = Eph(:,:,t);
                Eph_t=rt_find_eph(goObs.getGNSSeph(goGNSS.ID_GPS), goObs.getTime_Ref(t));
                for i=1:nRec
                    statistic = zeros(2,length(time_GPS)); % <-- VA DIMENSIONATO IN 3D!!!!?
                    ambiguity = 0;                         % <-- VA DIMENSIONATO IN 3D!!!!?
                    
                    
                    
                    Xhat_t_t=zeros(size(Xhat_t_t));
                    
                    [X_sat conf_sat(1:32,i)]=goGPS_LS_DD_code_phase(time_GPS(t), pos_M(:,t), pr1_R(:,t,i), pr1_M(:,t), pr2_R(:,t,i), pr2_M(:,t), ph1_R(:,t,i), ph1_M(:,t), ph2_R(:,t,i), ph2_M(:,t), snr_R(:,t,i), snr_M(:,t), Eph_t, SP3_time, SP3_coor, SP3_clck, iono, 1);

                    XR_DD(1,t,i)=Xhat_t_t(1);
                    XR_DD(2,t,i)=Xhat_t_t(3);
                    XR_DD(3,t,i)=Xhat_t_t(5);
                    
                    % compute elevation and atmospheric corrections from XR_DD coordinates
                    [phiR, lamR, hR] = cart2geod(XR_DD(1,t,i), XR_DD(2,t,i), XR_DD(3,t,i));
                    phiR = phiR * 180 / pi;
                    lamR = lamR * 180 / pi;
 
                    distM(~distM)=NaN;
                    azM(isnan(distM))=NaN;
                    elM(isnan(distM))=NaN;
                    azR(isnan(distM))=NaN;
                    elR(isnan(distM))=NaN;
                    distR(isnan(distM))=NaN;
                    
                    satCoord.az(:,1)= azM;
                    satCoord.el(:,1)= elM;
                    satCoord.dist(:,1)= distM;
                    
                    satCoord.az(:,i+1)= azR;
                    satCoord.el(:,i+1)= elR;
                    satCoord.dist(:,i+1)= distR;
                    
                    %computation of atmospheric errors of Master receiver
                    err_tropo(:,1) = tropo_error_correction(satCoord.el(:,1), hM);
                    err_iono(:,1) = iono_error_correction(phiM, lamM, satCoord.az(:,1), satCoord.el(:,1), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),1)=NaN;
                    
                    %computation of atmospheric errors of Rover receiver
                    err_tropo(:,i+1) = tropo_error_correction(satCoord.el(:,i+1), hR);
                    
                    %computation of ionospheric errors
                    err_iono(:,i+1) = iono_error_correction(phiR, lamR, satCoord.az(:,i+1), satCoord.el(:,i+1), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),i+1)=NaN;
                    
                end
                
            end
            
            % computing apriori attitude
            % --------------------------
            
            % global XYZ and geographic coordinates of barycenter
            %Xb_apriori=mean([XR_DD(:,t,1), XR_DD(:,t,2),XR_DD(:,t,3)],2);
            
            % the origin is the first point
            %------------------------------
            Xb_apriori=[XR_DD(1,t,1), XR_DD(2,t,1),XR_DD(3,t,1)]';
            [phi_b_apriori, lam_b_apriori, h_b_apriori] = cart2geod(Xb_apriori(1), Xb_apriori(2), Xb_apriori(3));
            %------------------------------
            
            % rotation matrix from local to global coordinates,
            % centered into the receiver barycenter
            R1t=[-sin(lam_b_apriori) cos(lam_b_apriori) 0; ...
                -sin(phi_b_apriori)*cos(lam_b_apriori) -sin(phi_b_apriori)*sin(lam_b_apriori) cos(phi_b_apriori); ...
                cos(phi_b_apriori)*cos(lam_b_apriori) cos(phi_b_apriori)*sin(lam_b_apriori) sin(phi_b_apriori)];
            
            
            % compute local East-North-Up coordinates from XYZ obtained
            % from CODE+PHASE DD
            XRl_apriori=NaN(3,nRec);
            for i=1:nRec
                XRl_apriori(1:3,i)=R1t*(XR_DD(1:3,t,i)-Xb_apriori);
            end
            % compute Euler parameters of the rotation from local to
            % instrumental RF
            % xR=Rt*XRl -> where Rt is a 3x3 rotation matrix containing
            % the Euler parameters r11, r21, r31 ; r12 ... (it's
            % transposed)
            
            % non funziona se cerco il baricentro, viene A linearmente
            % dipendente! per questo l'origine è sul primo punto
%             y0=xR(:);
%             A=NaN(nRec*3,9);
%             for i=1:nRec
%                 A((i-1)*3+1:i*3,1:9)=[XRl_apriori(1,i) XRl_apriori(2,i) XRl_apriori(3,i) 0 0 0 0 0 0 ; ...
%                     0 0 0 XRl_apriori(1,i) XRl_apriori(2,i) XRl_apriori(3,i) 0 0 0 ; ...
%                     0 0 0 0 0 0 XRl_apriori(1,i) XRl_apriori(2,i) XRl_apriori(3,i)];
%             end
%             euler_parameters=inv(A'*A)*A'*y0;
%                        
%             %Cxx_euler=(y0-A*euler_parameters)'*(y0-A*euler_parameters)/(size(A,1)-size(A,2)+1)*(inv(A'*A));
%             % non c'è ridondanza e il condizionamento di N è uno schifo
%             
%             roll_approx=mod(atan2(euler_parameters(8),euler_parameters(9)),2*pi);
%             yaw_approx=mod(atan2(-cos(roll_approx)*euler_parameters(2)+sin(roll_approx)*euler_parameters(3) , cos(roll_approx)*euler_parameters(5)-sin(roll_approx)*euler_parameters(6)),2*pi);
%             pitch_approx=mod(atan2(-euler_parameters(7),sin(roll_approx)*euler_parameters(8)+cos(roll_approx)*euler_parameters(9)),2*pi);
%             
            

            % block solution
            euler_parameters=[];
            A=NaN(nRec,3);
            for i=1:nRec
                A(i,1:3)=[XRl_apriori(1,i) XRl_apriori(2,i) XRl_apriori(3,i)];
            end
            y0=xR(:);
            euler_parameters(1:3,1)=inv(A'*A)*A'*y0(1:3:end);
            euler_parameters(4:6,1)=inv(A'*A)*A'*y0(2:3:end);            
            euler_parameters(7:9,1)=inv(A'*A)*A'*y0(3:3:end);  
            roll_approx=mod(atan2(euler_parameters(8),euler_parameters(9)),2*pi);
            yaw_approx=mod(atan2(-cos(roll_approx)*euler_parameters(2)+sin(roll_approx)*euler_parameters(3) , cos(roll_approx)*euler_parameters(5)-sin(roll_approx)*euler_parameters(6)),2*pi);
            pitch_approx=mod(atan2(-euler_parameters(7),sin(roll_approx)*euler_parameters(8)+cos(roll_approx)*euler_parameters(9)),2*pi);
                      
            attitude_approx=[roll_approx;pitch_approx;yaw_approx];
            
            
            % code + phase double differenecs with Xb and attitude
            XS=NaN(goGNSS.MAX_SAT,3);
            XS(conf_sat==1,1:3)=X_sat;
            sat=find(conf_sat==1);
            
            %actual pivot

            [null_max_elR pivot_index]=max(satCoord.el(sat,1));
            pivot = sat(pivot_index);
            

                        
            if (size(sat,1) >= 4) % & cond_num < cond_num_threshold)
                phase_1=1;
                               
                %loop is needed to improve the atmospheric error correction
                for i = 1 : 3
                    %if (phase == 1)
                        [Xb_apriori, N1_hat, cov_XR, cov_N1, PDOP, HDOP, VDOP, up_bound, lo_bound, posType, attitude_approx, XR_DD(:,t,:)] = LS_DD_code_phase_MR(Xb_apriori, XR_DD(:,t,:), pos_M(:,t), X_sat, pr1_R(sat,t,:), ph1_R(sat,t,:), snr_R(sat,t,:), pr1_M(sat,t), ph1_M(sat,t), snr_M(sat,t), satCoord.el(sat,2:nRec+1), satCoord.az(sat,1), err_tropo(sat,2:nRec+1), err_iono(sat,2:nRec+1), err_tropo(sat,1), err_iono(sat,1), pivot_index, phase_1, attitude_approx, geometry);
                    %else
                    %    [XR, N1_hat, cov_XR, cov_N1, PDOP, HDOP, VDOP, up_bound, lo_bound, posType] = LS_DD_code_phase(XR, XM, XS, pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(st), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase);
                    %end
                 
                    %[phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
                    %[azR(azR ~= 0), elR(elR ~= 0), distR(distR ~= 0)] = topocent(XR, XS);
                    
                    %err_tropo_R = tropo_error_correction(elR(elR ~= 0), hR);
                    %err_iono_R = iono_error_correction(phiR*180/pi, lamR*180/pi, azR(azR ~= 0), elR(elR ~= 0), time_rx, iono, []);
                end
                keyboard
            else
                if (~isempty(Xhat_t_t))
                    XR = Xhat_t_t([1,o1+1,o2+1]);
                    pivot = 0;
                else
                    return
                end
            end
            
            
            
            
            
        end
        
        

        
        
        
            
            
            
                           %% computation of linearized equations
                % (expression of XR as funtction of xR)
                
%                 syms s_Xb s_Yb s_Zb s_phi_b s_lam_b s_roll s_pitch s_yaw s_x s_y s_z s_XS s_YS s_ZS
%                 
%                 % rotation from local to global
%                 s_Rgl=[-sin(s_lam_b) cos(s_lam_b) 0; ...
%                     -sin(s_phi_b)*cos(s_lam_b) -sin(s_phi_b)*sin(s_lam_b) cos(s_phi_b); ...
%                     cos(s_phi_b)*cos(s_lam_b) cos(s_phi_b)*sin(s_lam_b) sin(s_phi_b)];
%                 
%                 % rotation from instrumental to local
%                 s_Rli=[cos(s_roll)*cos(s_pitch) -sin(s_pitch)*cos(s_yaw)+cos(s_roll)*sin(s_pitch)+sin(s_yaw) sin(s_roll)*sin(s_yaw)+cos(s_roll)*sin(s_pitch)*cos(s_yaw); ...
%                     sin(s_roll)+cos(s_pitch) sin(s_roll)*sin(s_pitch)*sin(s_yaw) sin(s_roll)*sin(s_pitch)*cos(s_yaw)-cos(s_roll)*sin(s_yaw); ...
%                     -sin(s_pitch) cos(s_pitch)*sin(s_yaw) cos(s_pitch)*cos(s_yaw)];
%                 
%                 
%                 s_X=[s_Xb; s_Yb; s_Zb]+s_Rgl*s_Rli*[s_x;s_y;s_z];                
%                 s_PR=sqrt((s_X(1)-s_XS)^2+(s_X(2)-s_YS)^2 + (s_X(3)-s_ZS)^2);       
%                 
%                 
%                 s_Ai=[diff(s_PR,s_Xb) diff(s_PR,s_Yb) diff(s_PR,s_Zb) diff(s_PR,s_roll) diff(s_PR,s_pitch) diff(s_PR,s_yaw)];
%                 
%                 F_Ai=inline(s_Ai); %% <-- colonne della matrice disegno corrispondenti alle incognite geometriche 
%                 F_PR=inline(s_PR); %% <-- parte geometrica dei termini noti
%                 
%                 
                
            
            
            
            
            
            
            
            
       
        
        
        
        
        
        
        function [X, dt, usableSat, satCoord, cov_X, var_dt, PDOP, HDOP, VDOP, cond_num] = LS_MR_C_SA(goObs, goIni)

            n_rec=goObs.getNumRec();
            
            
            
            %------------------------------------------------------------------------------------
            % APPROXIMATE POSITION for the rovers and the master
            %-----------------------------------------------------------------------------------
           
            [XM0 flag_M] = goObs.getX0_M();    %[sized 3x1]
            
            % get a Least Squares solution using only Code Stand Alone using Master and Multiple Rover Receivers
            
            %----------------------------------------------------------------------------------------------
            % FIRST ESTIMATE OF SATELLITE POSITIONS from Master observations
            %----------------------------------------------------------------------------------------------
            % satellites in view by the master at the first epoch
            sat_pr_init = find(goObs.getTrackedSat_pr(goGNSS.ID_GPS, 0, -1, 1, 1) ~= 0);
            % extract ephemerides for the first epoch
            Eph_1 = rt_find_eph (goObs.getGNSSeph(goGNSS.ID_GPS), goObs.getTime_Ref(1));
            XS = NaN(goGNSS.MAX_SAT,3);
            dtS = NaN(goGNSS.MAX_SAT,1);
            XS_tx = NaN(goGNSS.MAX_SAT,3);
            VS_tx = NaN(goGNSS.MAX_SAT,3);
            time_tx = NaN(goGNSS.MAX_SAT,1);
            no_eph = NaN(goGNSS.MAX_SAT,1);
            [XS_cut, dtS_cut, XS_tx_cut, VS_tx_cut, time_tx_cut, no_eph_cut] = satellite_positions(goObs.getTime_M(1), goObs.getGNSSpr_M(goGNSS.ID_GPS, 0, 1, 1), sat_pr_init, Eph_1, goObs.getGNSS_SP3time(), goObs.getGNSS_SP3coordinates(), goObs.getGNSS_SP3clock(), goObs.getSBAS(), zeros(length(sat_pr_init),1), goObs.getIono(), 0);
            XS(sat_pr_init,:) = XS_cut;
            dtS(sat_pr_init) = dtS_cut;
            XS_tx(sat_pr_init,:) = XS_tx_cut;
            VS_tx(sat_pr_init,:) = VS_tx_cut;
            time_tx(sat_pr_init) = time_tx_cut;
            no_eph(sat_pr_init) = no_eph_cut;
            %satellites with no ephemeris available
            usableSatM = (no_eph == 0);
            
            %----------------------------------------------------------------------------------------------
            % APPROXIMATE RECEIVER POSITION BY BANCROFT ALGORITHM
            %----------------------------------------------------------------------------------------------
            XR_approx = zeros(3,n_rec);
            dtR = zeros(1,n_rec);
            prR = reshape(goObs.getGNSSpr_R(goGNSS.ID_GPS, 0, 0, 1, 1),goGNSS.MAX_SAT,n_rec);
            for r=1:n_rec
                [XR_approx(:,r), dtR(:,r)] = goGNSS.getBancroftPos(XS(usableSatM,:), dtS(usableSatM), prR(usableSatM,r));
            end
            
            %----------------------------------------------------------------------------------------------
            % ELEVATION CUTOFF, SNR CUTOFF AND REMOVAL OF SATELLITES WITHOUT EPHEMERIS
            %----------------------------------------------------------------------------------------------
            snr = NaN(goGNSS.MAX_SAT,n_rec+1);
            snr(sat_pr_init,1) = goObs.getGNSSsnr_M(goGNSS.ID_GPS, sat_pr_init, 1, 1);
            for r=1:n_rec
                snr(sat_pr_init,r+1) = goObs.getGNSSsnr_R(goGNSS.ID_GPS, sat_pr_init, r, 1, 1);
            end
            satCoord = struct('az',zeros(goGNSS.MAX_SAT,n_rec+1),'el',zeros(goGNSS.MAX_SAT,n_rec+1),'dist',zeros(goGNSS.MAX_SAT,n_rec+1)); %first column: MASTER, further columns: ROVERS
            %master
            %topacentric satellite coordinates (azimuth, elevation, distance)
            [satCoord.az(:,1), satCoord.el(:,1), satCoord.dist(:,1)] = topocent(XM0, XS);
            %elevation cutoff, SNR cutoff and removal of satellites without ephemeris
            if (any(snr(:,1)))
                usableSatM = (satCoord.el(:,1) > goIni.getCutoff) & isfinite(snr(:,1)) & (snr(:,1) > goIni.getSnrThr) & (usableSatM);
            else
                usableSatM = (satCoord.el(:,1) > goIni.getCutoff) & (usableSatM);
            end
            %rover
            usableSatR = NaN(goGNSS.MAX_SAT, n_rec);
            for r=1:n_rec
                %topacentric satellite coordinates (azimuth, elevation, distance)
                [satCoord.az(:,r+1), satCoord.el(:,r+1), satCoord.dist(:,r+1)] = topocent(XR_approx(:,r), XS);
                %elevation cutoff, SNR cutoff and removal of satellites without ephemeris
                if (any(snr(:,r+1)))
                    usableSatR(:,r) = (satCoord.el(:,r+1) > goIni.getCutoff) & isfinite(snr(:,r+1)) & (snr(:,r+1) > goIni.getSnrThr) & (usableSatM);
                else
                    usableSatR(:,r) = (satCoord.el(:,r+1) > goIni.getCutoff) & (usableSatM);
                end
            end
            
            usableSat = [usableSatM usableSatR];
            %--------------------------------------------------------------------------------------------
            % LEAST SQUARES SOLUTION
            %--------------------------------------------------------------------------------------------
            % if I have at least 4 satellites for each receiver (master included) 
            if sum(sum(usableSat)>=4) == n_rec+1
                err_tropo = NaN(goGNSS.MAX_SAT, n_rec);
                err_iono = NaN(goGNSS.MAX_SAT, n_rec);
                %master
                %cartesian to geodetic conversion of ROVER coordinates
                [phiM, lamM, hM] = cart2geod(XM0(1), XM0(2), XM0(3));
                
                %radians to degrees
                phiM = phiM * 180 / pi;
                lamM = lamM * 180 / pi;
                
                %computation of tropospheric errors
                err_tropo(:,1) = tropo_error_correction(satCoord.el(:,1), hM);
                
                %computation of ionospheric errors
                err_iono(:,1) = iono_error_correction(phiM, lamM, satCoord.az(:,1), satCoord.el(:,1), goObs.getTime_Ref(1), goObs.getIono(), goObs.getSBAS());
                %rovers
                for r=1:n_rec
                    %cartesian to geodetic conversion of ROVER coordinates
                    [phiR(:,r), lamR(:,r), hR(:,r)] = cart2geod(XR_approx(1,r), XR_approx(2,r), XR_approx(3,r));
                    
                    %radians to degrees
                    phiR = phiR * 180 / pi;
                    lamR = lamR * 180 / pi;
                    
                    %computation of tropospheric errors
                    err_tropo(:,r+1) = tropo_error_correction(satCoord.el(:,r+1), hR(:,r));
                    
                    %computation of ionospheric errors
                    err_iono(:,r+1) = iono_error_correction(phiR(:,r), lamR(:,r), satCoord.az(:,r+1), satCoord.el(:,r+1), goObs.getTime_Ref(1), goObs.getIono(), goObs.getSBAS());
                end

                % inizialize dtR
                dtR_approx=zeros(1,n_rec+1);
                
                
                % get geometry
                [geometry ev_point]=goIni.getGeometry();

                
                %% obtaining improved receiver coordinates using constraints on distances with least squares

                % compute distances from object coordinates
                
                keyboard 
                % graph     
                baselines=[];
                
                %distances
                distances_3D=[];
                distances_2D=[];

                for i=1:n_rec-1
                    for j=i+1:n_rec
                        baselines=[baselines;i,j];                        
                        distances_3D=[distances_3D;goGNSS.compute_distance(geometry(:,i),geometry(:,j))];  
                        distances_2D=[distances_2D;goGNSS.compute_distance(geometry(1:2,i),geometry(1:2,j))];  
                    end
                end

                % pre-allocation                
                
                % y0: observation vector
                y0 = NaN(sum(sum(usableSat(:,2:end)))+length(distances_3D),1);
                
                %known term vector
                b = NaN(sum(sum(usableSat(:,2:end)))+length(distances_3D),1);
                
                %observation covariance matrix
                Q = zeros(sum(sum(usableSat(:,2:end)))+length(distances_3D));
                
                % A: design matrix..
                A = zeros(sum(sum(usableSat(:,2:end)))+length(distances_3D),4*n_rec);
                n_row_A=0;
                
                
                % pseuduoranges of all the receivers from all the visible satellites and known distances between receivers
                for i=1:n_rec
                    n_sat_i= sum(usableSat(:,i+1));
                    index_sat_i=find(usableSat(:,i+1)==1);
                    
                    %satellite-receivers geometry
                    A(n_row_A+1:n_row_A+n_sat_i,(i-1)*3+1:i*3)=[(XR_approx(1,i) - XS(index_sat_i,1)) ./ satCoord.dist(index_sat_i,i+1), ... %column for X coordinate
                        (XR_approx(2,i) - XS(index_sat_i,2)) ./ satCoord.dist(index_sat_i,i+1), ... %column for Y coordinate
                        (XR_approx(3,i) - XS(index_sat_i,3)) ./ satCoord.dist(index_sat_i,i+1)];
                    
                    %receiver clocks
                    A(n_row_A+1:n_row_A+n_sat_i,3*n_rec+i)=ones(n_sat_i,1);
                    
                    y0(n_row_A+1:n_row_A+n_sat_i,1) = prR(index_sat_i,i);
                    b(n_row_A+1:n_row_A+n_sat_i,1) = satCoord.dist(index_sat_i,i+1) - goGNSS.V_LIGHT.*dtS(index_sat_i) + err_tropo(index_sat_i,i+1) + err_iono(index_sat_i,i+1) + goGNSS.V_LIGHT.*dtR_approx(i+1);
                    
                    Q(n_row_A+1:n_row_A+n_sat_i,n_row_A+1:n_row_A+n_sat_i)= cofactor_matrix_SA(satCoord.el(index_sat_i,i+1), snr(index_sat_i,i+1));
                    
                    n_row_A=n_row_A+n_sat_i;
                end
                

                % constraints on distances
                for i=1:size(distances_3D,1);
                    % indexes of the two receivers
                    rec_i = baselines(i,1);
                    rec_j = baselines(i,2);
                    dist_ij_approx=goGNSS.compute_distance(XR_approx(:,rec_i),XR_approx(:,rec_j));
                    

                    A(n_row_A+1,(rec_i-1)*3+1:rec_i*3) = [(XR_approx(1,rec_i)-XR_approx(1,rec_j))/dist_ij_approx ...
                        +(XR_approx(2,rec_i)-XR_approx(2,rec_j))/dist_ij_approx ...
                        +(XR_approx(3,rec_i)-XR_approx(3,rec_j))/dist_ij_approx];
                    A(n_row_A+1,(rec_j-1)*3+1:rec_j*3) = -[(XR_approx(1,rec_i)-XR_approx(1,rec_j))/dist_ij_approx ...
                        +(XR_approx(2,rec_i)-XR_approx(2,rec_j))/dist_ij_approx ...
                        +(XR_approx(3,rec_i)-XR_approx(3,rec_j))/dist_ij_approx];
                    
                    y0(n_row_A+1,1)= distances_3D(i);
                    b(n_row_A+1,1) = dist_ij_approx;  
                    
                    Q(n_row_A+1,n_row_A+1) = 10^-4;
                    n_row_A = n_row_A+1;
                end


                %Least-Squares adjustement
                %normal matrix
                N = (A'*(Q^-1)*A);
                
                %least squares solution
                x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
                XR_hat  = XR_approx + reshape(x_hat(1:n_rec*3),3,n_rec);
                dtR_hat(2:n_rec+1) = (dtR_approx(2:end)+[x_hat(3*n_rec+1:end)]') ./ goGNSS.V_LIGHT; %#ok<NBRAK>



                % number of observations (n) and unknown (m)
                [n,m] = size(A);
                
                %estimation of the variance of the observation error
                y_hat = A*x_hat + b;
                v_hat = y0 - y_hat;
                sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);
                

                
                % try with OOLO
                %[x_oolo Cxx_oolo y_hat_oolo Cyy_hat_oolo v_hat_oolo Cuu_hat_oolo sigma02_hat_oolo index_outlier_oolo v_hat_all_oolo]=goGNSS.OLOO_General(A, y0-b, Q,ones(n,1), 0.95,0);
                
                
                %% compute apriori values for rotations 
                
                % global XYZ and geographic coordinates of barycenter
                Xb_hat=mean(XR_hat,2);
                [phi_b_hat, lam_b_hat, h_b_hat] = cart2geod(Xb_hat(1), Xb_hat(2), Xb_hat(3));
                
                % rotation matrix from local to global coordinates,
                % centered into the receiver barycenter
                R1t=[-sin(lam_b_hat) cos(lam_b_hat) 0; ...
                    -sin(phi_b_hat)*cos(lam_b_hat) -sin(phi_b_hat)*sin(lam_b_hat) cos(phi_b_hat); ...
                    cos(phi_b_hat)*cos(lam_b_hat) cos(phi_b_hat)*sin(lam_b_hat) sin(phi_b_hat)];
                
                
                % compute local East-North-Up coordinates from XYZ obtained
                % from the constrained LS
                XRl_hat=NaN(3,n_rec);
                for i=1:n_rec
                    XRl_hat(1:3,i)=R1t*(XR_hat(1:3,i)-Xb_hat);
                end
                
                % instrumental RS coordinates
                % barycenter definition
                xb=mean(geometry,2);
                % barycentric instrumental RS coordinates  
                xR=geometry-repmat(xb,1,n_rec);
                
                % compute Euler parameters of the rotation from local to
                % instrumental RF
                % xR=Rt*XRl -> where Rt is a 3x3 rotation matrix containing
                % the Euler parameters r11, r21, r31 ; r12 ... (it's
                % transposed)
                y0=xR(:);
                A=NaN(n_rec*3,9);
                for i=1:n_rec
                    A((i-1)*3+1:i*3,1:9)=[XRl_hat(1,i) XRl_hat(2,i) XRl_hat(3,i) 0 0 0 0 0 0 ; ...
                        0 0 0 XRl_hat(1,i) XRl_hat(2,i) XRl_hat(3,i) 0 0 0 ; ...
                        0 0 0 0 0 0 XRl_hat(1,i) XRl_hat(2,i) XRl_hat(3,i)];
                end
                
                euler_parameters=inv(A'*A)*A'*y0;                
                roll_approx=mod(atan2(euler_parameters(8),euler_parameters(9)),2*pi);
                yaw_approx=mod(atan2(-cos(roll)*euler_parameters(2)+sin(roll)*euler_parameters(3) , cos(roll)*euler_parameters(5)-sin(roll)*euler_parameters(6)),2*pi);
                pitch_approx=mod(atan2(-euler_parameters(7),sin(roll)*euler_parameters(8)+cos(roll)*euler_parameters(9)),2*pi);
                
                % force the apriori attitude to be planar
                roll_approx=0;
                pitch_approx=0;    
                

                %% computation of linearized equations
                % (expression of XR as funtction of xR)
                
                syms s_Xb s_Yb s_Zb s_phi_b s_lam_b s_roll s_pitch s_yaw s_x s_y s_z s_XS s_YS s_ZS
                
                % rotation from local to global
                s_Rgl=[-sin(s_lam_b) cos(s_lam_b) 0; ...
                    -sin(s_phi_b)*cos(s_lam_b) -sin(s_phi_b)*sin(s_lam_b) cos(s_phi_b); ...
                    cos(s_phi_b)*cos(s_lam_b) cos(s_phi_b)*sin(s_lam_b) sin(s_phi_b)];
                
                % rotation from instrumental to local
                s_Rli=[cos(s_roll)*cos(s_pitch) -sin(s_pitch)*cos(s_yaw)+cos(s_roll)*sin(s_pitch)+sin(s_yaw) sin(s_roll)*sin(s_yaw)+cos(s_roll)*sin(s_pitch)*cos(s_yaw); ...
                    sin(s_roll)+cos(s_pitch) sin(s_roll)*sin(s_pitch)*sin(s_yaw) sin(s_roll)*sin(s_pitch)*cos(s_yaw)-cos(s_roll)*sin(s_yaw); ...
                    -sin(s_pitch) cos(s_pitch)*sin(s_yaw) cos(s_pitch)*cos(s_yaw)];
                
                
                s_X=[s_Xb; s_Yb; s_Zb]+s_Rgl*s_Rli*[s_x;s_y;s_z];                
                s_PR=sqrt((s_X(1)-s_XS)^2+(s_X(2)-s_YS)^2 + (s_X(3)-s_ZS)^2);       
                
                
                s_Ai=[diff(s_PR,s_Xb) diff(s_PR,s_Yb) diff(s_PR,s_Zb) diff(s_PR,s_roll) diff(s_PR,s_pitch) diff(s_PR,s_yaw)];
                
                F_Ai=inline(s_Ai); %% <-- colonne della matrice disegno corrispondenti alle incognite geometriche 
                F_PR=inline(s_PR); %% <-- parte geometrica dei termini noti
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
figure
hold on
axis equal
plot3([xR(1,:) xR(1,1)], [xR(2,:) xR(2,1)], [xR(3,:) xR(3,1)],'-r');
plot3([XRl(1,:) XRl(1,1)], [XRl(2,:) XRl(2,1)] , [XRl(3,:) XRl(3,1)],'-b');
grid on

% distanze a posteriori
distances_3D_post=[];
distances_2D_post=[];
for i=1:n_rec-1
    for j=i+1:n_rec
        distances_3D_post=[distances_3D_post;goGNSS.compute_distance(XR_hat(:,i),XR_hat(:,j))];
        distances_2D_post=[distances_2D_post;goGNSS.compute_distance(XR_hat(1:2,i),XR_hat(1:2,j))];
    end
end

              
                %%LS solution for ROVER receivers
                %[XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code_nRec(XR, XS, prR, snr(:,2:end), satCoord.el(:,2:end), satCoord.dist(:,2:end), dtS, err_tropo(:,2:end), err_iono(:,2:end), usableSat(:,2:end), geometry);
                

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                for r=1:n_rec
                    %satellite topocentric coordinates (azimuth, elevation, distance)
                    [satCoord.az(:,r), satCoord.el(:,r), satCoord.dist(:,r)] = topocent(XR(:,r), XS);
                end
            else
                %empty variables
                dtR = [];
                err_tropo = [];
                err_iono  = [];
                
                cov_XR = [];
                var_dtR = [];
                
                PDOP = -9999;
                HDOP = -9999;
                VDOP = -9999;
                cond_num = [];
            end
            
        end
        
    end
    
    methods (Static, Access = 'private')
        % function to compute distance between points of known coordinates
        % a and b should be vectors with 2 or 3 components
        function dist=compute_distance(a,b)
            dist= sqrt(sum((a-b).^2));       
        end
        
        
        % optimized Leave One Out for outlier detection in LS approach
        function [xfin Cxx yfin Cyy ufin Cuu s2fin imax um]=OLOO_General(A, y, Q, Dim, significance,rms_threshold)
            %Purpose:   perform LS on blocks of correlated observations
            %           identify one (block) outlier
            %           reject it
            %           re-estimate unknowns
            %           according to the theory in "##"
            %1.0: Ludovico Biagi, Stefano Caldera, 06.02.2012
            %input:
            %   A: design matrix
            %   y: observations vector
            %   Q: cofactor matrix
            %   Dim: dimension of each block of correlated observations,
            %   significance: test significance
            %output
            %   xfin: final parameters estimates
            %   Cxx: covariance matrix
            %   yfin: final observation estimates (outliers rejected)
            %   Cyy: covariance matrix
            %   ufin: estimated residuals (outlier rejected)
            %   Cuu= covariance matrix
            %   s2fin: final estimated variance
            %   imax: index of the rejected block
            %   um: vector of residuals of the intial global solution

            
            m=length(Dim);
            
            
            %fprintf(' ******** OOLO v.1.0 *********\n');
            %fprintf(' *  Optimized Leave One Out  * \n');
            %fprintf(' ***************************** \n\n');
            
            [M,n]=size(A);
            
            Qm=cell(m,1);
            Am=cell(m,1);
            Bm=cell(m,1);
            Cm=cell(m,1);
            Km=cell(m,1);
            Kminv=cell(m,1);
            Wm=cell(m,1);
            ym=cell(m,1);
            um=cell(m,1);
            N=zeros(n,n);
            xsi=0;
            s2cap=0;
            wm=cell(m,1);
            s2m=zeros(m,1);
            Qw=cell(m,1);
            Qwinv=cell(m,1);
            F=zeros(m,1);
            Flim=zeros(m,1);
            
            Row(1)=1;
            for i=1:m
                if i>1
                    Row(i)=Row(i-1)+Dim(i-1);
                end
                Qm{i}=Q(Row(i):Row(i)+Dim(i)-1,Row(i):Row(i)+Dim(i)-1);
                Am{i}=A(Row(i):Row(i)+Dim(i)-1,:);
                Wm{i}=inv(Qm{i});
                ym{i}=y(Row(i):Row(i)+Dim(i)-1);
            end
            
            
            %% compute the global solution
            for i=1:m
                N = N+Am{i}'*Wm{i}*Am{i};
                xsi =xsi+Am{i}'*Wm{i}*ym{i};
            end
            Ninv=N\eye(length(N));
            xcap=Ninv*xsi;
            
            for i=1:m
                um{i}=ym{i}-Am{i}*xcap;
                s2cap=s2cap+um{i}'*Wm{i}*um{i};
            end
            s2cap =s2cap/(M-n);
            
            
            %fprintf('       -> S0  : %.3e\n',s2cap^.5);
            %fprintf('       -> Chi2: %.2f, Chi2Lim: %.2f \n',chi2cap,chi2lim);
            %fprintf('       -> MaxU: %.3e, in Block: %d \n',u0_max,u0_max_ib);
            %fprintf('       -> MaxV: %.3e, in Block: %d \n',v0_max,v0_max_ib);
            
            stop=0;
            if M-n>2 && s2cap>rms_threshold^2
                
                %% start outliers rejection
                for i=1:m
                    Im=eye(n);
                    Bm{i}=Ninv*Am{i}';
                    Cm{i}=Am{i}*Bm{i};
                    Km{i}=Qm{i}-Cm{i};
                    if sum(sum(Km{i}== zeros(Dim(i))))==Dim(i)*Dim(i)
                        %fprintf('WARNING!! Cluster %d cannot be checked. It is essential to the estimation\n\n',i);
                        stop=1;
                        break
                    else
                        %Kminv{i}=inv(Km{i});
                        Kminv{i}=Km{i}\eye(length(Km{i}));
                        wm{i}=Qm{i}*Kminv{i}*um{i};
                        if M-n-Dim(i)<1
                            %fprintf('WARNING!! Cluster %d cannot be checked. Redudancy= %d\n\n',i, M-n-Dim(i));
                            stop=1;
                            break
                        else
                            s2m(i)=((M-n)*s2cap-um{i}'*Kminv{i}*um{i})/(M-n-Dim(i));
                            Qw{i}=Qm{i}+Bm{i}'*(Im+Am{i}'*Kminv{i}*Bm{i}')*Am{i}';
                            %Qwinv{i}=inv(Qw{i});
                            Qwinv{i}=Qw{i}\eye(length(Qw{i}));
                            deg1=Dim(i);
                            deg2=M-n-Dim(i);
                            F(i)=wm{i}'*Qwinv{i}*wm{i}/(Dim(i)*s2m(i));
                            Flim(i)=finv(significance,deg1,deg2);
                        end
                    end
                end
                if stop==0
                    %% apply final solution
                    % find maximum F(i)/Flim(i)
                    imax=find(abs(F./Flim)==max(abs(F./Flim)));
                    
                    if (abs(F(imax))<Flim(imax))
                        % no outlier
                        %fprintf('  -> No outlier found!\n');
                        %fprintf('        -> Chi quadro sper   : %.2f\n',s2cap/(0.5)^2*(M-n));
                        imax=0;
                        xfin=xcap;
                        dum=y;
                        Afin=A;
                        Qfin=Q;
                        s2fin=s2cap;
                        Ninvfin=Ninv;
                    else
                        %     % if the maximum ratio exceedes the threshold, the block is eliminated from the solution
                        %     %fprintf('  -> BLOCK %d found as outlier!\n',imax);
                        xfin=xcap-Bm{imax}*Kminv{imax}*um{imax};
                        dum=[y(1:Row(imax)-1);y(Row(imax)+Dim(imax):length(y))];
                        Afin1=A(1:Row(imax)-1,:);
                        Afin2=A(Row(imax)+Dim(imax):length(y),:);
                        Afin=[Afin1;Afin2];
                        Qfin11=Q(1:Row(imax)-1,1:Row(imax)-1);
                        Qfin12=Q(1:Row(imax)-1,Row(imax)+Dim(imax):length(y));
                        Qfin21=Q(Row(imax)+Dim(imax):length(y),1:Row(imax)-1);
                        Qfin22=Q(Row(imax)+Dim(imax):length(y),Row(imax)+Dim(imax):length(y));
                        Qfin=[Qfin11 Qfin12;Qfin21 Qfin22];
                        s2fin=s2m(imax);
                        Ninvfin=Ninv+Bm{imax}*Kminv{imax}*Bm{imax}';
                        
                        %fprintf('       -> S0 with block %2d   : %.3e\n',imax,s2cap^.5);
                        %fprintf('       -> S0 without block %2d: %.3e\n',imax,s2fin^.5);
                        
                    end
                    
                    yfin=Afin*xfin;
                    ufin=dum-yfin;
                    Cxx=s2fin*Ninvfin;
                    Cyy=Afin*Cxx*Afin';
                    Cuu=s2fin*Qfin-Cyy;
                else
                    imax=0;
                    xfin=xcap;
                    dum=y;
                    Afin=A;
                    Qfin=Q;
                    s2fin=s2cap;
                    Ninvfin=Ninv;
                    yfin=Afin*xfin;
                    ufin=dum-yfin;
                    Cxx=s2fin*Ninvfin;
                    Cyy=Afin*Cxx*Afin';
                    Cuu=s2fin*Qfin-Cyy;
                end
            else
                imax=0;
                xfin=xcap;
                dum=y;
                Afin=A;
                Qfin=Q;
                s2fin=s2cap;
                Ninvfin=Ninv;
                yfin=Afin*xfin;
                ufin=dum-yfin;
                Cxx=s2fin*Ninvfin;
                Cyy=Afin*Cxx*Afin';
                Cuu=s2fin*Qfin-Cyy;
            end
        end
        
        
        
        % Bancroft algorithm for the computation of ground coordinates
        % having at least 4 visible satellites.
        function pos = bancroft(matB)
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
            %
            % Copyright (C) Kai Borre
            % Kai Borre 04-30-95, improved by C.C. Goad 11-24-96
            %
            % Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
            %----------------------------------------------------------------------------------------------
            
            pos = zeros(4,1);   % Init position of the receiver
            
            nSat = size(matB,1);
            
            for iter = 1:2
                % Compute correctionon XS coordinates
                % due to the rotation of the Earth
                for s = 1:nSat
                    x = matB(s,1);      % X coordinate of the satellite S
                    y = matB(s,2);      % Y coordinate of the satellite S
                    if iter == 1
                        traveltime = 0.072; % [s]
                    else
                        z = matB(s,3);  % Z coordinate of the satellite S
                        % Compute approximate distance^2 between R and S
                        rho = (x-pos(1))^2+(y-pos(2))^2+(z-pos(3))^2;
                        traveltime = sqrt(rho)/goGNSS.V_LIGHT;
                    end
                    angle = traveltime*goGNSS.OMEGAE_DOT;
                    cosa = cos(angle);
                    sina = sin(angle);
                    % Applying correction (rotation)
                    matB(s,1) =	cosa*x + sina*y;
                    matB(s,2) = -sina*x + cosa*y;
                end % i-loop
                
                % Criptical way to implement Bancroft
                if nSat > 4
                    BBB = (matB'*matB)\matB';
                else
                    BBB = inv(matB);
                end
                e = ones(nSat,1);
                alpha = zeros(nSat,1);
                for s = 1:nSat
                    alpha(s) = lorentz(matB(s,:)',matB(s,:)')/2;
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
                for s = 1:2
                    possible_pos(:,s) = r(s)*BBBe+BBBalpha;
                    possible_pos(4,s) = -possible_pos(4,s);
                end
                
                abs_omc = zeros(2,1);
                for s =1:nSat
                    for i = 1:2
                        c_dt = possible_pos(4,i);
                        calc = norm(matB(s,1:3)' -possible_pos(1:3,i))+c_dt;
                        omc = matB(s,4)-calc;
                        abs_omc(i) = abs(omc);
                    end
                end; % j-loop
                
                % discrimination between roots
                if abs_omc(1) > abs_omc(2)
                    pos = possible_pos(:,2);
                else
                    pos = possible_pos(:,1);
                end
            end;
        end
    end
    
end
