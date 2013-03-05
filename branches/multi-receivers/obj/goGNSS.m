% =========================================================================
%   OBJECT goGNSS
% =========================================================================
%
% DESCRIPTION
%   Object to manage a useful functions / standard parameters of goGPS
%
% EXAMPLE
%   go = goGNSS();
%
% LIST of CONSTANT
%
%   V_LIGHT         velocity of light in the void                 [m/s]
%   F1              GPS carriers frequencies F1                   [1/s]
%   F2              GPS carriers frequencies F2                   [1/s]
%   A               ellipsoid semi-major axis                     [m]
%   F               ellipsoid flattening
%   E               eccentricity
%   GM              gravitational constant (mass of Earth)        [m^3/s^2]
%   OMEGAE_DOT      angular velocity of the Earth rotation        [rad/s]
%
% LIST of METHODS
%
%   [pos] = bancroft(B_pass)
%
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti
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
        % velocity of light in the void
        V_LIGHT = 299792458; % [m/s]
        
        % GPS carriers frequencies
        F1 = 1575420000; % [1/s]
        F2 = 1227600000; % [1/s]
        
        % GPS carriers wavelengths
        LAMBDA1 = goGNSS.V_LIGHT / goGNSS.F1; % [m]
        LAMBDA2 = goGNSS.V_LIGHT / goGNSS.F2; % [m]
        
        % ellipsoid semi-major axis
        ELL_A = 6378137; % [m]
        
        % ellipsoid flattening
        ELL_F = 1/298.257222101;
        
        % eccentricity
        ELL_E = sqrt(1-(1-goGNSS.ELL_F)^2);
         
        % gravitational constant (mass of Earth)
        GM = 3.986004418e14; % [m^3/s^2]
        
        % angular velocity of the Earth rotation
        OMEGAE_DOT = 7.2921151467e-5; % [rad/s]
    end
    
    % Creator (empty)
    methods
        function obj = goGNSS()
        end
    end
    
    % Useful function (static)
    methods (Static, Access = 'public')   
        
        % Function to get a bancroft solution from one receiver
        % 
        % Typical old way of calling the function
        % index = find(no_eph == 0);
        % pos = getBancroftPos(XS(index,:), dtS(index), pseudorange(index))
        function pos = getBancroftPos(XS, dtS, prR)
            matB = [XS(:,:), prR(:) + goGNSS.V_LIGHT * dtS(:)]; % Bancroft matrix
            pos = obj.bancroft(matB);            
        end        
    end
    
    methods (Static, Access = 'private')
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
            %----------------------------------------------------------------------------------------------            pos = zeros(4,1);
            
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
