% Copyright (c) 2019, Amir Khodabandeh and Peter J.G. Teunissen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Proper referencing must be done in any publications which use the 
% 	  source code and its functionalities:
% 		1. Teunissen P.J.G.(2019): A New GLONASS FDMA Model. GPS Solutions, doi:10.1007/s10291-019-0889-0.
% 		2. Teunissen P.J.G. and Khodabandeh A.(2019): GLONASS Ambiguity Resolution. GPS Solutions, doi:10.1007/s10291-019-0890-7.	
% 		
%     * Redistributions of the source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
% 	  
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.  
% 
% THIS SOFTWARE ROUTINE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS ROUTINE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [L, iL, Z, iZ] = GLONASS_L(kappa, opt, Dm)
%
%        [L, iL, Z, iZ] = GLONASS_L(kappa, opt, Dm);
%
% This routine allows you to construct the design matrix of the integer-estimable GLONASS FDMA model (REFS 1, 2). 
% It computes the GLONASS lower-triangular ambiguity design matrix L and its inverse, and optionally also the 
% GLONASS ambiguity-defining Z-transformation and its inverse.
%
% INPUT:
%
%    kappa : m-vector contaning channel numbers of GLONASS visible
%            satellites (m: the number of tracked satellites)
%
%    opt   : option for having supplementary matrices next to matrix L:
%            opt = 1 ==> only computes L as an output (Default)
%            opt = 2 ==> also computes supplementary matrices
%
%    Dm    : the m x (m-1) between-satellite differencing matrix (optional)
%            matrix Dm is only applicable for opt = 2 
%
%
% OUTPUT:
%       L  : the (m-1) x (m-1) lower-triangular matrix L
%      iL  : the (m-1) x (m-1) inverse-matrix of L
%       Z  : the corresponding m x m admissible integer transformation Z
%      iZ  : the m x m inverse-matrix of Z
%
%--------------------------------------------------------------------------
% DATE     : March-2019                                             
% Authors  : Amir Khodabandeh & Peter J.G. Teunissen                                          
%            Department of Infrastrucutre Engineering, 
%            The University of Melbourne, and  
%            Department of Geoscience and Remote Sensing,
%            Delft University of Technology
%--------------------------------------------------------------------------
%
% REFERENCES: 
%
%  1. Teunissen P.J.G.(2019). A New GLONASS FDMA Model. GPS Solutions, doi:10.1007/s10291-019-0889-0.
%  2. Teunissen P.J.G. and Khodabandeh A.(2019). GLONASS Ambiguity Resolution. GPS Solutions, doi:10.1007/s10291-019-0890-7.
%
%
%-----------------------------EXAMPLE--------------------------------
%
% input  ==>   kappa = [2, -4, -7, 3, -1, 0]';
%              opt   = 2;
% (optional)    Dm   = [-1    -1    -1    -1    -1
%                        1     0     0     0     0
%                        0     1     0     0     0
%                        0     0     1     0     0
%                        0     0     0     1     0
%                        0     0     0     0     1];
%--------------------------------------------------------------------
%              
% output ==> L = [ 0.0021   -0.0000   -0.0000         0         0
%                  0.0032    0.5012   -0.0000         0         0
%                 -0.0004    0.1665    0.3330         0         0
%                  0.0011   -0.5002   -1.0004    1.0004         0
%                  0.0007   -0.3333   -0.6667         0    1.0000]
%--------------------------------------------------------------------
%           iL = [474.3329    0.0000    0.0000    0.0000         0
%                  -2.9958    1.9951   -0.0000   -0.0000         0
%                   1.9972   -0.9975    3.0032    0.0000         0
%                   0.0000    0.0000    3.0032    0.9996         0
%                   0.0000    0.0000    2.0021    0.0000    1.0000]
%--------------------------------------------------------------------
%            Z =  [ 1    -475   -950    0     0  2850
%                   1    -474   -948    0     0  2844
%                   1    -473   -947    0     0  2841
%                   1    -475   -950    0     0  2851
%                   1    -475   -950    1     0  2847
%                   1    -475   -950    0     1  2848]
%--------------------------------------------------------------------
%            iZ = [-474   475     0     0     0     0
%                     1    -3     2     0     0     0
%                    -4     2    -1     3     0     0
%                    -4     0     0     3     1     0
%                    -3     0     0     2     0     1
%                    -1     0     0     1     0     0]
%--------------------------------------------------------------------

%%============================BEGIN PROGRAM==============================%%

%---------------------initialization------------------------------  
ze  = 2848;
a   = ze + kappa;
m   = length(kappa);

%Initializing the differencing matrix Dm
if(nargin < 2);  opt = 1;Dm = eye(m); Dm(1,:)= -1; Dm(:,1) = []; end

%In case Dm is absent, the default matrix is chosen
if (nargin < 3)
  disp('"The default differencing matrix Dm is chosen!"');
  Dm = eye(m); Dm(1,:)= -1; Dm(:,1) = [];
end

%check whether the option is specified correctly
if opt ~= 1; opt = 2; end

%-----------------------------------------------------------------


%Check whether the input channel numbers kappa are integer
if kappa ~= round(kappa)
   error ('Channel numbers ''kappa'' have to be integers!');
end

%Check whether the differencing matrix is given in its correct form
if size(Dm,1) < size(Dm,2)
   error ('The differencing matrix Dm must be of size m x (m-1)');
end

%-----------------------------------------------------------------

g  = a(1);
switch opt
    
    case 1     %default (option 1)
        Z  = []; 
        iZ = [];     
        L  = zeros(m-1);
        for i = 1:m-1
        
            %--computing the GCDs and Bezout coefficients alpha & beta
            [g(i+1),alpha] = gcd(a(i+1),g(i));
            alpha = -alpha;
            
            L(i,i) = g(i+1)/(g(i)*a(i+1));
            
            for j = i+1:m-1
                                                            
                L(j,i) = -(alpha*(a(j+1)-a(1)))/(g(i)*a(j+1));
                
            end
            
        end
        L  = ze*L; 
        iL = L\eye(m-1); 
        
    case 2     %option 2: also delivering the admissible integer transformations Z and iZ
        
        Z  = 1; 
        iZ = 1;      
        for i = 1:m-1
        
            %--computing the GCDs and Bezout coefficients alpha & beta
            [g(i+1),alpha, beta] = gcd(a(i+1),g(i));
            alpha = -alpha;
            
            %--computing the Z-matrix recursively
            hm  = (g(i)/g(i+1))*Z(:,end);
            aux = zeros(1,i); aux(end) = beta-(g(1)/g(i)); 
            wm  = Z(1,:) + aux;
            Zt  = Z * blkdiag(eye(i-1),alpha);
            Z   = [[Zt, hm];[wm, a(i+1)/g(i+1)]]; 
            
            %--computing the INVERSE of Z-matrix recursively
            v1    = zeros(i,1); v1(1)   = 1;
            vm    = zeros(i,1); vm(end) = 1;
            aux   = -((a(i+1)-a(1))/g(i+1))*iZ(end,:) - (g(i)/g(i+1))*v1';
            iZt   = [iZ(1:end-1,:);aux];
            ihm   = (g(i)/g(i+1))*vm;
            aux   = (g(i+1)+alpha*(a(i+1)-a(1)))/g(i);
            iwm   = alpha*v1' + aux*iZ(end,:);
            iZ    = [[iZt, ihm];[iwm,-alpha]];
            
        end
        iA       = diag(1./a);
        L        = ze*Dm'*iA*Z;
        L(:,end) = [];
              iL = L\eye(m-1);
    
        
end



return