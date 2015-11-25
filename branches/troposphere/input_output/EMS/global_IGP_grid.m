function [grid_igp] = global_IGP_grid

% SYNTAX:
%   [grid_igp] = global_IGP_grid;
%
% INPUT:
%
% OUTPUT:
%   grid_igp = global IGP grid, composed by the various bands
%
% DESCRIPTION:
%   Global IGP grid production tool.
%
%   NOTE: bands are from 0 to 10, but this tool only uses bands 3,4,5,6,9,
%         i.e. those that cover Europe

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
%----------------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------------

%each band is composed by a [nx4] matrix
% n:        number of points in the band (e.g. 201 or 192)
% column 1: number of the band
% column 2: number of the IGP point in the band (see tab. A-14 pag. 31 appendix A, RTCA229C)
% column 3: latitude
% column 4: longitude
%
% the global matrix is composed by the various band matrices

%%%%%%%%%%
% BAND 3 %
%%%%%%%%%%

c13 = 3 * ones(201,1);
c23 = [1 : 201]';
%lat
c33a = [-75 -65 -55:5:55 65 75]'; 
c33b = [-55:5:55]';
c33c = [-85 -75 -65 -55:5:55 65 75]';
%assemble
c33 = [c33a; c33b; c33c; c33b; c33a; c33b; c33a; c33b];

%lon
c43a = -60 * ones(27,1); 
c43b = -55 * ones(23,1); 
c43c = -50 * ones(28,1); 
c43d = -45 * ones(23,1); 
c43e = -40 * ones(27,1); 
c43f = -35 * ones(23,1); 
c43g = -30 * ones(27,1); 
c43h = -25 * ones(23,1); 
%assemble
c43 = [c43a; c43b; c43c; c43d; c43e; c43f; c43g; c43h];

grid3 = [c13 c23 c33 c43];


%%%%%%%%%%
% BAND 4 %
%%%%%%%%%%

c14 = 4 * ones(201,1);
c24 = [1 : 201]';
%lat
c34a = [-75 -65 -55:5:55 65 75]'; 
c34b = [-55:5:55]';
c34c = [-75 -65 -55:5:55 65 75 85]';
%assemble
c34 = [c34a; c34b; c34a; c34b; c34c; c34b; c34a; c34b];

%lon
c44a = -20 * ones(27,1); 
c44b = -15 * ones(23,1); 
c44c = -10 * ones(27,1); 
c44d =  -5 * ones(23,1); 
c44e =   0 * ones(28,1); 
c44f =   5 * ones(23,1); 
c44g =  10 * ones(27,1); 
c44h =  15 * ones(23,1); 
%assemble
c44 = [c44a; c44b; c44c; c44d; c44e; c44f; c44g; c44h];

grid4 = [c14 c24 c34 c44];


%%%%%%%%%%
% BAND 5 %
%%%%%%%%%%

c15 = 5 * ones(201,1);
c25 = [1 : 201]';
%lat
c35a = [-75 -65 -55:5:55 65 75]'; 
c35b = [-55:5:55]';
c35c = [-85 -75 -65 -55:5:55 65 75]';
%assemble
c35 = [c35a; c35b; c35a; c35b; c35c; c35b; c35a; c35b];

%lon
c45a = 20 * ones(27,1); 
c45b = 25 * ones(23,1); 
c45c = 30 * ones(27,1); 
c45d = 35 * ones(23,1); 
c45e = 40 * ones(28,1); 
c45f = 45 * ones(23,1); 
c45g = 50 * ones(27,1); 
c45h = 55 * ones(23,1); 
%assemble
c45 = [c45a; c45b; c45c; c45d; c45e; c45f; c45g; c45h];

grid5 = [c15 c25 c35 c45];


%%%%%%%%%%
% BAND 6 %
%%%%%%%%%%

c16 = 6 * ones(201,1);
c26 = [1 : 201]';
%lat
c36a = [-75 -65 -55:5:55 65 75]'; 
c36b = [-55:5:55]';
c36c = [-75 -65 -55:5:55 65 75 85]';
%assemble
c36 = [c36a; c36b; c36a; c36b; c36a; c36b; c36c; c36b];

%lon
c46a = 60 * ones(27,1); 
c46b = 65 * ones(23,1); 
c46c = 70 * ones(27,1); 
c46d = 75 * ones(23,1); 
c46e = 80 * ones(27,1); 
c46f = 85 * ones(23,1); 
c46g = 90 * ones(28,1); 
c46h = 95 * ones(23,1); 
%assemble
c46 = [c46a; c46b; c46c; c46d; c46e; c46f; c46g; c46h];

grid6 = [c16 c26 c36 c46];


%%%%%%%%%%
% BAND 9 %
%%%%%%%%%%

c19 = 9 * ones(192,1);
c29 = [1 : 192]';
%lat
c39a = 60 * ones(72,1);
c39b = 65 * ones(36,1);
c39c = 70 * ones(36,1);
c39d = 75 * ones(36,1);
c39e = 85 * ones(12,1);
%assemble
c39 = [c39a; c39b; c39c; c39d; c39e];

%lon
c49a = [-180 : 5 : 175]';
c49b = [-180 : 10 : 170]';
c49c = [-180 : 30 : 150]';
%assemble
c49 = [c49a; c49b; c49b; c49b; c49c];

grid9 = [c19 c29 c39 c49];


%%%%%%%%%%%%%%%%%%%
% GLOBAL IGP GRID %
%%%%%%%%%%%%%%%%%%%

grid_igp = [grid3; grid4; grid5; grid6; grid9];
