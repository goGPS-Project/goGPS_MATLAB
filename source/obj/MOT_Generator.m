%   CLASS MOT_Generator
% =========================================================================
%
% DESCRIPTION
%   Class to manage Minor Ocean Tides
%
% EXAMPLE
%   ls = LS();
%


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro, Andrea Gatti ...
%  Contributors:      Giulio Tagliaferro, Andrea Gatti ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef MOT_Generator < handle
    
    properties(Constant)
        N_TIDES = 342;
        % Doodson numbers for the 342 tides
        DOODS_MOT = [2, 0, 0, 0, 0, 0;   2, 2,-2, 0, 0, 0;   2,-1, 0, 1, 0, 0;  ...
            2, 2, 0, 0, 0, 0;   2, 2, 0, 0, 1, 0;   2, 0, 0, 0,-1, 0; ...
            2,-1, 2,-1, 0, 0;   2,-2, 2, 0, 0, 0;   2, 1, 0,-1, 0, 0; ...
            2, 2,-3, 0, 0, 1;   2,-2, 0, 2, 0, 0;   2,-3, 2, 1, 0, 0; ...
            2, 1,-2, 1, 0, 0;   2,-1, 0, 1,-1, 0;   2, 3, 0,-1, 0, 0; ...
            2, 1, 0, 1, 0, 0;   2, 2, 0, 0, 2, 0;   2, 2,-1, 0, 0,-1; ...
            2, 0,-1, 0, 0, 1;   2, 1, 0, 1, 1, 0;   2, 3, 0,-1, 1, 0; ...
            2, 0, 1, 0, 0,-1;   2, 0,-2, 2, 0, 0;   2,-3, 0, 3, 0, 0; ...
            2,-2, 3, 0, 0,-1;   2, 4, 0, 0, 0, 0;   2,-1, 1, 1, 0,-1; ...
            2,-1, 3,-1, 0,-1;   2, 2, 0, 0,-1, 0;   2,-1,-1, 1, 0, 1; ...
            2, 4, 0, 0, 1, 0;   2,-3, 4,-1, 0, 0;   2,-1, 2,-1,-1, 0; ...
            2, 3,-2, 1, 0, 0;   2, 1, 2,-1, 0, 0;   2,-4, 2, 2, 0, 0; ...
            2, 4,-2, 0, 0, 0;   2, 0, 2, 0, 0, 0;   2,-2, 2, 0,-1, 0; ...
            2, 2,-4, 0, 0, 2;   2, 2,-2, 0,-1, 0;   2, 1, 0,-1,-1, 0; ...
            2,-1, 1, 0, 0, 0;   2, 2,-1, 0, 0, 1;   2, 2, 1, 0, 0,-1; ...
            2,-2, 0, 2,-1, 0;   2,-2, 4,-2, 0, 0;   2, 2, 2, 0, 0, 0; ...
            2,-4, 4, 0, 0, 0;   2,-1, 0,-1,-2, 0;   2, 1, 2,-1, 1, 0; ...
            2,-1,-2, 3, 0, 0;   2, 3,-2, 1, 1, 0;   2, 4, 0,-2, 0, 0; ...
            2, 0, 0, 2, 0, 0;   2, 0, 2,-2, 0, 0;   2, 0, 2, 0, 1, 0; ...
            2,-3, 3, 1, 0,-1;   2, 0, 0, 0,-2, 0;   2, 4, 0, 0, 2, 0; ...
            2, 4,-2, 0, 1, 0;   2, 0, 0, 0, 0, 2;   2, 1, 0, 1, 2, 0; ...
            2, 0,-2, 0,-2, 0;   2,-2, 1, 0, 0, 1;   2,-2, 1, 2, 0,-1; ...
            2,-1, 1,-1, 0, 1;   2, 5, 0,-1, 0, 0;   2, 1,-3, 1, 0, 1; ...
            2,-2,-1, 2, 0, 1;   2, 3, 0,-1, 2, 0;   2, 1,-2, 1,-1, 0; ...
            2, 5, 0,-1, 1, 0;   2,-4, 0, 4, 0, 0;   2,-3, 2, 1,-1, 0; ...
            2,-2, 1, 1, 0, 0;   2, 4, 0,-2, 1, 0;   2, 0, 0, 2, 1, 0; ...
            2,-5, 4, 1, 0, 0;   2, 0, 2, 0, 2, 0;   2,-1, 2, 1, 0, 0; ...
            2, 5,-2,-1, 0, 0;   2, 1,-1, 0, 0, 0;   2, 2,-2, 0, 0, 2; ...
            2,-5, 2, 3, 0, 0;   2,-1,-2, 1,-2, 0;   2,-3, 5,-1, 0,-1; ...
            2,-1, 0, 0, 0, 1;   2,-2, 0, 0,-2, 0;   2, 0,-1, 1, 0, 0; ...
            2,-3, 1, 1, 0, 1;   2, 3, 0,-1,-1, 0;   2, 1, 0, 1,-1, 0; ...
            2,-1, 2, 1, 1, 0;   2, 0,-3, 2, 0, 1;   2, 1,-1,-1, 0, 1; ...
            2,-3, 0, 3,-1, 0;   2, 0,-2, 2,-1, 0;   2,-4, 3, 2, 0,-1; ...
            2,-1, 0, 1,-2, 0;   2, 5, 0,-1, 2, 0;   2,-4, 5, 0, 0,-1; ...
            2,-2, 4, 0, 0,-2;   2,-1, 0, 1, 0, 2;   2,-2,-2, 4, 0, 0; ...
            2, 3,-2,-1,-1, 0;   2,-2, 5,-2, 0,-1;   2, 0,-1, 0,-1, 1; ...
            2, 5,-2,-1, 1, 0;   1, 1, 0, 0, 0, 0;   1,-1, 0, 0, 0, 0; ...
            1, 1,-2, 0, 0, 0;   1,-2, 0, 1, 0, 0;   1, 1, 0, 0, 1, 0; ...
            1,-1, 0, 0,-1, 0;   1, 2, 0,-1, 0, 0;   1, 0, 0, 1, 0, 0; ...
            1, 3, 0, 0, 0, 0;   1,-2, 2,-1, 0, 0;   1,-2, 0, 1,-1, 0; ...
            1,-3, 2, 0, 0, 0;   1, 0, 0,-1, 0, 0;   1, 1, 0, 0,-1, 0; ...
            1, 3, 0, 0, 1, 0;   1, 1,-3, 0, 0, 1;   1,-3, 0, 2, 0, 0; ...
            1, 1, 2, 0, 0, 0;   1, 0, 0, 1, 1, 0;   1, 2, 0,-1, 1, 0; ...
            1, 0, 2,-1, 0, 0;   1, 2,-2, 1, 0, 0;   1, 3,-2, 0, 0, 0; ...
            1,-1, 2, 0, 0, 0;   1, 1, 1, 0, 0,-1;   1, 1,-1, 0, 0, 1; ...
            1, 4, 0,-1, 0, 0;   1,-4, 2, 1, 0, 0;   1, 0,-2, 1, 0, 0; ...
            1,-2, 2,-1,-1, 0;   1, 3, 0,-2, 0, 0;   1,-1, 0, 2, 0, 0; ...
            1,-1, 0, 0,-2, 0;   1, 3, 0, 0, 2, 0;   1,-3, 2, 0,-1, 0; ...
            1, 4, 0,-1, 1, 0;   1, 0, 0,-1,-1, 0;   1, 1,-2, 0,-1, 0; ...
            1,-3, 0, 2,-1, 0;   1, 1, 0, 0, 2, 0;   1, 1,-1, 0, 0,-1; ...
            1,-1,-1, 0, 0, 1;   1, 0, 2,-1, 1, 0;   1,-1, 1, 0, 0,-1; ...
            1,-1,-2, 2, 0, 0;   1, 2,-2, 1, 1, 0;   1,-4, 0, 3, 0, 0; ...
            1,-1, 2, 0, 1, 0;   1, 3,-2, 0, 1, 0;   1, 2, 0,-1,-1, 0; ...
            1, 0, 0, 1,-1, 0;   1,-2, 2, 1, 0, 0;   1, 4,-2,-1, 0, 0; ...
            1,-3, 3, 0, 0,-1;   1,-2, 1, 1, 0,-1;   1,-2, 3,-1, 0,-1; ...
            1, 0,-2, 1,-1, 0;   1,-2,-1, 1, 0, 1;   1, 4,-2, 1, 0, 0; ...
            1,-4, 4,-1, 0, 0;   1,-4, 2, 1,-1, 0;   1, 5,-2, 0, 0, 0; ...
            1, 3, 0,-2, 1, 0;   1,-5, 2, 2, 0, 0;   1, 2, 0, 1, 0, 0; ...
            1, 1, 3, 0, 0,-1;   1,-2, 0, 1,-2, 0;   1, 4, 0,-1, 2, 0; ...
            1, 1,-4, 0, 0, 2;   1, 5, 0,-2, 0, 0;   1,-1, 0, 2, 1, 0; ...
            1,-2, 1, 0, 0, 0;   1, 4,-2, 1, 1, 0;   1,-3, 4,-2, 0, 0; ...
            1,-1, 3, 0, 0,-1;   1, 3,-3, 0, 0, 1;   1, 5,-2, 0, 1, 0; ...
            1, 1, 2, 0, 1, 0;   1, 2, 0, 1, 1, 0;   1,-5, 4, 0, 0, 0; ...
            1,-2, 0,-1,-2, 0;   1, 5, 0,-2, 1, 0;   1, 1, 2,-2, 0, 0; ...
            1, 1,-2, 2, 0, 0;   1,-2, 2, 1, 1, 0;   1, 0, 3,-1, 0,-1; ...
            1, 2,-3, 1, 0, 1;   1,-2,-2, 3, 0, 0;   1,-1, 2,-2, 0, 0; ...
            1,-4, 3, 1, 0,-1;   1,-4, 0, 3,-1, 0;   1,-1,-2, 2,-1, 0; ...
            1,-2, 0, 3, 0, 0;   1, 4, 0,-3, 0, 0;   1, 0, 1, 1, 0,-1; ...
            1, 2,-1,-1, 0, 1;   1, 2,-2, 1,-1, 0;   1, 0, 0,-1,-2, 0; ...
            1, 2, 0, 1, 2, 0;   1, 2,-2,-1,-1, 0;   1, 0, 0, 1, 2, 0; ...
            1, 0, 1, 0, 0, 0;   1, 2,-1, 0, 0, 0;   1, 0, 2,-1,-1, 0; ...
            1,-1,-2, 0,-2, 0;   1,-3, 1, 0, 0, 1;   1, 3,-2, 0,-1, 0; ...
            1,-1,-1, 0,-1, 1;   1, 4,-2,-1, 1, 0;   1, 2, 1,-1, 0,-1; ...
            1, 0,-1, 1, 0, 1;   1,-2, 4,-1, 0, 0;   1, 4,-4, 1, 0, 0; ...
            1,-3, 1, 2, 0,-1;   1,-3, 3, 0,-1,-1;   1, 1, 2, 0, 2, 0; ...
            1, 1,-2, 0,-2, 0;   1, 3, 0, 0, 3, 0;   1,-1, 2, 0,-1, 0; ...
            1,-2, 1,-1, 0, 1;   1, 0,-3, 1, 0, 1;   1,-3,-1, 2, 0, 1; ...
            1, 2, 0,-1, 2, 0;   1, 6,-2,-1, 0, 0;   1, 2, 2,-1, 0, 0; ...
            1,-1, 1, 0,-1,-1;   1,-2, 3,-1,-1,-1;   1,-1, 0, 0, 0, 2; ...
            1,-5, 0, 4, 0, 0;   1, 1, 0, 0, 0,-2;   1,-2, 1, 1,-1,-1; ...
            1, 1,-1, 0, 1, 1;   1, 1, 2, 0, 0,-2;   1,-3, 1, 1, 0, 0; ...
            1,-4, 4,-1,-1, 0;   1, 1, 0,-2,-1, 0;   1,-2,-1, 1,-1, 1; ...
            1,-3, 2, 2, 0, 0;   1, 5,-2,-2, 0, 0;   1, 3,-4, 2, 0, 0; ...
            1, 1,-2, 0, 0, 2;   1,-1, 4,-2, 0, 0;   1, 2, 2,-1, 1, 0; ...
            1,-5, 2, 2,-1, 0;   1, 1,-3, 0,-1, 1;   1, 1, 1, 0, 1,-1; ...
            1, 6,-2,-1, 1, 0;   1,-2, 2,-1,-2, 0;   1, 4,-2, 1, 2, 0; ...
            1,-6, 4, 1, 0, 0;   1, 5,-4, 0, 0, 0;   1,-3, 4, 0, 0, 0; ...
            1, 1, 2,-2, 1, 0;   1,-2, 1, 0,-1, 0;   0, 2, 0, 0, 0, 0; ...
            0, 1, 0,-1, 0, 0;   0, 0, 2, 0, 0, 0;   0, 0, 0, 0, 1, 0; ...
            0, 2, 0, 0, 1, 0;   0, 3, 0,-1, 0, 0;   0, 1,-2, 1, 0, 0; ...
            0, 2,-2, 0, 0, 0;   0, 3, 0,-1, 1, 0;   0, 0, 1, 0, 0,-1; ...
            0, 2, 0,-2, 0, 0;   0, 2, 0, 0, 2, 0;   0, 3,-2, 1, 0, 0; ...
            0, 1, 0,-1,-1, 0;   0, 1, 0,-1, 1, 0;   0, 4,-2, 0, 0, 0; ...
            0, 1, 0, 1, 0, 0;   0, 0, 3, 0, 0,-1;   0, 4, 0,-2, 0, 0; ...
            0, 3,-2, 1, 1, 0;   0, 3,-2,-1, 0, 0;   0, 4,-2, 0, 1, 0; ...
            0, 0, 2, 0, 1, 0;   0, 1, 0, 1, 1, 0;   0, 4, 0,-2, 1, 0; ...
            0, 3, 0,-1, 2, 0;   0, 5,-2,-1, 0, 0;   0, 1, 2,-1, 0, 0; ...
            0, 1,-2, 1,-1, 0;   0, 1,-2, 1, 1, 0;   0, 2,-2, 0,-1, 0; ...
            0, 2,-3, 0, 0, 1;   0, 2,-2, 0, 1, 0;   0, 0, 2,-2, 0, 0; ...
            0, 1,-3, 1, 0, 1;   0, 0, 0, 0, 2, 0;   0, 0, 1, 0, 0, 1; ...
            0, 1, 2,-1, 1, 0;   0, 3, 0,-3, 0, 0;   0, 2, 1, 0, 0,-1; ...
            0, 1,-1,-1, 0, 1;   0, 1, 0, 1, 2, 0;   0, 5,-2,-1, 1, 0; ...
            0, 2,-1, 0, 0, 1;   0, 2, 2,-2, 0, 0;   0, 1,-1, 0, 0, 0; ...
            0, 5, 0,-3, 0, 0;   0, 2, 0,-2, 1, 0;   0, 1, 1,-1, 0,-1; ...
            0, 3,-4, 1, 0, 0;   0, 0, 2, 0, 2, 0;   0, 2, 0,-2,-1, 0; ...
            0, 4,-3, 0, 0, 1;   0, 3,-1,-1, 0, 1;   0, 0, 2, 0, 0,-2; ...
            0, 3,-3, 1, 0, 1;   0, 2,-4, 2, 0, 0;   0, 4,-2,-2, 0, 0; ...
            0, 3, 1,-1, 0,-1;   0, 5,-4, 1, 0, 0;   0, 3,-2,-1,-1, 0; ...
            0, 3,-2, 1, 2, 0;   0, 4,-4, 0, 0, 0;   0, 6,-2,-2, 0, 0; ...
            0, 5, 0,-3, 1, 0;   0, 4,-2, 0, 2, 0;   0, 2, 2,-2, 1, 0; ...
            0, 0, 4, 0, 0,-2;   0, 3,-1, 0, 0, 0;   0, 3,-3,-1, 0, 1; ...
            0, 4, 0,-2, 2, 0;   0, 1,-2,-1,-1, 0;   0, 2,-1, 0, 0,-1; ...
            0, 4,-4, 2, 0, 0;   0, 2, 1, 0, 1,-1;   0, 3,-2,-1, 1, 0; ...
            0, 4,-3, 0, 1, 1;   0, 2, 0, 0, 3, 0;   0, 6,-4, 0, 0, 0]
        % id of the 11 pricipal tides
        ID_11_PRINCIPAL = [1 2 3 4 110 111 112 113 264 265 266];
        % 342 tides amplitude
        T_AMP = [.632208, .294107, .121046, .079915, .023818,-.023589, .022994, ...
            .019333,-.017871, .017192, .016018, .004671,-.004662,-.004519, ...
            .004470, .004467, .002589,-.002455,-.002172, .001972, .001947, ...
            .001914,-.001898, .001802, .001304, .001170, .001130, .001061, ...
            -.001022,-.001017, .001014, .000901,-.000857, .000855, .000855, ...
            .000772, .000741, .000741,-.000721, .000698, .000658, .000654, ...
            -.000653, .000633, .000626,-.000598, .000590, .000544, .000479, ...
            -.000464, .000413,-.000390, .000373, .000366, .000366,-.000360, ...
            -.000355, .000354, .000329, .000328, .000319, .000302, .000279, ...
            -.000274,-.000272, .000248,-.000225, .000224,-.000223,-.000216, ...
            .000211, .000209, .000194, .000185,-.000174,-.000171, .000159, ...
            .000131, .000127, .000120, .000118, .000117, .000108, .000107, ...
            .000105,-.000102, .000102, .000099,-.000096, .000095,-.000089, ...
            -.000085,-.000084,-.000081,-.000077,-.000072,-.000067, .000066, ...
            .000064, .000063, .000063, .000063, .000062, .000062,-.000060, ...
            .000056, .000053, .000051, .000050, .368645,-.262232,-.121995, ...
            -.050208, .050031,-.049470, .020620, .020613, .011279,-.009530, ...
            -.009469,-.008012, .007414,-.007300, .007227,-.007131,-.006644, ...
            .005249, .004137, .004087, .003944, .003943, .003420, .003418, ...
            .002885, .002884, .002160,-.001936, .001934,-.001798, .001690, ...
            .001689, .001516, .001514,-.001511, .001383, .001372, .001371, ...
            -.001253,-.001075, .001020, .000901, .000865,-.000794, .000788, ...
            .000782,-.000747,-.000745, .000670,-.000603,-.000597, .000542, ...
            .000542,-.000541,-.000469,-.000440, .000438, .000422, .000410, ...
            -.000374,-.000365, .000345, .000335,-.000321,-.000319, .000307, ...
            .000291, .000290,-.000289, .000286, .000275, .000271, .000263, ...
            -.000245, .000225, .000225, .000221,-.000202,-.000200,-.000199, ...
            .000192, .000183, .000183, .000183,-.000170, .000169, .000168, ...
            .000162, .000149,-.000147,-.000141, .000138, .000136, .000136, ...
            .000127, .000127,-.000126,-.000121,-.000121, .000117,-.000116, ...
            -.000114,-.000114,-.000114, .000114, .000113, .000109, .000108, ...
            .000106,-.000106,-.000106, .000105, .000104,-.000103,-.000100, ...
            -.000100,-.000100, .000099,-.000098, .000093, .000093, .000090, ...
            -.000088, .000083,-.000083,-.000082,-.000081,-.000079,-.000077, ...
            -.000075,-.000075,-.000075, .000071, .000071,-.000071, .000068, ...
            .000068, .000065, .000065, .000064, .000064, .000064,-.000064, ...
            -.000060, .000056, .000056, .000053, .000053, .000053,-.000053, ...
            .000053, .000053, .000052, .000050,-.066607,-.035184,-.030988, ...
            .027929,-.027616,-.012753,-.006728,-.005837,-.005286,-.004921, ...
            -.002884,-.002583,-.002422, .002310, .002283,-.002037, .001883, ...
            -.001811,-.001687,-.001004,-.000925,-.000844, .000766, .000766, ...
            -.000700,-.000495,-.000492, .000491, .000483, .000437,-.000416, ...
            -.000384, .000374,-.000312,-.000288,-.000273, .000259, .000245, ...
            -.000232, .000229,-.000216, .000206,-.000204,-.000202, .000200, ...
            .000195,-.000190, .000187, .000180,-.000179, .000170, .000153, ...
            -.000137,-.000119,-.000119,-.000112,-.000110,-.000110, .000107, ...
            -.000095,-.000095,-.000091,-.000090,-.000081,-.000079,-.000079, ...
            .000077,-.000073, .000069,-.000067,-.000066, .000065, .000064, ...
            -.000062, .000060, .000059,-.000056, .000055,-.000051];
    end
    
    methods
        function this = MOT_Generator()
            % Initialisation of the variables
            %
            % SYNTAX
            %   this = MOT_Generator()
        end
    end

    methods (Static)
        function [freq] = getDoodsonFr(dod_num, gps_time)
            % get the frequency of the tides based on its doodson number
            %
            % SYNTAX:
            %   [freq] = this.getDoodsonFr(dod_num, gps_time)
            jd = gps_time.getJD;
            t = (jd - 2451545)/36525; % century from JD 2000.0
            %  Find frequencies of Delauney variables (in cycles/day), and from these
            %  the same for the Doodson arguments
            
            fd1 =  0.0362916471 +  0.0000000013*t;
            fd2 =  0.0027377786;
            fd3 =  0.0367481951 - 0.0000000005*t;
            fd4 =  0.0338631920 - 0.0000000003*t;
            fd5 = -0.0001470938 +  0.0000000003*t;
            dd = zeros(6,1);
            dd(1) = 1.0 - fd4;
            dd(2) = fd3 +  fd5;
            dd(3) = dd(2) - fd4;
            dd(4) = dd(2) - fd1;
            dd(5) = -fd5;
            dd(6) = dd(3) - fd2;
            
            %  End of intialization (likely to be called only once)
            %  Compute phase and frequency of the given tidal constituent
            
            freq = dod_num*dd;
            
        end
        
        function [phase] = getDoodsonPh(dod_num, gps_time)
            % get the phase of the tides based on its doodson number
            %
            % SYNTAX:
            %   [freq] = this.getDoodsonPh(dod_num, gps_time)
            jd = gps_time.getJD;
            t = (jd - 2451545)/36525; % century from JD 2000.0
            [~,~,sod] = gps_time.getDOY; %fractiona part of the day <_ check might be not true
            day_fr = sod / 86400;
            f1 =     134.9634025100 +  t.*( 477198.8675605000 +  t.*(      0.0088553333 +   t.*(      0.0000143431 +   t.*(     -0.0000000680 ))));
            f2 =     357.5291091806 +  t.*(  35999.0502911389 +  t.*(     -0.0001536667 +   t.*(      0.0000000378 +   t.*(     -0.0000000032 ))));
            f3 =      93.2720906200 +  t.*( 483202.0174577222 +  t.*(     -0.0035420000 +   t.*(     -0.0000002881 +   t.*(      0.0000000012 ))));
            f4 =     297.8501954694 +  t.*( 445267.1114469445 +  t.*(     -0.0017696111 +   t.*(      0.0000018314 +   t.*(     -0.0000000088 ))));
            f5 =     125.0445550100 +  t.*(  -1934.1362619722 +  t.*(      0.0020756111 +   t.*(      0.0000021394 +   t.*(     -0.0000000165 ))));
            
            % Convert to Doodson (Darwin) variables
            d = zeros(6, gps_time.length);
            d(1,:) = 360.0 * day_fr - f4;
            d(2,:) = f3 +  f5;
            d(3,:) = d(2,:) - f4';
            d(4,:) = d(2,:) - f1';
            d(5,:) = -f5';
            d(6,:) = d(3,:) - f2';
            
            %  End of intialization (likely to be called only once)
            %  Compute phase and frequency of the given tidal constituent
            
            phase = dod_num * d;
            
            % Adjust phases so that they fall in the positive range 0 to 360
            phase = mod(phase, 360);
            if phase < 0
                phase = phase + 360;
            end
        end
        
        function phs_out = getDoodsonPhs(gps_time)
            % get all phases
            %
            % SYNTAX:
            %     phs_out = this.getDoodsonPhs(gps_time)
            phs_out = zeros(MOT_Generator.N_TIDES,gps_time.length);
            for t = 1 : MOT_Generator.N_TIDES
                phs_out(t,:) = MOT_Generator.getDoodsonPh(MOT_Generator.DOODS_MOT(t,:), gps_time);
            end
        end
        
        function [amp_out, ph_out, f_out] = admitanceInterp(amp_in, ph_in, gps_time)
            % given the local apmlitude and phese component at the
            % location for the 11 main component retun the inteprolated
            % amplitude and phase fro all the 342 componenet
            % NOTE: Porting of the IERS ADMINT.F function
            %
            % SYNTAX:
            %    [amp_out, ph_out, f_out]= this.admitanceInterp(amp_in, ph_in)
            admittance_real = zeros(MOT_Generator.N_TIDES,1);
            admittance_im = zeros(MOT_Generator.N_TIDES,1);
            f_out = zeros(MOT_Generator.N_TIDES,1);
            for t = 1 : MOT_Generator.N_TIDES
                f_out(t) = MOT_Generator.getDoodsonFr(MOT_Generator.DOODS_MOT(t,:), gps_time);
            end
            % compute the admittance real and imaginary part (wathever it
            % means ..) for the 11 princiapl components
            admittance_real(MOT_Generator.ID_11_PRINCIPAL) = amp_in .* cosd(ph_in) ./ abs(MOT_Generator.T_AMP(MOT_Generator.ID_11_PRINCIPAL));
            admittance_im(MOT_Generator.ID_11_PRINCIPAL) = amp_in .* sind(ph_in) ./ abs(MOT_Generator.T_AMP(MOT_Generator.ID_11_PRINCIPAL));
            lid_11_principal = false(MOT_Generator.N_TIDES,1);
            lid_11_principal(MOT_Generator.ID_11_PRINCIPAL) = true;
            for i = 0:2
                id_fr_part = MOT_Generator.DOODS_MOT(:,1) == i;
                id_int_trg = id_fr_part & ~lid_11_principal;
                lid_int_source = id_fr_part & lid_11_principal;
                %                 min_f = min(f_out(id_fr_part));
                %                 max_f = max(f_out(id_fr_part));
                [min_f, id_min] = min(f_out(lid_int_source));
                id_int_source = find(lid_int_source);
                id_min = id_int_source(id_min);
                [max_f, id_max] = max(f_out(lid_int_source));
                id_max = id_int_source(id_max);
                admittance_real(id_int_trg) = interp1(f_out(lid_int_source),admittance_real(lid_int_source),f_out(id_int_trg),'spline');
                admittance_real(f_out < min_f & id_fr_part) = admittance_real(id_min);
                admittance_real(f_out > max_f & id_fr_part) = admittance_real(id_max);
                admittance_im(id_int_trg) = interp1(f_out(lid_int_source),admittance_im(lid_int_source),f_out(id_int_trg),'spline');
                admittance_im(f_out < min_f & id_fr_part) = admittance_im(id_min);
                admittance_im(f_out > max_f & id_fr_part) = admittance_im(id_max);
            end
            amp_out = MOT_Generator.T_AMP' .* sqrt(admittance_real.^2 + admittance_im.^2);
            ph_out = MOT_Generator.getDoodsonPhs(gps_time);
            ph_out(MOT_Generator.DOODS_MOT(:,1) == 0) = ph_out(MOT_Generator.DOODS_MOT(:,1) == 0) + 180;
            ph_out(MOT_Generator.DOODS_MOT(:,1) == 1) = ph_out(MOT_Generator.DOODS_MOT(:,1) == 1) + 90;
            ph_out = rem(ph_out - atan2(admittance_im, admittance_real)*180/pi,360);
                        
            ph_out(ph_out > 180) = ph_out(ph_out > 180) - 360;
        end
        
        function [disp] = computeDisplacement(amp_in, ph_in, gps_time)
            % given amp and phase of the 11 major tides compute the
            % displacement using all 342 tides
            %
            % SYNTAX:
            %      [disp] = this.computeDisplacement(amp_in, ph_in, gps_time)
            disp = zeros(gps_time.length, 3);
            %[~,~,sod] = gps_time.getDOY;
            sod = gps_time - gps_time.first;
            %ph_dood = getDoodsonPhs(this, gps_time);
            for i = 1 : 3
                [amp_out, ph_out, f_out] = MOT_Generator.admitanceInterp(amp_in(i,:),  ph_in(i,:), gps_time.first);
                for t = 1 : MOT_Generator.N_TIDES
                    %ph = ph_dood(t,:) + ph_out(t);
                    ph = sod'/43200*180*f_out(t) + ph_out(t);
                    disp (:,i) = disp (:,i) + amp_out(t) .* cos(ph/180*pi)';
                end
            end
        end
    end
end
