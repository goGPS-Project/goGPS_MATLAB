
%   CLASS Colormap
% =========================================================================
%
% DESCRIPTION
%   Class to manages utilities
%
% EXAMPLE
%   % set of static utilities
%   Colormap.diffAndPred
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, Giulio Tagliaferro, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Colormap < handle
    properties (Constant)
    end
    
    methods (Static)
    end
    
    methods (Static, Access = private)
        function [r, g, b] = lab2RGB(L, a, b)
            % Function to move from CIE 1976 Lab 2 RGB color space
            %
            % INPUT
            %   L, a, b     arrays of color coordinates in th Lab space
            %
            % OUTPUT
            %   r, g, b     arrays of primary color components
            %
            % SYNTAX
            %   [r, g, b] = Colormap.Lab2RGB(L, a, b)
            %   [rgb] = Colormap.Lab2RGB(L, a, b)
            
            y = ( L + 16 ) ./ 116;
            x = a ./ 500 + y;
            z = y - b ./ 200;
            
            pos = (x.^3 > 0.008856);
            x(pos) = x(pos).^3;
            x(~pos) = (x(~pos) - 16 / 116 ) ./ 7.787;

            pos = (y.^3 > 0.008856);
            y(pos) = y(pos).^3;
            y(~pos) = (y(~pos) - 16 / 116 ) ./ 7.787;
                        
            pos = (z.^3 > 0.008856);
            z(pos) = z(pos).^3;
            z(~pos) = (z(~pos) - 16 / 116 ) ./ 7.787;
            
            % Observer= 2?, Illuminant= D65
            ref_x =  95.047;
            ref_y = 100.000;
            ref_z = 108.883;
            
            X = ref_x .* x;
            Y = ref_y .* y;
            Z = ref_z .* z;
            
                                 % D65 white point
            x = X ./ 100;        % X from 0 to  95.0456
            y = Y ./ 100;        % Y from 0 to 100.0000
            z = Z ./ 100;        % Z from 0 to 108.8754
            
            r = x .*  3.240479 + y .* -1.537150 + z .* -0.498535;
            g = x .* -0.969256 + y .*  1.875992 + z .*  0.041556;
            b = x .*  0.055648 + y .* -0.204043 + z .*  1.057311;
            
            pos = ( r > 0.0031308 );
            r(pos) = 1.055 .* ( r(pos) .^ ( 1 / 2.4 ) ) - 0.055;
            r(~pos) = 12.92 .* r(~pos);
            
            pos = ( g > 0.0031308 );
            g(pos) = 1.055 .* ( g(pos) .^ ( 1 / 2.4 ) ) - 0.055;
            g(~pos) = 12.92 .* g(~pos);
            
            pos = ( b > 0.0031308 );
            b(pos) = 1.055 .* ( b(pos) .^ ( 1 / 2.4 ) ) - 0.055;
            b(~pos) = 12.92 .* b(~pos);
            
            r = max(0,min(r .* 255,255));
            g = max(0,min(g .* 255,255));
            b = max(0,min(b .* 255,255));
            
            if ((nargout == 1) || (nargout == 0))
                r = [r,g,b];
            end            
        end
    end
end
