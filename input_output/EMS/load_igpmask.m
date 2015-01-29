function [iodi_mask, band_mask, igp_mask, n_bands_mask] = load_igpmask(MT, msg)

% SYNTAX:
%   [iodi_mask, band_mask, igp_mask, n_bands_mask] = load_igpmask(MT, msg);
%
% INPUT:
%   MT  = message types [vector]
%   msg = EGNOS messages strings [matrix]
%
% OUTPUT:
%   iodi_mask = IODI of the IGP masks [vector]
%   band_mask = band number [vector]
%   igp_mask  = band mask [matrix]
%   n_bands_mask = number of the transmitted bands [scalar]
%
% DESCRIPTION:
%   Load the Iono Grid Point (IGP) masks referring to each band (MT 18).

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

r_MT = find(MT == 18);

%decode the MT18 - IGP mask
for i = 1 : length(r_MT)
    
    [n_bands18, band18, iodi18, igp_mask18] = ems2igpmask(msg(r_MT(i),:)); %#ok<ASGLU>
    
    %n_bands(i,1) = n_bands18;
    band(i,1)  = band18;
    iodi(i,1)  = iodi18;
    mask(i,:) = igp_mask18;
   
end

%check which bands are transmitted
bands = unique(band);

b = 1;

%asociate IGP masks to each band
for i = 1 : length(bands)
    
    %message referring to each band
    r_band = find(band == bands(i));
    
    %store the first IODI and the first IGP mask for each band
    iodi_mask(b,1) = iodi(r_band(1),1);
    igp_mask(b, :) = mask(r_band(1),:);
    band_mask(b,1) = band(r_band(1),1);
    
    for j = 1 : length(r_band)-1
        
        if (~isequal(iodi(r_band(j),1), iodi(r_band(j+1),1)))
            
            %if two subsequent IODI are different, then store the new IODI
            %and the new IGP mask
            b = b + 1;
            iodi_mask(b,1) = iodi(r_band(j+1),1);
            igp_mask(b, :) = mask(r_band(j+1),:);
            band_mask(b,1) = band(r_band(j+1),1);
        end
    end
    
    %next band
    b = b + 1;
end
    
n_bands_mask = length(bands);    

%DEBUG
%disp(' ')
%disp(['IODI: ', num2str(iodi_mask') ]) 
%disp(['bands: ', num2str(band_mask') ]) 
