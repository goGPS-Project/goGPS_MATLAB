function [iodp_mask, prn_mask] = load_prnmask(MT, msg)

% SYNTAX:
%   [iodp_mask, prn_mask] = load_prnmask(MT, msg);
%
% INPUT:
%   MT  = message type (MT) [vector]
%   msg = EGNOS messages strings [matrix]
%
% OUTPUT:
%   iodp_mask = PRN mask IODPs [vector]
%   prn_mask  = PRN mask [matrix]
%
% DESCRIPTION:
%   Load the PRN mask.

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

r_MT = find(MT == 1);

%decode the MT1 - PRN mask data
[iodp_mask1, prn_mask1] = ems2prnmask(msg(r_MT(1),:));
   
%store the first IODP and PRN mask
iodp_mask(1,1) = iodp_mask1;
prn_mask(1,:) = prn_mask1;

%counter for different IODPs
i = 1;

for j = 2 : length(r_MT) 
    
    [iodp_temp, mask_temp] = ems2prnmask(msg(r_MT(j),:));
     
    if ~isequal( iodp_mask(i,1), iodp_temp )
        
        %if two subsequent IODPs are different, then store the new IODP and
        %the new PRN mask
        i = i + 1;
        iodp_mask(i,1) = iodp_temp;
        prn_mask(i,:)  = mask_temp;
    end
end
