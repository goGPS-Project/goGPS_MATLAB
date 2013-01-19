function rtplot_snr(snr)

% SYNTAX:
%   rtplot_snr (snr);
%
% INPUT:
%   snr = signal-to-noise ratio
%
% DESCRIPTION:
%   Real time bar graph of signal-to-noise ratio.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

subplot(2,3,6)

sat = find(snr > 0);
snr = snr(sat);

if ~isempty(sat)
    barh(1:1:size(sat),snr)
    axis([0 60 0 size(sat,1)+1])
    set(gca,'XTick',[0 10 20 30 40 50 60])
    set(gca,'YTickLabel',sat)
    grid on
end

title('signal-to-noise ratio')

%-------------------------------------------------------------------------------