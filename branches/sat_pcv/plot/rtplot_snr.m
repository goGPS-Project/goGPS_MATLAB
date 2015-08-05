function rtplot_snr(snr, Eph, SP3)

% SYNTAX:
%   rtplot_snr (snr, Eph, SP3);
%
% INPUT:
%   snr = signal-to-noise ratio
%   Eph = matrix containing 33 navigation parameters for each satellite
%   SP3 = structure containing precise ephemeris data
%
% DESCRIPTION:
%   Real time bar graph of signal-to-noise ratio.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

if (isempty(SP3))
    eph_avail = Eph(30,:);
    sys = Eph(31,:);
    prn = Eph(1,:);
else
    eph_avail = SP3.avail';
    sys = SP3.sys';
    prn = SP3.prn';
end

sat = find(snr > 0);
snr = snr(sat);

[~, idx1, idx2] = intersect(eph_avail, sat);
sat = sat(idx2);
snr = snr(idx2);

sys = char(sys(idx1));
prn = prn(idx1);

sysprn = {[sys' num2str(prn')]};

if ~isempty(sat)
    barh(1:1:size(sat),snr)
    axis([0 60 0 size(sat,1)+1])
    set(gca,'XTick',[0 10 20 30 40 50 60])
    set(gca,'YTick',1:numel(prn))
    set(gca,'YTickLabel',sysprn)
    set(gca,'FontSize',9)
    grid on
end

title('signal-to-noise ratio','FontSize',10)

%-------------------------------------------------------------------------------