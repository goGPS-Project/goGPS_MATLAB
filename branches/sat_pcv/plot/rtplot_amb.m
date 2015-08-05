function rtplot_amb (t, delta, stima_amb, sigma_amb, cs)

% SYNTAX:
%   rtplot_amb (t, delta, stima_amb, sigma_amb, cs)
%
% INPUT:
%   t = survey time (t=1,2,...)
%   delta =
%   stima_amb =
%   sigma_amb =
%   cs = boolean variable for cycle-slip
%
% DESCRIPTION:
%   Real-time plot of ambiguities.

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

if (t == 1)

   sat = find(stima_amb ~= 0);                  % satellites in view (not pivot)
   nsat = length(sat);                          % number of satellites in view

   dt = (1 : delta)';

   for i = 1 : nsat
      subfig(i) = subplot(5,3,i+6);

      %assessed N combination
      plot(t, stima_amb(sat(i)), 'b.-');
      hold on; grid on;

      %acceptability range
      plot(t, stima_amb(sat(i)) + sigma_amb(sat(i)),'r:');
      plot(t, stima_amb(sat(i)) - sigma_amb(sat(i)),'r:');
      hold off

      %satellite id
      set(subfig(i),'UserData',sat(i));

      %axes and title
      ax = axis;
      axis([dt(1) dt(delta) floor(ax(3)) ceil(ax(4))]);
      title(['SATELLITE ',num2str(sat(i))]);
   end
   
else
    
    nst = length(stima_amb);
    
    b1 = zeros(nst,delta);
    b2 = zeros(nst,delta);
    b3 = zeros(nst,delta);
    
    subfig = get(1,'Children');
    i = 1;
    for j = 1 : length(subfig)
        subfigTitle = get(get(subfig(j),'Title'),'String');
        if (length(subfigTitle) < 9 | ~strcmp(subfigTitle(1:9),'SATELLITE'))
            handleOff(i) = subfig(j);
            set(subfig(j),'HandleVisibility','off');
            i = i + 1;
        end
    end
    
    subfig = get(1,'Children');
    
    tLim = get(subfig(1),'XLim');
    dt = (tLim(1) : tLim(2))';
    
    for i = 1 : length(subfig)
        sat = get(subfig(i),'UserData');
        subobj = get(subfig(i),'Children');
        
        tData = get(subobj(end),'XData')';
        b1(sat,tData-dt(1)+1) = get(subobj(end),'YData');
        b2(sat,tData-dt(1)+1) = get(subobj(end-1),'YData') - b1(sat,tData-dt(1)+1);        
        
        if (length(subobj) > 3)
            tData = get(subobj(1),'XData')';
            b3(sat,tData-dt(1)+1) = get(subobj(1),'YData');
        end
    end
    
    if (t <= delta)
        b1(:,t) = stima_amb;
        b2(:,t) = sigma_amb;
        b3(:,t) = stima_amb .* cs;
        dt = (1 : delta)';
    else
        b1(:,1:end-1) = b1(:,2:end);
        b2(:,1:end-1) = b2(:,2:end);
        b3(:,1:end-1) = b3(:,2:end);
        b1(:,end) = stima_amb;
        b2(:,end) = sigma_amb;
        b3(:,end) = stima_amb .* cs;
        dt = (t-delta+1 : t)';
    end
    
    %----------------------------------------------------------------------------
    
    clf                                          % delete previous sub-figures
    
    sat = find(sum(b1,2) ~= 0);					% satellites in view (not pivot)
    nsat = length(sat);                          % number of satellites in view
    
    for i = 1 : nsat
        
        subfig(i) = subplot(5,3,i+6);
        
        j = find(b1(sat(i),:) ~= 0);
        k = find(b3(sat(i),:) ~= 0);
        
        %assessed N combination
        plot(dt(j), b1(sat(i),j), 'b.-');
        hold on; grid on;
        
        %acceptability range
        plot(dt(j), b1(sat(i),j) + b2(sat(i),j),'r:');
        plot(dt(j), b1(sat(i),j) - b2(sat(i),j),'r:');
        
        %cycle-slips
        plot(dt(k), b3(sat(i),k),'g.');
        hold off
        
        %satellite id
        set(subfig(i),'UserData',sat(i));
        
        yt = get(subfig(i),'YTick');
        set(subfig(i),'YTickLabel', sprintf('%.1f|',yt))
        
        %axes and title
        ax = axis;
        axis([dt(1) dt(delta) floor(ax(3)) ceil(ax(4))]);
        title(['SATELLITE ',num2str(sat(i))]);
    end
    
    for i = 1 : length(handleOff)
        set(handleOff(i),'HandleVisibility','on');
    end
end

%-------------------------------------------------------------------------------