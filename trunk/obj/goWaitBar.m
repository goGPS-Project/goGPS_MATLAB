% =========================================================================
%   OBJECT goWaitBar
% =========================================================================
%
% DESCRIPTION:
%   Object to show and manage a waitbar
%
% EXAMPLE:
%   goWB = goWaitBar(10);
%   goWB.go();
%
% LIST of METHODS
%
%  GENERIC ------------------------------------------------------------
%
%   init(obj, nSteps, msg)  Init the waitbar (and display it
%   close(obj)              Close the window
%
%  DISPLAY ------------------------------------------------------------
%   
%   go(obj, step)           Just update the waitbar
%
%   goMsg(obj, step)        Update the waitbar and accept a message to
%                           display within the window
%
%   goTime(obj, step)       Update the waitbar and estimate the remaining
%                           computational time supposing a linear trend
%
%   titleUpdate(obj, msg)   Change the title of the waitbar
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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
classdef goWaitBar < handle
            
    properties (GetAccess = 'public', SetAccess = 'public')
        h  = [];        % handle of the waitbar        
        t0 = 0;         % start time of process
        nSteps = 0;     % number of step of the bar
        lastStep = 0;   % Last step done
    end
    
    methods

        % Creator
        function obj = goWaitBar(nSteps, msg)
            if nargin == 1
                msg = 'Please wait...';
            end
            obj.nSteps = nSteps;
            obj.t0 = tic;
            fs = get(0,'defaultTextFontSize');
            set(0,'defaultTextFontSize', 12);
            obj.h = waitbar(0, msg);
            titleHandle = get(findobj(obj.h,'Type','axes'),'Title');
            set(titleHandle,'FontSize',12) ;
            set(0,'defaultTextFontSize', fs);
            drawnow;
        end
        
        % Just update the waitbar
        function go(obj,step)
            if nargin == 2
                obj.lastStep = min(obj.nSteps, step);
            else
                obj.lastStep = min(obj.lastStep + 1,obj.nSteps);
            end
            waitbar(min(1,obj.lastStep/obj.nSteps));
        end
        
        % Update the waitbar and accept a message to display within the window
        function goMsg(obj, step, msg)
            if nargin == 3
                obj.lastStep = min(obj.nSteps, step);
            else
                msg = step;
                obj.lastStep = min(obj.lastStep + 1,obj.nSteps);
            end
            obj.lastStep = min(obj.nSteps, obj.lastStep);
            waitbar(min(obj.lastStep/obj.nSteps),obj.h, msg);
        end        

        % Update the waitbar and estimate the remaining computational time supposing a linear trend
        function goTime(obj,step)
            if nargin == 2
                obj.lastStep = min(obj.nSteps, step);
            else
                obj.lastStep = min(obj.lastStep + 1,obj.nSteps);
            end
            % Elapsed time:
            t1= toc(obj.t0);
            hh = floor(t1/3600);
            mm = floor((t1-hh*3600)/60);
            ss = floor(t1-hh*3600-mm*60);
            elapsedTime = [num2str(hh,'%02d') ':' num2str(mm,'%02d') ':' num2str(ss,'%02d')];
            % Remaining Time
            t1 = t1/obj.lastStep*(obj.nSteps-obj.lastStep);
            hh = floor(t1/3600);
            mm = floor((t1-hh*3600)/60);
            ss = floor(t1-hh*3600-mm*60);
            remainingTime = [num2str(hh,'%02d') ':' num2str(mm,'%02d') ':' num2str(ss,'%02d')];
            waitbar(min(1,obj.lastStep/obj.nSteps),obj.h,[' Elapsed time                ' elapsedTime 10 ' Remaining time            ' remainingTime] );
        end        
       
        % Change the title of the waitbar
        function titleUpdate(obj, msg)
            set(obj.h,'Name', msg);
        end
        
        % Shift the waitbar downward (e.g. to make room for processing plots)
        function shiftDown(obj, shf)
            if (nargin < 2)
                shf = 120;
            end
            pos = get(obj.h,'Position');
            if (pos(2) > shf + pos(4))
                pos(2) = pos(2) - shf;
            end
            set(obj.h,'Position',pos);
        end
        
        % Close the window
        function close(obj)
            close(obj.h);
        end
    end
end