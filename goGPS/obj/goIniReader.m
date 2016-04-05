% =========================================================================
%   OBJECT goIniReader => for goGPS 
%   father: IniReader
% =========================================================================
%
% DESCRIPTION:
%   Object to read an ini file
%
% INI EXAMPLE:
%   [Section 1]
%       array = [1 2 3]
%       string = "I'm a string"
%       stringCellArray = ["I" "am" "a string cell" "array"]
%   ; I'm a commen
%   # I'm also a comment
%   [Section 2]
%       number = 10
%   [Section n]     # comment
%       double = 10.4
%
% REQUIRES:
%   cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
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
classdef goIniReader < IniReader
    
    properties (GetAccess = 'private', SetAccess = 'private')
        cutoff = 10;
        csThr = 3;
        snrThr = 0;
    end
    
    methods
        function obj = goIniReader(fileName, verbosity)
            % INI creator hinerited by iniReader
            if (nargin < 1)
                fileName = '';
            end
            if isempty(fileName)
                fileName = '';
            end
            obj.setFileName(fileName);
            if (nargin == 2)
                obj.setVerbosityLev(verbosity);
            end
        end
    end
    
    methods (Access = 'public')

        
        function nR = getNumRec(obj)
            % Get the number of receiver in the INI
            nR = obj.getData('Receivers','nRec');
        end
        
        function rate = getCaptureRate(obj)
            % Get the sampling rate (Hz) to be used in monitor mode
            rate = obj.getData('Monitor', 'rate');
            if isempty(rate)
                rate = 1;
            end
        end
        
        function setCaptureRate(obj, rate)
            % Set the sampling rate (Hz) to be used in monitor mode
            rate = str2double(rate(1:2));
            obj.addData('Monitor', 'rate', rate);
        end
                
        function [dataPath fileName nR] = getRecFiles(obj)
            % Get the dataPath containing the files, file names and number of available receivers
            nR = getNumRec(obj);
            
            dataPath = obj.getData('Receivers','data_path');
            if (isempty(dataPath))
                dataPath = '';
            end
            
            fileName = obj.getData('Receivers','file_name');
            if (isempty(fileName))
                fileName = '';
            else
                % If missing, get the number of receiver from the file name list
                if isempty(nR)
                    nR = length(fileName);
                end
            end
            
            % If I have only one receiver fileName is not a cell => convert it
            if nR == 1
                tmp = fileName;
                fileName = cell(1);
                fileName{1} = tmp;
            end
        end
        
        function [geometry ev_point] = getGeometry(obj)
            % Get the receiver coordinates in instrumental RF
            geometry = zeros(3,obj.getNumRec());
            for r = 1:obj.getNumRec()
               geometry(:,r) = obj.getData('Antennas RF',['XYZ_ant' num2str(r)]);
            end
            ev_point = obj.getData('Antennas RF','XYZ_ev_point');
        end
        
        function cutoff = getCutoff(obj)
            % Get the minimum angle of acceptance for a satellite
            cutoff = obj.getData('Generic','cutoff');
            if (isempty(cutoff))
                cutoff = obj.cutoff;
            end            
        end
        
        function csThr = getCsThr(obj)
            % Get the threshold to identify a cycle slip
            csThr = obj.getData('Generic','csThr');
            if (isempty(csThr))
                csThr = obj.csThr;
            end            
        end
                
        function snrThr = getSnrThr(obj)
            % Get the minimum SNR threshold acceptable
            snrThr = obj.getData('Generic','snrThr');
            if (isempty(snrThr))
                snrThr = obj.snrThr;
            end            
        end
        
        function timeStep = getTimeStep(obj)
            % Get time increment (e.g. v = obs(t+timeStep)-obs(t)
            timeStep = obj.getData('Variometric','timeStep');
            if (isempty(timeStep))
                timeStep = 1;
            end    
        end
        
    end
    
end
