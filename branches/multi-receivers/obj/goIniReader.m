% =========================================================================
%   OBJECT goIniReader -> for goGPS 
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
% LIST of METHODS
%
% IniReader(fileName, verbosity)                creator
% delete(obj)                                   distructor (empty)
%
% FILE ----------------------------------------------------------------
%
%  setFileName(obj, fileName)                   set the complete path of the ini file
%  fileName = getFileName(obj)                  get the complete path of the ini file
%
% INI SYNTAX ----------------------------------------------------------
%
%  setCommentChar(obj, character)               define the characters that define comment lines
%
% R/W MODES OF THE FILE -----------------------------------------------
%
%  setRW(obj,rwMode)                            set mode for open the ini file
%  rwMode = getRW(obj)                          get the current mode of the file
%  bool = readMode(obj)                         is it in read mode?
%  bool = writeMode(obj)                        is it in write mode?
%
%  bool = getReadStatus(obj)                    return if the file has been already parsed
%
%  readFile(obj)                                force reading of the File
%  update(obj, filename, force)                 update the object when needed:
%                                                - filename changed 
%                                                - force flag == 1
%                                                - INI not yet read
%
% MANAGING LOGGING ----------------------------------------------------
%
%  setColorMode(obj, bool)                      set useage of colors in text output
%  bool = getColorMode(obj)                     get useage of colors in text output
%  setVerbosityLev(obj, verbosity)              set level of verbosity
%  verbosity = getVerbosityLev(obj)             get level of verbosity
%
% GETTER OF THE PARAMETERS --------------------------------------------
%  
%  isS = isSection(obj, section)                get the presence of a section 
%  isK = isKey(obj, section, key)               get the presence of a key 
%  sectionList = getSections(obj)               get the list of available sections
%  keyList = getKeys(obj, <section>)            get the list of the keys available
%  data = getData(obj, <section>, key)          get the value of a specified key
%
% MODIFIER FUNCTIONS --------------------------------------------------
%
%  addSection(obj, newSection)                  add a new section to the object
%  addKey(obj, section, key, data)              add a new keys to the object
%  rmSection(obj, section)                      remove a section from the object IniReader
%  rmKey(obj, section, key)                     remove a key from the object IniReader
%  editKey(obj, section, key, data)             edit a key in the object IniReader
%
% CONTAIN FUNCTION ----------------------------------------------------
%
% isPresent = containsSection(obj, section)     search a section in the object IniReader
% isPresent = containsKey(obj, section, key)    search a key in the object IniReader
%
% VISUALIZATION FUNCTIONS ---------------------------------------------
%
%  listSections(obj, <colorMode>)               list the list of available sections
%  listKeys(obj, <section>, <colorMode>)        list the list of the keys available
%  showData(obj, <section>, <colorMode>)        list the data contained in the ini
%
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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
        snrThr = 0;
    end
    
    methods (Access = 'public')

        % Get the number of receiver in the INI
        function nR = getNumRec(obj)
            nR = obj.getData('Receivers','nRec');
        end
        
        % Get the minimum angle of acceptance for a satellit
        function cutoff = getCutoff(obj)
            cutoff = ini.getData('Generic','cutoff');
            if (isempty(cutoff))
                cutoff = obj.cutoff;
            end            
        end
        
        % Get the minimum SNR threshold acceptable
        function snrThr = getSnrThr(obj)
            snrThr = ini.getData('Generic','snrThr');
            if (isempty(snrThr))
                snrThr = obj.snrThr;
            end            
        end
        
    end
    
end
