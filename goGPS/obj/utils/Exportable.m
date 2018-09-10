%   CLASS Exportable
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that to import / export objects
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 4 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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

classdef Exportable < handle

    methods (Access = public)
    end

    methods

    end

    methods (Abstract)
        %copy = getCopy(this);
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function out = toStruct(this)
            log = Core.getLogger();
            cm = log.getColorMode();
            log.setColorMode(false);
            warning off
            if (numel(this) == 0)
                out = [];
            end
            for i = 1 : numel(this)
                tmp = struct(this(i));
                if isfield(tmp, 'PreviousInstance__')
                    tmp = rmfield(tmp,'PreviousInstance__');
                end
                prp = fieldnames(tmp);
                out(i) = tmp;
                
                for p = 1 : numel(prp)
                    if isobject(out(i).(prp{p}))
                        % Manage goGPS singletons
                        try
                            out(i).(prp{p}) = this(i).(prp{p}).toStruct;
                        catch ex
                            log.addWarning(ex.message)
                        end
                    end
                end
            end
            warning on
            log.setColorMode(cm);
        end
        
        function importFromStruct(this, fields)
            log = Core.getLogger();

            prp = fieldnames(fields);
            for p = 1 : numel(prp)
                if isobject(this.(prp{p}))
                    try
                        this.(prp{p}) = this.(prp{p}).fromStruct(fields.(prp{p}));
                    catch ex
                        log.addWarning(ex.message)
                    end
                else
                    this.(prp{p}) = fields.(prp{p});
                end
            end
        end
    end

end
