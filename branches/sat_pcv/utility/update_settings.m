function [out] = update_settings(settings_dir_path, field, value)

% SYNTAX:
%   [out] = update_settings(settings_dir_path, field, value);
%
% INPUT:
%   settings_dir_path = path to settings folder
%   field = name of the field to be added
%   value = default value for the field (-1 to remove the field; new field 
%           name to rename the field)
%
% OUTPUT:
%   out = update outcome
%
% DESCRIPTION:
%   Utility to update goGPS settings file.

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

out = 0;

%directory containing settings files
settings_dir = dir(settings_dir_path);

%check the number of files contained in the directory
nmax = size(settings_dir,1);

j = 0;
for i = 1 : nmax

    %read the name of the i-th file
    got = getfield(settings_dir,{i,1},'name');

    %get the number of characters in the filename
    fn_length = size(got,2);

    %check that the filename has the ".mat" extension
    if (fn_length >= 4 & strcmp(got(fn_length - 3 : fn_length), '.mat'))
        
        j = j+1;
        %load the settings file
        load([settings_dir_path '/' got]);
        %check if settings were loaded
        try
            if isstruct(state)
                if ischar(value)
                    %rename "field" by using "value"
                    state.(value) = state.(field);
                    state = rmfield(state,field);
                else
                    if(value ~= -1)
                        %add the new field to 'state' struct
                        % (if it does not already exist)
                        if (~isfield(state,field))
                            state.(field) = value;
                        end
                    else
                        %remove the specified field from the 'state' struct
                        state = rmfield(state,field);
                    end
                end
                %save the new state in settings file
                save([settings_dir_path '/' got], 'state');
            end
        catch
        end
    end
end

% if (j == 0)
%     error(['No settings file found in folder ' settings_dir_path '.']);
% else
%     fprintf('%d settings files processed.\n', j);
% end

if (j ~= 0)
    out = 1;
end
