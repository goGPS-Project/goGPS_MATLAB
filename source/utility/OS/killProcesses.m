% killProcesses kills processes that match given keywords and, optionally, have exceeded a CPU time limit.
%
% INPUT
%   keywords      - cell array of strings representing keywords to filter the processes by
%   max_time      - (optional) maximum CPU time in seconds; processes exceeding this will be killed
%
% OUTPUT
%   None
%
% This function uses the getProcesses function to retrieve a list of processes
% matching the specified keywords. If max_time is provided, only processes with a CPU time
% greater than max_time will be targeted for termination.
%
% SYNTAX
%   killProcesses(keywords)
%   killProcesses(keywords, max_time)
%
% EXAMPLE USAGE
%   killProcesses({'notepad', 'matlab'})
%   This will kill all processes that have 'notepad' or 'matlab' in their command line or path.
%
%   killProcesses({'python'}, 3600)
%   This will kill all 'python' processes running longer than 3600 seconds (1 hour).
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

function killProcesses(keywords, max_time)

    % Check if the max_time argument is provided
    if nargin < 2
        max_time = 0; % If not provided, set it to 0
    end
    
    % Retrieve the list of running processes matching the keywords
    processes = getProcesses(keywords);
    
    % Loop through the processes and kill those that match the criteria
    for i = 1:length(processes)
        process = processes(i);
        if ~isempty(process) && process.cpu_time > max_time
            % Kill the process
            if ispc
                killCommand = sprintf('taskkill /PID %d /F', process.pid);
            elseif isunix
                killCommand = sprintf('kill -9 %d', process.pid);
            else
                error('Unsupported operating system.');
            end
    
            % Execute the kill command
            [status, cmdout] = system(killCommand);
            if status == 0
                fprintf('Process with PID %d killed successfully.\n', process.pid);
            else
                fprintf('Failed to kill process with PID %d: %s\n', process.pid, cmdout);
            end
        end
    end
end
