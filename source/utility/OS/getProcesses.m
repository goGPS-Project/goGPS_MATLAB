% getProcesses retrieves a list of running system processes that match a set of keywords
%
% INPUT
%   keywords      - cell array of strings representing keywords to filter the processes by
%
% OUTPUT
%   processes     - structure array containing detailed information about each process that matches the keywords
%
% The output structure has the following fields:
%   name          - name of the executable
%   params        - parameters used to run the executable
%   full_path     - full path of the executable
%   user          - user under which the process is running
%   pid           - process ID
%   mem           - memory usage of the process in MB
%   cpu           - CPU usage of the process in percentage (not available on Windows)
%   cpu_time      - CPU time used by the process in seconds
%
% The function supports both Windows and Unix-based systems (macOS, Linux). It filters processes
% to include only those that contain all the specified keywords in their command line or path.
%
% Note: On Windows, the 'cpu' field is not available, and 'user' is not provided by the WMIC command used.
%
% SYNTAX
%   processes = getProcesses(keywords)
%
% EXAMPLE USAGE
%   processes = getProcesses({'matlab', 'python'})
%   This will return processes that have both 'matlab' and 'python' in their command line or path.
%

%--------------------------------------------------------------------------
%   ______   ______ _______ _    _ _______ (TM)
%   |_____] |_____/ |______  \  /  |_____|
%   |_____] |    \_ |______   \/   |     |
%                                           v0.9999.1
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti
%
%  This software is based on GReD's goGPS software
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%-------------------------------------------------------------------------------
%
%    GReD Geomatics Research & Development - Proprietary Use License
%    Version 1.0, - July, 2023
%
%    This software and its associated code, functions, and tools are the property
%    of GReD (Geomatics Research & Development). All rights reserved.
%
%    Redistribution and use of this software in source and binary forms, with or
%    without modification, are strictly prohibited without explicit permission
%    from GReD, except as permitted by law.
%
%    You may not distribute, sublicense, sell, or otherwise transfer any part of
%    this software, code, functions, or tools, nor grant usage rights to any
%    third party without explicit written permission from GReD.
%
%    Unauthorized use, reproduction, or distribution of any part of this software
%    will be considered a breach of the GReD Geomatics Research & Development
%    Proprietary Use License and may result in legal consequences.
%
%    For licensing inquiries or permissions, please contact GReD at:
%    info@g-red.eu
%
%    For the avoidance of doubt, any part of the software not explicitly
%    mentioned under this license remains under the terms of its
%    respective licence
%
%-------------------------------------------------------------------------------
% 1000010 1010010 1000101 1010110 1000001
%--------------------------------------------------------------------------

function processes = getProcesses(keywords)    
    if nargin < 1
        error('You must provide at least one keyword to filter by.');
    end
    
    if ~iscell(keywords)
        keywords = {keywords};
    end
    
    if ispc
        % Windows systems
        command = 'tasklist /v /fo csv';
        for i = 1:length(keywords)
            command = sprintf('%s | findstr "%s"', command, keywords{i});
        end
    elseif isunix
        % Unix-based systems (includes macOS and Linux)
        command = 'ps aux';
        for i = 1:length(keywords)
            command = sprintf('%s | grep "%s"', command, keywords{i});
        end
        % Exclude the grep command itself from the results
        command = [command ' | grep -v grep'];
    else
        error('Unsupported operating system.');
    end
    
    [status, cmdout] = system(command);
    
    %if status ~= 0
    %    error('Failed to retrieve the list of processes.');
    %end
    
    % Parse the command output and populate the structure
    processes = parseProcessesOutput(cmdout, ispc);
end

function processes = parseProcessesOutput(cmdout, is_windows)
    % Initialize an empty structure array
    processes = struct('name', {}, 'params', {}, 'full_path', {}, 'user', {}, ...
                       'pid', {}, 'mem', {}, 'cpu', {}, 'cpu_time', {});
    
    % Split the output into lines
    lines = regexp(cmdout, '\r?\n', 'split');
    lines = lines(~cellfun('isempty',lines)); % Remove empty lines
    
    for i = 1:length(lines)
        if is_windows
            % Parse CSV formatted line from Windows WMIC output
            tokens = regexp(lines{i}, ',', 'split');
            if length(tokens) < 5
                continue; % Skip if the line doesn't have enough fields
            end
            % Populate the structure
            processes(i).full_path = tokens{2};
            processes(i).params = ''; % WMIC does not provide params separately
            processes(i).name = getFileNameFromPath(tokens{2});
            processes(i).pid = str2double(tokens{4});
            processes(i).mem = str2double(tokens{5}) / 1024^2; % Convert bytes to MB
            processes(i).cpu_time = str2double(tokens{3}) / 1000; % Convert milliseconds to seconds
            processes(i).user = ''; % WMIC does not provide user in this format
            processes(i).cpu = NaN; % WMIC does not provide CPU usage in this format
        else
            % Parse line from Unix 'ps' output
            tokens = regexp(lines{i}, '\s+', 'split');
            if length(tokens) < 11
                continue; % Skip if the line doesn't have enough fields
            end
            % Populate the structure
            processes(i).user = tokens{1};
            processes(i).pid = str2double(tokens{2});
            processes(i).cpu = str2double(tokens{3});
            processes(i).mem = str2double(tokens{4}) * 1024; % Convert KB to MB
            processes(i).cpu_time = convertCpuTimeToSeconds(tokens{10});
            [processes(i).full_path, processes(i).params] = parseUnixCommandLine(tokens(11:end));
            processes(i).name = getFileNameFromPath(processes(i).full_path);
        end
    end
end

function [full_path, params] = parseUnixCommandLine(command_line_tokens)
    % This function attempts to separate the full path and parameters from a Unix command line.
    % It handles spaces in paths by checking for the existence of the file.
    full_path = '';
    params = '';
    for i = 1:length(command_line_tokens)
        potential_path = strjoin(command_line_tokens(1:i), ' ');
        if exist(potential_path, 'file') == 2
            full_path = potential_path;
            params = strjoin(command_line_tokens(i+1:end), ' ');
            break;
        end
    end
    if isempty(full_path)
        full_path = command_line_tokens{1}; % Fallback to the first token
        params = strjoin(command_line_tokens(2:end), ' ');
    end
end

function file_name = getFileNameFromPath(path)
    % Extract the file name from a full path.
    [~, name, ext] = fileparts(path);
    file_name = [name, ext];
end

function seconds = convertCpuTimeToSeconds(cpu_time)
    % Convert CPU time format to seconds
    time_parts = regexp(cpu_time, ':', 'split');
    seconds = 0;
    if length(time_parts) == 3
        seconds = str2double(time_parts{1}) * 3600 + ... % hours to seconds
                  str2double(time_parts{2}) * 60 + ...   % minutes to seconds
                  str2double(time_parts{3});             % seconds
    elseif length(time_parts) == 2
        seconds = str2double(time_parts{1}) * 60 + ...   % minutes to seconds
                  str2double(time_parts{2});             % seconds
    else
        seconds = str2double(cpu_time);                  % seconds
    end
end