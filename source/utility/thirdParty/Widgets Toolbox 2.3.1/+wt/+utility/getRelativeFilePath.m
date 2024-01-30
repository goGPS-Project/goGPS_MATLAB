function RelPath = getRelativeFilePath(FullPath, RootPath, FlagRequireSubdir)
% getRelativeFilePath - Utility to return a relative file path within a root
% 
% This function will find the relative file path, given a full absolute
% path and a root path
%
% Syntax:
%       RelPath = wt.utility.getRelativeFilePath(FullPath, RootPath)
%
% Inputs:
%       FullPath - the full absolute path to a file or folder
%       RootPath - the root folder to get a relative path for
%       FlagRequireSubdir - optional flag indicating whether FullPath must
%           be a subdirectory of RootPath [(true)|false]
%
% Outputs:
%       RelPath - the relative path
%
% Examples:
%
%     >> FullPath = 'C:\Program Files\MATLAB\R2016b\toolbox'
%     >> RootPath = 'C:\Program Files\MATLAB'
%     >> RelPath = wt.utility.getRelativeFilePath(FullPath, RootPath)
% 
%     RelPath =
%          \R2016b\toolbox
% 
% Notes:
%   If FullPath is not a subdirectory of RootPath and FlagRequireSubdir is
%   false, the path will contain parent directory separators "..\"
%

% Copyright 2016-2020 The MathWorks Inc.
% ---------------------------------------------------------------------


% Validate inputs
arguments
    FullPath (1,1) string
    RootPath (1,1) string
    FlagRequireSubdir (1,1) logical = true
end


% Is RootPath empty?
if ~strlength(FullPath)
    
    RelPath = RootPath;
    
elseif ~strlength(RootPath)
    
    RelPath = FullPath;
    
else
    
    % Remove trailing filesep
    if endsWith(RootPath,filesep)
        RootPath{1}(end) = '';
    end
    if endsWith(FullPath,filesep)
        FullPath{1}(end) = '';
    end
    
    % Split the paths apart
    RootParts = strsplit(RootPath,filesep);
    FullParts = strsplit(FullPath,filesep);
    
    % Find where the paths diverge
    idx = 1;
    SmallestPath = min(numel(RootParts), numel(FullParts));
    while idx<=SmallestPath && strcmpi(RootParts(idx),FullParts(idx))
        idx = idx+1;
    end
    
    % Is the specified path outside of the root directory?
    NumAbove = max(numel(RootParts) - idx + 1, 0);
    if FlagRequireSubdir && NumAbove>0
        error('The specified path:\n\t"%s"\nis not a subdirectory of the root path:\n\t"%s"',...
            char(FullPath), char(RootPath) );
    else
        % In case full path is above the RootPath, add ".." paths
        ParentPaths = string(repmat(['..' filesep],1,NumAbove));
        
        % Form the relative path
        RelPath = filesep + fullfile(ParentPaths, FullParts{idx:end});
    end
    
    % What if paths are still the same?
    if isempty(RelPath)
        RelPath = "." + filesep;
    end
    
end %if isempty(RootPath)

