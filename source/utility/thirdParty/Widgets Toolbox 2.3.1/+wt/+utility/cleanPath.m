function path = cleanPath(path)
% cleanPath - Utility to clean and standardize a file/folder path
%
% This function will clean and standardize a file or folder path. It
% removes leading/trailing whitespace and removes any file separator from
% the trailing end of the path.
%
% Syntax:
%       path = wt.utility.cleanPath(path)
%
% Inputs:
%       path - the path to a file or folder
%
% Outputs:
%       path - the cleaned path
%
% Examples:
%
%     >> path = "   C:\Program Files\MATLAB\" %Note leading space and trailing separator
%     >> path = wt.utility.cleanPath(path)
%
%     path =
%
%         "C:\Program Files\MATLAB"
%

% Copyright 2020-2021 The MathWorks Inc.
% ---------------------------------------------------------------------

% File separator - in case of regional variants
fsep = regexptranslate("escape",filesep);

% Pattern for regexp
fsepOpts = join(["\\","/",fsep],"|");
pattern = "^\s+|(?<=\S)(\s|" + fsepOpts + ")+$";

% Perform replacement
path = regexprep(path,pattern,"");
