function when(func_name)
%% this function is designed to check the date when the input function was introduced by MATLAB
% The input variable can be either a string or a cell containing different
% strings. The information are gathered from web. Check the following examples:
%
% Example 1:
% >> when('rand')
% ## rand is a built-in function (Introduced before R2006a)
%
% Example 2:
% func_name = {'rand','plot','grid','findstr','weboptions'};
% when(func_name)
% ## rand is a built-in function (Introduced before R2006a)
% ## plot is a built-in function (Introduced before R2006a)
% ## grid is a Matlab function or an ordinary m-file (Introduced before R2006a)
% ## findstr is a built-in function (Introduced before R2006a)
% ## weboptions is a Matlab function or an ordinary m-file (Introduced in R2014b)
%
%
% -------------------------------------------------
% code by: Reza Ahmadzadeh (reza.ahmadzadeh@iit.it)
% -------------------------------------------------
%
if ischar(func_name)
    checkFunctionName(func_name);
elseif iscell(func_name)
    for ii = 1:size(func_name,2)
        fname = func_name{1,ii};
        checkFunctionName(fname);
    end
else
    disp('Error! The input should be either a string or a cell-string.');
end
end

function checkFunctionName(fname)
A = exist(fname);
switch A
    case 0
        w = ' does not exist';
    case 1
        w = ' is a variable in the workspace';
    case 2
        w = ' is a Matlab function or an ordinary m-file';
    case 3
        w = ' is a MEX-file.';
    case 4
        w = ' is a Simulink model or library';
    case 5
        w = ' is a built-in function';
    case 6
        w = ' is a P-file.';
    case 7
        w = ' is a folder';
    case 8
        w = ' is a class';
    otherwise
        w = ' is nonsense';
end

if A == 2 || A == 5
    url = ['http://mathworks.com/help/matlab/ref/' fname '.html'];
    [str,status] = urlread(url);
    if isempty(str)
        url = ['http://mathworks.com/help/simulink/slref/' fname '.html'];
        [str,status] = urlread(url);
    end
    if status == 0
        disp(['Connection error or no online documentation found for [' fname '].']);
        return;
    end
    idx = strfind(str,'release_introduced');
    if isempty(idx)
        str3 = 'No information found';
    else
        str1 = str(1,idx+24:idx+47);
        idx2 = strfind(str1,'R20');
        str3 = str1(1,1:idx2+5);
    end
    strr = ['## [' fname ']' w ' (' str3 ').'];
else
    strr = ['## [' fname ']' w '.'];
end
disp(strr);
end



