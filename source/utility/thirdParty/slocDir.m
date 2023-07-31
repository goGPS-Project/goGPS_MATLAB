%SLOCDIR Counts number source lines of code.
%   LOCSUM = SLOCDIR(VARARGIN) returns the line count for directory tree.  
%   If there are multiple functions in one file, subfunctions are not 
%   counted separately, but rather together.
%
%   This function calls SLOC to count lines of code on an individual file.
%   SLOC was obtained from the File Exchange on MATLAB Central, and it
%   was written by Raymond Norris of MathWorks Inc.
%
%   The following statistics are calculated for a directory tree
%       - sum of lines of code from all mfiles in the directory tree
%       - sum of lines of code from all mfiles in parent directory 
%       - sum of lines of code from all mfiles in subdirectories
%       - lines of code in each mfile in directory tree
%
%   The output of this function can either be to a text file or the 
%   MATLAB command prompt or both.  
%
%   LOCSUM = SLOCDIR(DIRECTORY)
%   Searches directory tree for and outputs results to command line
%
%   LOCSUM = SLOCDIR(DIRECTORY, OUTPUTFILE)
%   DIRECTORY and it's subfolders are searched.  The results are output 
%   to the text file, OUTPUTFILE, and to the command line.
%
%   LOCSUM = SLOCDIR(DIRECTORY, OUTPUTFILE, CMDWINDOW_FLAG)
%   DIRECTORY and it's subfolders are searched.  The results are output 
%   to the text file, OUTPUTFILE.  If CMDWINDOW_FLAG is true, the 
%   results are also output to the command line.
%
%
%   Examples:
%   ========
%       locsum = slocdir(�c:\work�);
%       locsum = slocDir('C:\work', 'workspace_SLOC.txt');
%       locsum = slocDir('C:\work', 'workspace_SLOC.txt', false);
%
%
%   David Roberts
%   MITRE Corporation/CAAASD 2009
%   Revision: 1.0 Date: 2009/04/15 
function locSum = slocDir(varargin)
    global subDirectoryCellArr
    global fid
    global fileOutput
    global commandLineOutput
    global MAXNUMOFSUBDIRS
    
    MAXNUMOFSUBDIRS = 1000;
    
    directoryStr =    varargin{1};
    subDirectoryCellArr = cell(1,MAXNUMOFSUBDIRS);
    currCellIndex       = 0;
    fileOutput          = 0;
    commandLineOutput   = 1;
    
    if nargin > 1
        fid = fopen(varargin{2},'w');
        if (fid ==-1)
            fileOutput = 0;
        else
            fileOutput = 1;
            if nargin > 2
                commandLineOutput = varargin{3};
            end
        end        
    end
    locSum = slocSubDir(directoryStr);
    
    outputLocStats(directoryStr, locSum);
    
    % fill subDirectoryCellArr with all paths of subdirectories    
    currCellIndex = walkIn(directoryStr, currCellIndex);
    
    % search subdirectories and count LOC of mfiles
    for ii = 1:currCellIndex
        loc     = slocSubDir(subDirectoryCellArr{ii});   
        locSum  = loc + locSum;
    end
    
    % output total LOC for all mfile contained in directoryStr and its subfolders
    outputLocStats(directoryStr, locSum);
end
function currCellIndex = walkIn(root, currCellIndex)    
    
    global subDirectoryCellArr
       
    if root(end) ~= filesep
        root(end+1) = filesep;
    end
    
    dirListing  = dir([root '*']);    
    
    for ii=1:length(dirListing)
        if dirListing(ii).isdir
            if dirListing(ii).name(1) ~= '.'
                % Add this folder path to listing of subdirectories
                currCellIndex                       = currCellIndex +1;
                subDirectoryCellArr{currCellIndex}  = [root dirListing(ii).name];
            
                % Look in this folder for more subdirectories
                currCellIndex = walkIn(subDirectoryCellArr{currCellIndex}, currCellIndex);                        
            end
        end
    end
end
function loc = slocSubDir(directoryStr)
    loc     = 0;            
    
    if directoryStr(end) ~= filesep
        directoryStr(end+1) = filesep;
    end
    
    mfiles  = dir([directoryStr '*.m']);    
    for ii = 1:length(mfiles)
        pathname    = [directoryStr mfiles(ii).name];
        locFile     = sloc(pathname);
        loc         = loc + locFile;
        
        outputLocStats(pathname, locFile);
    end    
end
function outputLocStats(pathname, locSum)
    global fid
    global fileOutput
    global commandLineOutput
    if fileOutput
        fprintf(fid,'%s,%d \n',pathname, locSum);
    end    
    if commandLineOutput
        disp(sprintf('%s : %d LOC',pathname, locSum));        
    end
end
