function writemat (fid,text,matrix,format)
%WRITEMAT: Write the contents of a matrix to a file
%
% Input arguments:
%    fid:    File identifier:
%            0: Do not write at all
%            1: Write to screen (matlab default output)
%            n: Write to file 'n', which should be open
%    text  : Descriptive text
%    matrix: The matrix to be printed
%    format: Format to be used, optional, default = "%15.5f"

% ----------------------------------------------------------------------
% File.....: writemat.m
% Date.....: 19-MAR-1999
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

if nargin < 4; format = '%15.5f'; end;

if fid ~= 0;

   fprintf (fid,'%c',[text ':']);
   fprintf (fid,'\n\n');

   for i = 1:size(matrix,1);

      for j = 1:size(matrix,2);

         if size(format,1) < size (matrix,2);
            fprintf (fid,format(1,:),matrix(i,j));
         else;
            fprintf (fid,format(j,:),matrix(i,j));
         end;

      end;

      fprintf (fid,'\n');

   end;

   fprintf (fid,'\n');

end;

% ----------------------------------------------------------------------
% End of routine: writemat
% ----------------------------------------------------------------------
