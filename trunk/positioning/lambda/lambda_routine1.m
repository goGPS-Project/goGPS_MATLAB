function [afixed,sqnorm,Qahat,Z] = lambda_routine1 (afloat,Qahat,ncands,fidlog)
%LAMBDA1: Integer ambiguity estimation using LAMBDA (extended version)
%
% This routine performs an integer ambiguity estimation using the 
% LAMBDA-method, as developed by the Delft University of Technology, 
% Mathematical Geodesy and Positioning.
%
% This version (lambda1) is intended to be used mainly for research.
% It has more options than the original code, such as output of
% intermediate results, optional number of candidates to be returned
% and so on. If you are not interested in these options, you can use 
% "lambda2"
%
% Input arguments:
%    afloat: Float ambiguities (\hat{a}, must be a column!)
%    Qahat : Variance/covariance matrix of ambiguities (Q_\hat{a})
%    ncands: Number of candidates to be returned (default 2)
%    fidlog: File-id for intermediate results (default 0, no output)
%
% Output arguments:
%    afixed: Estimated integers
%    sqnorm: Distance between candidates and float ambiguity vector
%    Qzhat : Decorrelated variance/covariance matrix
%    Z     : Transformation matrix

% ----------------------------------------------------------------------
% File.....: lambda1.m
% Date.....: 19-MAY-1999
% Version..: 2.0b
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
% Initialise
% ----------------------------------------------------------------------

if nargin < 3; ncands = 2; end;
if nargin < 4; fidlog = 0; end;

n      = size (Qahat,1);
afixed = zeros(n,ncands);
sqnorm = zeros(1,ncands);

if fidlog ~= 0;
   if fidlog == 1; clc; end;
   fprintf (fidlog, '\n');
   fprintf (fidlog, ' -----------------------------------------------------------------------\n');
   fprintf (fidlog, ' "LAMBDA" - Least Squares Ambiguity Decorrelation Adjustment            \n');
   fprintf (fidlog, ' Rel. 2.0b d.d. 19 MAY 1999                                             \n');
   fprintf (fidlog, ' (c) Department of Mathematical Geodesy and Positioning (MGP)           \n');
   fprintf (fidlog, '     Faculty of Civil Engineering and Geosciences                       \n');
   fprintf (fidlog, '     Delft University of Technology, The Netherlands                    \n');
   fprintf (fidlog, ' -----------------------------------------------------------------------\n');
   fprintf (fidlog, '\n');
end;

% ----------------------------------------------------------------------
% Display the input-arguments
% ----------------------------------------------------------------------

writemat (fidlog,'Original float ambiguities',afloat);
writemat (fidlog,'Original variance/covariance matrix (Qahat)',Qahat);
writemat (fidlog,'Number of candidates requested',ncands);

% ----------------------------------------------------------------------
% --- Perform tests on the input, these tests are rather extensive   ---
% --- and can be switched off (by commenting-out) if necessary       ---
% --- Tests: Is the Q-matrix symmetric?                              ---
% ---        Is the Q-matrix positive-definite?                      ---
% ---        Do the Q-matrix and ambiguity-vector have identical     ---
% ---        dimensions?                                             ---
% ---        Is the ambiguity vector a column?                       ---
% ----------------------------------------------------------------------

if ~isequal(Qahat-Qahat'<1E-12,ones(size(Qahat)));
  error ('Variance/covariance matrix is not symmetric!');
end;

if sum(eig(Qahat)>0) ~= size(Qahat,1);
  error ('Variance/covariance matrix is not positive definite!');
end;

if length(afloat) ~= size(Qahat,1);
  error (['Variance/covariance matrix and vector of ambiguities do not have' ...
	  ' identical dimensions!']);
end;

if size(afloat,2) ~= 1;
  error ('Ambiguity-vector should be a column-vector');
end;

% ----------------------------------------------------------------------
% Make estimates in 'a' between -1 and +1 by subtracting an
% integer number, store the increments in 'incr' (= shift the centre
% of the ellipsoid over the grid by an integer translation)
% ----------------------------------------------------------------------

incr   = afloat - rem(afloat,1);
afloat = rem(afloat,1);
writemat (fidlog,'shifted ambiguities - increments',[afloat incr]);

% ----------------------------------------------------------------------
% Compute the Z-transformation based on L and D of Q, ambiguities
% are transformed according to \hat{z} = Z^T\hat{a}, Q is transformed
% according to Q_\hat{z} = Z^T * Q_\hat{a} * Z
% ----------------------------------------------------------------------

[Qz,Z,L,D,z] = decorrel (Qahat,afloat);

writemat (fidlog,'Decorrelated Q-matrix',Qz);
writemat (fidlog,'Transformation matrix (Z)',Z);
writemat (fidlog,'L-matrix, as used in ''search''',L);
writemat (fidlog,'D-matrix, as used in ''search''',D);
writemat (fidlog,'Transformed reduced ambiguities',z);

% ----------------------------------------------------------------------
% If the number of requested candidates is no more than "n+1", with "n"
% the dimension of Q_\hat(a), we can compute a suitable Chi^2 such that
% we have the requested number of candidates at minimum; use an eps
% to make sure the candidates are inside the ellipsoid and not exactly
% on the border.
%
% If the requested number of canditates is larger, we have to "guess"
% the volume of the search-space, it can not be guarenteed the requested
% number of canditates will actually be found.
% ----------------------------------------------------------------------

Chi2 = chistart (D,L,afloat,ncands);
writemat (fidlog,'Chi-squared, size of the search-ellipsoid',Chi2);

% ----------------------------------------------------------------------
% Find the requested number of candidates (search)
% ----------------------------------------------------------------------

[afixed,sqnorm,ierr] = lsearch (afloat,L,D,Chi2,ncands);

if ierr == 1;
   fprintf (1,'Warning: Not enough candidates were found!!\n');
   fprintf (1,'         Requested number.: %2d\n',ncands);
end;

writemat (fidlog,'Fixed transformed ''reduced'' ambiguities',afixed);
writemat (fidlog,'Corresponding squared norms',sqnorm);

% ----------------------------------------------------------------------
% Finally, perform the back-transformation and add the increments
% ----------------------------------------------------------------------

afixed = (afixed' * inv(Z))';
afixed = round(afixed + repmat(incr,1,ncands));
writemat (fidlog,'Final fixed ambiguities',afixed);

% ----------------------------------------------------------------------
% End of routine: lambda1
% ----------------------------------------------------------------------
