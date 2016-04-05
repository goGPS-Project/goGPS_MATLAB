function [afixed,sqnorm,Qahat,Z,D,L] = lambda_routine2 (afloat,Qahat)
%LAMBDA2: Integer ambiguity estimation using LAMBDA (basic version)
%
% This routine performs an integer ambiguity estimation using the 
% LAMBDA-method, as developed by the Delft University of Technology, 
% Mathematical Geodesy and Positioning.
%
% Input arguments:
%    afloat: Float ambiguities (\hat{a}, must be a column!)
%    Qahat : Variance/covariance matrix of ambiguities (Q_\hat{a})
%
% Output arguments:
%    afixed: Estimated integers
%    sqnorm: Distance between candidates and float ambiguity vector
%    Qahat : Decorrelated variance/covariance matrix
%    Z     : Transformation matrix

% ----------------------------------------------------------------------
% File.....: lambda2.m
% Date.....: 19-MAY-1999
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
% Initalise
% ----------------------------------------------------------------------

ncands = 2;
n      = size (Qahat,1);
afixed = zeros(n,ncands);
sqnorm = zeros(1,ncands);

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

% ----------------------------------------------------------------------
% Compute the Z-transformation based on L and D of Q, ambiguities
% are transformed according to \hat{z} = Z^T\hat{a}, Q is transformed
% according to Q_\hat{z} = Z^T * Q_\hat{a} * Z
% ----------------------------------------------------------------------

[Qahat,Z,L,D,afloat] = decorrel_v2 (Qahat,afloat);

% ----------------------------------------------------------------------
% Compute a suitable Chi^2 such that we have the requested number of 
% candidates at minimum; use an 'eps' to make sure the candidates are 
% inside the ellipsoid and not exactly on the border.
% ----------------------------------------------------------------------

Chi2 = chistart_v2 (D,L,afloat,ncands);

% ----------------------------------------------------------------------
% Find the requested number of candidates (search)
% ----------------------------------------------------------------------

[afixed,sqnorm,ierr] = lsearch_v2 (afloat,L,D,Chi2,ncands);
if ierr == 1; error ('Not enough candidates were found!!'); end;

% ----------------------------------------------------------------------
% --- Perform the back-transformation and add the increments
% ----------------------------------------------------------------------

afixed = (afixed' * (Z)^(-1))';
afixed = round(afixed + repmat(incr,1,ncands));

% ----------------------------------------------------------------------
% End of routine: lambda2
% ----------------------------------------------------------------------
