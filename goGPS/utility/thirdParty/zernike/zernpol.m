function z = zernpol(n,m,r,nflag)
%ZERNPOL Radial Zernike polynomials of order N and frequency M.
%   Z = ZERNPOL(N,M,R) returns the radial Zernike polynomials of 
%   order N and frequency M, evaluated at R.  N is a vector of
%   positive integers (including 0), and M is a vector with the
%   same number of elements as N.  Each element k of M must be a 
%   positive integer, with possible values M(k) = 0,2,4,...,N(k)
%   for N(k) even, and M(k) = 1,3,5,...,N(k) for N(k) odd.  R is 
%   a vector of numbers between 0 and 1.  The output Z is a matrix 
%   with one column for every (N,M) pair, and one row for every
%   element in R.
%
%   Z = ZERNPOL(N,M,R,'norm') returns the normalized Zernike poly-
%   nomials.  The normalization factor Nnm = sqrt(2*(n+1)) is
%   chosen so that the integral of (r * [Znm(r)]^2) from r=0 to 
%   r=1 is unity.  For the non-normalized polynomials, Znm(r=1)=1
%   for all [n,m].
%
%   The radial Zernike polynomials are the radial portion of the
%   Zernike functions, which are an orthogonal basis on the unit
%   circle.  The series representation of the radial Zernike 
%   polynomials is
%
%          (n-m)/2
%            __
%    m      \       s                                          n-2s
%   Z(r) =  /__ (-1)  [(n-s)!/(s!((n-m)/2-s)!((n+m)/2-s)!)] * r
%    n      s=0
%
%   The following table shows the first 12 polynomials.
%
%       n    m    Zernike polynomial    Normalization
%       ---------------------------------------------
%       0    0    1                        sqrt(2)
%       1    1    r                           2
%       2    0    2*r^2 - 1                sqrt(6)
%       2    2    r^2                      sqrt(6)
%       3    1    3*r^3 - 2*r              sqrt(8)
%       3    3    r^3                      sqrt(8)
%       4    0    6*r^4 - 6*r^2 + 1        sqrt(10)
%       4    2    4*r^4 - 3*r^2            sqrt(10)
%       4    4    r^4                      sqrt(10)
%       5    1    10*r^5 - 12*r^3 + 3*r    sqrt(12)
%       5    3    5*r^5 - 4*r^3            sqrt(12)
%       5    5    r^5                      sqrt(12)
%       ---------------------------------------------
%
%   Example:
%
%       % Display three example Zernike radial polynomials
%       r = 0:0.01:1;
%       n = [3 2 5];
%       m = [1 2 1];
%       z = zernpol(n,m,r);
%       figure
%       plot(r,z)
%       grid on
%       legend('Z_3^1(r)','Z_2^2(r)','Z_5^1(r)','Location','NorthWest')
%
%   See also ZERNFUN, ZERNFUN2.

% A note on the algorithm.
% ------------------------
% The radial Zernike polynomials are computed using the series
% representation shown in the Help section above. For many special
% functions, direct evaluation using the series representation can
% produce poor numerical results (floating point errors), because
% the summation often involves computing small differences between
% large successive terms in the series. (In such cases, the functions
% are often evaluated using alternative methods such as recurrence
% relations: see the Legendre functions, for example). For the Zernike
% polynomials, however, this problem does not arise, because the
% polynomials are evaluated over the finite domain r = (0,1), and
% because the coefficients for a given polynomial are generally all
% of similar magnitude.
% 
% ZERNPOL has been written using a vectorized implementation: multiple
% Zernike polynomials can be computed (i.e., multiple sets of [N,M]
% values can be passed as inputs) for a vector of points R.  To achieve
% this vectorization most efficiently, the algorithm in ZERNPOL
% involves pre-determining all the powers p of R that are required to
% compute the outputs, and then compiling the {R^p} into a single
% matrix.  This avoids any redundant computation of the R^p, and 
% minimizes the sizes of certain intermediate variables.
%
%   Paul Fricker 11/13/2006


% Check and prepare the inputs:
% -----------------------------
if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
    error('zernpol:NMvectors','N and M must be vectors.')
end

if length(n)~=length(m)
    error('zernpol:NMlength','N and M must be the same length.')
end

n = n(:);
m = m(:);
length_n = length(n);

if any(mod(n-m,2))
    error('zernpol:NMmultiplesof2','All N and M must differ by multiples of 2 (including 0).')
end

if any(m<0)
    error('zernpol:Mpositive','All M must be positive.')
end

if any(m>n)
    error('zernpol:MlessthanN','Each M must be less than or equal to its corresponding N.')
end

if any( r>1 | r<0 )
    error('zernpol:Rlessthan1','All R must be between 0 and 1.')
end

if ~any(size(r)==1)
    error('zernpol:Rvector','R must be a vector.')
end

r = r(:);
length_r = length(r);

if nargin==4
    isnorm = ischar(nflag) & strcmpi(nflag,'norm');
    if ~isnorm
        error('zernpol:normalization','Unrecognized normalization flag.')
    end
else
    isnorm = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Zernike Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the required powers of r:
% -----------------------------------
rpowers = [];
for j = 1:length(n)
    rpowers = [rpowers m(j):2:n(j)];
end
rpowers = unique(rpowers);

% Pre-compute the values of r raised to the required powers,
% and compile them in a matrix:
% -----------------------------
if rpowers(1)==0
    rpowern = arrayfun(@(p)r.^p,rpowers(2:end),'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
    rpowern = [ones(length_r,1) rpowern];
else
    rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
end

% Compute the values of the polynomials:
% --------------------------------------
z = zeros(length_r,length_n);
for j = 1:length_n
    s = 0:(n(j)-m(j))/2;
    pows = n(j):-2:m(j);
    for k = length(s):-1:1
        p = (1-2*mod(s(k),2))* ...
                   prod(2:(n(j)-s(k)))/          ...
                   prod(2:s(k))/                 ...
                   prod(2:((n(j)-m(j))/2-s(k)))/ ...
                   prod(2:((n(j)+m(j))/2-s(k)));
        idx = (pows(k)==rpowers);
        z(:,j) = z(:,j) + p*rpowern(:,idx);
    end
    
    if isnorm
        z(:,j) = z(:,j)*sqrt(2*(n(j)+1));
    end
end

% EOF zernpol