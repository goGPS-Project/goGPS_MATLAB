function z = zernfun(n,m,r,theta,nflag)
%ZERNFUN Zernike functions of order N and frequency M on the unit circle.
%   Z = ZERNFUN(N,M,R,THETA) returns the Zernike functions of order N
%   and angular frequency M, evaluated at positions (R,THETA) on the
%   unit circle.  N is a vector of positive integers (including 0), and
%   M is a vector with the same number of elements as N.  Each element
%   k of M must be a positive integer, with possible values M(k) = -N(k)
%   to +N(k) in steps of 2.  R is a vector of numbers between 0 and 1,
%   and THETA is a vector of angles.  R and THETA must have the same
%   length.  The output Z is a matrix with one column for every (N,M)
%   pair, and one row for every (R,THETA) pair.
%
%   Z = ZERNFUN(N,M,R,THETA,'norm') returns the normalized Zernike
%   functions.  The normalization factor sqrt((2-delta(m,0))*(n+1)/pi),
%   with delta(m,0) the Kronecker delta, is chosen so that the integral
%   of (r * [Znm(r,theta)]^2) over the unit circle (from r=0 to r=1,
%   and theta=0 to theta=2*pi) is unity.  For the non-normalized
%   polynomials, max(Znm(r=1,theta))=1 for all [n,m].
%
%   The Zernike functions are an orthogonal basis on the unit circle.
%   They are used in disciplines such as astronomy, optics, and
%   optometry to describe functions on a circular domain.
%
%   The following table lists the first 15 Zernike functions.
%
%       n    m    Zernike function             Normalization
%       ----------------------------------------------------
%       0    0    1                              1/sqrt(pi)
%       1    1    r * cos(theta)                 2/sqrt(pi)
%       1   -1    r * sin(theta)                 2/sqrt(pi)
%       2    2    r^2 * cos(2*theta)             sqrt(6/pi)
%       2    0    (2*r^2 - 1)                    sqrt(3/pi)
%       2   -2    r^2 * sin(2*theta)             sqrt(6/pi)
%       3    3    r^3 * cos(3*theta)             sqrt(8/pi)
%       3    1    (3*r^3 - 2*r) * cos(theta)     sqrt(8/pi)
%       3   -1    (3*r^3 - 2*r) * sin(theta)     sqrt(8/pi)
%       3   -3    r^3 * sin(3*theta)             sqrt(8/pi)
%       4    4    r^4 * cos(4*theta)             sqrt(10/pi)
%       4    2    (4*r^4 - 3*r^2) * cos(2*theta) sqrt(10/pi)
%       4    0    6*r^4 - 6*r^2 + 1              sqrt(5/pi)
%       4   -2    (4*r^4 - 3*r^2) * sin(2*theta) sqrt(10/pi)
%       4   -4    r^4 * sin(4*theta)             sqrt(10/pi)
%       ----------------------------------------------------
%
%   Example 1:
%
%       % Display the Zernike function Z(n=5,m=1)
%       x = -1:0.01:1;
%       [X,Y] = meshgrid(x,x);
%       [theta,r] = cart2pol(X,Y);
%       idx = r<=1;
%       z = nan(size(X));
%       z(idx) = zernfun(5,1,r(idx),theta(idx));
%       figure
%       pcolor(x,x,z), shading interp
%       axis square, colorbar
%       title('Zernike function Z_5^1(r,\theta)')
%
%   Example 2:
%
%       % Display the first 10 Zernike functions
%       x = -1:0.01:1;
%       [X,Y] = meshgrid(x,x);
%       [theta,r] = cart2pol(X,Y);
%       idx = r<=1;
%       z = nan(size(X));
%       n = [0  1  1  2  2  2  3  3  3  3];
%       m = [0 -1  1 -2  0  2 -3 -1  1  3];
%       Nplot = [4 10 12 16 18 20 22 24 26 28];
%       y = zernfun(n,m,r(idx),theta(idx));
%       figure('Units','normalized')
%       for k = 1:10
%           z(idx) = y(:,k);
%           subplot(4,7,Nplot(k))
%           pcolor(x,x,z), shading interp
%           set(gca,'XTick',[],'YTick',[])
%           axis square
%           title(['Z_{' num2str(n(k)) '}^{' num2str(m(k)) '}'])
%       end
%
%   See also ZERNPOL, ZERNFUN2.

%   Paul Fricker 2/28/2012

% Check and prepare the inputs:
% -----------------------------
if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
    error('zernfun:NMvectors','N and M must be vectors.')
end

if length(n)~=length(m)
    error('zernfun:NMlength','N and M must be the same length.')
end

n = n(:);
m = m(:);
if any(mod(n-m,2))
    error('zernfun:NMmultiplesof2', ...
          'All N and M must differ by multiples of 2 (including 0).')
end

if any(m>n)
    error('zernfun:MlessthanN', ...
          'Each M must be less than or equal to its corresponding N.')
end

if any( r>1 | r<0 )
    error('zernfun:Rlessthan1','All R must be between 0 and 1.')
end

if ( ~any(size(r)==1) ) || ( ~any(size(theta)==1) )
    error('zernfun:RTHvector','R and THETA must be vectors.')
end

r = r(:);
theta = theta(:);
length_r = length(r);
if length_r~=length(theta)
    error('zernfun:RTHlength', ...
          'The number of R- and THETA-values must be equal.')
end

% Check normalization:
% --------------------
if nargin==5 && ischar(nflag)
    isnorm = strcmpi(nflag,'norm');
    if ~isnorm
        error('zernfun:normalization','Unrecognized normalization flag.')
    end
else
    isnorm = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Zernike Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the required powers of r:
% -----------------------------------
m_abs = abs(m);
rpowers = [];
for j = 1:length(n)
    rpowers = [rpowers m_abs(j):2:n(j)];
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
z = zeros(length_r,length(n));
for j = 1:length(n)
    s = 0:(n(j)-m_abs(j))/2;
    pows = n(j):-2:m_abs(j);
    for k = length(s):-1:1
        p = (1-2*mod(s(k),2))* ...
                   prod(2:(n(j)-s(k)))/              ...
                   prod(2:s(k))/                     ...
                   prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                   prod(2:((n(j)+m_abs(j))/2-s(k)));
        idx = (pows(k)==rpowers);
        z(:,j) = z(:,j) + p*rpowern(:,idx);
    end
    
    if isnorm
        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1)/pi);
    end
end
% END: Compute the Zernike Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the Zernike functions:
% ------------------------------
idx_pos = m>0;
idx_neg = m<0;

if any(idx_pos)
    z(:,idx_pos) = z(:,idx_pos).*cos(theta*m_abs(idx_pos)');
end
if any(idx_neg)
    z(:,idx_neg) = z(:,idx_neg).*sin(theta*m_abs(idx_neg)');
end

% EOF zernfun