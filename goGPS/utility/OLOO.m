function [imax, xfin, s2fin, ufin, Cxx, uout] = OLOO(A, y, Q)
% SYNTAX:
%   [imax, xfin, s2fin, ufin, Cxx, uout] = OLOO(A, y, Q)
%
% INPUT:
%   A: design matrix
%   y: observations vector
%   Q: cofactor matrix
%
% OUTPUT:
%   imax:  index of the rejected blocks
%   x_fin: estimated parameters without outlier
%   s2fin: a posteriori sigma without outlier
%   ufin:  estimated residuals without outlier
%   Cxx:   parameters covariance of the final solution
%   uout:  residual of outlier observation
%
% DESCRIPTION:
%   perform LS on blocks of correlated observations
%   identify one (block) outlier
%   reject it
%   re-estimate unknowns
%   according to the theory in "L. Biagi and S. Caldera. An efficient leave one block out approach to identify outliers.Journal of Applied Geodesy, Volume 7, Issue 1, pages 11..19, 2013"
%
%   this version is optimized to manage l.o.o. of 1 observation at time, no blocks!
%
% CREDITS:
%   1.0: Stefano Caldera, 22.05.2014
%   1.1: Stefano Caldera, Andrea Gatti 08.12.2016 ( speedup improvements )

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Stefano Caldera 22.05.2014
%  Contributors:     Stefano Caldera,
%                    Andrea Gatti 08.12.2016 ( speedup improvements )
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

% when FTABLE is undefined, redefine it
global FTABLE; if isempty(FTABLE); FTABLE = finv(0.9995,1,1:size(A,1))'; end
if (size(FTABLE,1) < size(A,1)); FTABLE = finv(0.9995,1,1:size(A,1))'; end
n_blocks = length(y);
[m, n] = size(A);           % m: number of observations, n: number of unknowns
uout = NaN;

if m - n > 0

    % compute the global solution
    if isdiag(Q) % Q is tipically diagonal -> let's use this information
        Qm = diag(Q);
        Qm = Qm ./ min(Qm);
        invQ = diag(1./Qm);
        At_invQ = A' * invQ;
        Ninv = inv(At_invQ * A);
        xcap = Ninv * At_invQ * y;
        um = y - A * xcap;
        s2cap = um' * invQ * um/(m-n);
    else
        Q = Q ./ (min(diag(Q)));
        At_invQ = A'/Q;
        Ninv = inv(At_invQ * A);
        xcap = Ninv * At_invQ * y;
        um = y - A * xcap;
        s2cap = um' / Q * um/(m-n);
        Qm = diag(Q);
    end
    % convert A in a sparse matrix, faster and lighter

    use_sparse_approach = false;

    if sum(A(1,:)==0) > size(A,2) / 0.577     % if A is sparse enough
        use_sparse_approach  = true;
        A = sparse(A);
    end
    if m - n > 1
        %% start outliers rejection
        Im = eye(n);
        Bm = Ninv * A';
        Cm = diag(A * Bm);
        Km = Qm - Cm;
        Kminv = Km.^-1;
        wm = Qm .* Kminv .* um;
        s2m = ((m-n) .* s2cap - um .* Kminv .* um) ./ (m - n - 1);

        % original loop
        if use_sparse_approach
            % modified loop to exploit the sparse property of the A matrix
            Qw = zeros(n_blocks,1);
            tmp = (A .* repmat(Kminv,1,n));
            for i = 1 : n_blocks
                idOk = (tmp(i,:) ~= 0);
                Qw(i) = Qm(i) + Bm(idOk,i)' * (Im(idOk,idOk) + tmp(i,idOk)' * Bm(idOk,i)')  * A(i,idOk)';
            end
        else
            Qw = zeros(n_blocks,1);
            for i = 1 : n_blocks
                Qw(i) = Qm(i) + Bm(:,i)' * (Im + A(i,:)' * Kminv(i) * Bm(:,i)') * A(i,:)';
                %Qw3=Bm(:,i)'*(Im+A(i,:)'*Kminv(i)*Bm(:,i)');
            end
        end

        %toc
        Qwinv = Qw.^-1;
        deg2 = m - n - 1;

        F = wm .* Qwinv .* wm ./ s2m;
        Flim = FTABLE(deg2);


        %% apply final solution
        % find maximum F(i)/Flim(i)
        [Fmax, imax] = max(abs(F ./ Flim));

        if (Fmax < 1)
            % no outlier
            imax = 0;
            xfin = xcap;
            s2fin = s2cap;
            Ninvfin = Ninv;
            ufin = um;
            Cxx = s2fin * Ninvfin;

        else
            % if the maximum ratio exceedes the threshold, the observation is eliminated from the solution
            uout = um(imax);
            xfin = xcap - Bm(:, imax) * Kminv(imax) * um(imax);
            yfin = y;
            yfin(imax) = [];

            Afin = A;
            Afin(imax,:) = [];

            s2fin = s2m(imax);

            Ninvfin = Ninv + Bm(:, imax) * Kminv(imax) * Bm(:, imax)';
            ufin = yfin - Afin * xfin;
            Cxx = s2fin * Ninvfin;

        end

    else
        imax = 0;
        xfin = xcap;
        s2fin = s2cap;
        Cxx = s2cap * Ninv;
        ufin = um;
    end

else
    % detection is not possibile
    imax = 0;
    xfin = [];
    Cxx = [];
    ufin = [];
    s2fin = NaN;
end

