function [y_splined, x_spline, s_weights, y_splined_ext] = splinerMat(x, y, dxs, reg_factor, x_ext)
% SYNTAX:
%   [y_splined, x_spline, s_weights, y_splined_ext] = splinerMat(x, y, dxs, reg_factor, <x_ext>)
%
% EXAMPLE:
%   [y_splined, x_spline, s_weights, y_splined_ext] = splinerMat(x, y, 4, 0, x_ext);
%   [y_splined, x_spline, s_weights, y_splined_ext] = splinerMat(x, [y y_var], 4, 0, x_ext);
%
% INPUT:
%   x            [n x 1] observation time
%   y            [n x 1] observation value
%                [n x 2] observation value and variances (in this case the variances of the data are taken into account)
%   dxs          [1 x 1] spline base size
%   reg_factor   [1 x 1] regularization factor on the first derivative
%   x_ext        [m x 1] points in which to compute the interpolation
%
% OUTPUT:
%   y_splined     [n x 1] interpolated observation (on the observation epochs)
%   x_spline      [o x 1] center of the (o) splines
%   s_weights     [o x 1] weights of the splines
%   y_splined_ext [m x 1] spline interpolated in x_ext positions
%
% DESCRIPTION:
%   Interpolate with cubic splines a given dataset
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro
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

    if (nargin < 4)
        reg_factor = 0;
    end
    if isempty(x)
        x = 1:size(y, 1);
    end
    
    inan = isnan(y);
    if (size(y,2) == 2)
        inan = inan(:,1) | inan(:,2);
    end
    
    x = double(x);
    y = double(y);

    x(inan) = [];
    y(inan,:) = [];
    
    [x, id] = sort(x);
    y = y(id,:);
    if ~isempty(y)
        if ((nargin == 3) || (reg_factor == 0))
            if (size(y,2) == 2)
                [y_splined, x_spline, s_weights] = spliner_v51(x,y(:,1),y(:,2),dxs);
            else
                [y_splined, x_spline, s_weights] = spliner_v5(x,y,dxs);
            end
        else
            if (size(y,2) == 2)
                [y_splined, x_spline, s_weights] = spliner_v51R(x,y(:,1),y(:,2),dxs, reg_factor);
            else
                [y_splined, x_spline, s_weights] = spliner_v5R(x,y,dxs, reg_factor);
            end
        end
    else
        y_splined = y;
        x_spline = [];
        s_weights = [];
    end

    % Interpolation => using spline to predict in different coordinates
    % (not present in the C version)
    if (nargin == 5)
        if ~isempty(x_spline)
            mask = (isnan(s_weights));
            if (length(mask) > 2)
                mask = mask | [mask(2:end); 0] | [0; mask(1:end-1)];
            end
            if any(mask)
                s_weights = interp1(x_spline(~mask),s_weights(~mask),x_spline);
            end
            if (size(x_ext,1)==1)
                x_ext = x_ext';
            end
            y_splined_ext = zeros(size(x_ext,1),1);
            for s = 1:length(x_spline)
                tau = round((x_ext-repmat(x_spline(s),length(x_ext),1))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
                y_splined_ext = y_splined_ext + s_weights(s)*cubicSpline(tau);
            end
        else
            y_splined_ext = nan(numel(x_ext),1);
        end
    else
        y_splined_ext = y_splined;
    end
    tmp = nan(numel(inan),1);
    tmp(~inan) = y_splined;
    y_splined = tmp;
end

% No Regularization + variances
function [y_splined, x_spline, s_weights] = spliner_v51(x, y, y_var, dxs)
    nObs = length(x);

    % size of the intervall to interpolate
    x_span = x(nObs) - x(1);
    x0 = x(1);
    x = x - x0;
    
    % compute the number of splines needed for the interpolation
    n_splines = ceil(x_span/dxs) + 3;

    % compute spline centers
    x_spline = zeros(n_splines,1);
    s_weights = [];
    s_center = x(1) - (((n_splines-3)*dxs-x_span)/2) - dxs;
    for i = 1:n_splines
        x_spline(i) = s_center+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(n_splines,1);
    N = sparse(n_splines,n_splines);
    TN = zeros(n_splines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    %tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    n_skip = 0;              % number of spline to "skip" because ain't intersecting an observation
    used_obs = 0;            % number of observation used in building the N matrix
    y_splined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i)-x_spline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            used_obs = used_obs+1;
            n_skip = 0;
            i = i+1;
        else
            cur_local_spline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            n_skip = n_skip +1;
            if (n_skip < 4)
                if (n_skip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    iQ = sparse(diag(1./y_var(first_obs+skips(cur_spline):i-1)));
                    N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) = sparse(N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) + A2'*iQ*A2);

                    % Computing TN
                    TN(cur_local_spline:cur_local_spline+3) = TN(cur_local_spline:cur_local_spline+3) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % If I skip more than 3 times the spline solutions are independent,
                % I can start solving my filtering for the first i points

                s_par = [];
                if (used_obs < size(N,2))
                    fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
                    R = sparse(eye(cur_local_spline)*1e-9);
                    s_par = (N(1:cur_local_spline,1:cur_local_spline)+R)\TN(1:cur_local_spline);
                else
                    s_par = (N(1:cur_local_spline,1:cur_local_spline))\TN(1:cur_local_spline);
                end
                s_weights = [s_weights; s_par];

                for s = first_spline:cur_spline-3
                    y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) = y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * s_par(s-first_spline+1:s+3-first_spline+1);
                end

                % find the next spline whose domain intersect the next observation
                tau = (x(i)-x_spline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-x_spline(cur_spline))/dxs;
                    if (tau > 2),
                        s_weights(cur_spline) = nan;
                    end
                end
                first_spline = cur_spline;
                first_obs = i;
                skips(cur_spline) = 0;

                A = zeros(nObs-i, 4);
                N = sparse(n_splines-cur_spline+1, n_splines-cur_spline+1);
                TN = zeros(n_splines-cur_spline+1,1);
                used_obs = 0;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (n_skip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        iQ = sparse(diag(1./y_var(first_obs+skips(cur_spline):i-1)));
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = sparse(N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*iQ*A2);

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    s_par = [];
    if (used_obs < size(N,2))
        fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = speye(size(N,2), size(N,2))*1e-9;
        s_par = (N+R)\TN;
    else
        s_par = N\TN;
    end
    s_weights = [s_weights; s_par];

    mask = false(size(y_splined));
    for s = first_spline:cur_spline
        if ((skips(s)<skips(s+1))) && (s < 2 || (skips(s) > skips(s-1)))
            y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) = y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * s_par(s-first_spline+1:s+3-first_spline+1);
            mask(skips(s)+first_obs:skips(s+1)+first_obs-1) = true;
        end
    end
    y_splined(~mask) = nan;
    
    x_spline = x_spline + x0;
end

% Regularization + variances
function [y_splined, x_spline, s_weights] = spliner_v51R(x, y, y_var, dxs, reg_factor)
    nObs = length(x);

    y_var(y_var < reg_factor/2) = reg_factor/2;
    % size of the intervall to interpolate
    x_span = x(nObs) - x(1);
    x0 = x(1);
    x = x - x0;
    
    % compute the number of splines needed for the interpolation
    n_splines = ceil(x_span/dxs) + 3;

    % compute spline centers
    x_spline = zeros(n_splines,1);
    s_weights = [];
    s_center = x(1) - (((n_splines-3)*dxs-x_span)/2) - dxs;
    for i = 1:n_splines
        x_spline(i) = s_center+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(n_splines,1);
    N = sparse(n_splines,n_splines);
    TN = zeros(n_splines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    n_skip = 0;              % number of spline to "skip" because ain't intersecting an observation
    used_obs = 0;            % number of observation used in building the N matrix
    y_splined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i)-x_spline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            used_obs = used_obs+1;
            n_skip = 0;
            i = i+1;
        else
            cur_local_spline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            n_skip = n_skip +1;
            if (n_skip < 4)
                if (n_skip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    iQ = sparse(double(diag(1./y_var(first_obs + skips(cur_spline):i-1))));
                    N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) =   sparse(N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) + A2'*iQ*A2);

                    % Computing TN
                    TN(cur_local_spline:cur_local_spline+3) = TN(cur_local_spline:cur_local_spline+3) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % find the next spline whose domain intersect the next observation
                tau = (x(i)-x_spline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-x_spline(cur_spline))/dxs;
                end
                skips(cur_spline) = skips(cur_spline-1);
                used_obs = -1e10;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (n_skip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        iQ = sparse(diag(double(1./y_var(first_obs+skips(cur_spline):i-1))));
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = sparse(N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*iQ*A2);

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    s_par = [];
    if (size(N,2)>2)
        if reg_factor == 0
            fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
            reg_factor = 1e-9;
        end
        R = sparse(eye(size(N,2))-diag(ones(size(N,2)-1,1),1)-diag(ones(size(N,2)-1,1),-1) + diag([0; ones(size(N,2)-2,1); 0]));
        R = R*reg_factor;
        s_par = (N+R)\TN;
    else
        if (used_obs < size(N,2))
            if reg_factor == 0
                fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
                reg_factor = 1e-9;
            end
            R = sparse(eye(size(N,2))*reg_factor);
            s_par = (N+R)\TN;
        else
            s_par = N\TN;
        end
    end
    s_weights = [s_weights; s_par];

    mask = false(size(y_splined));
    for s = first_spline:cur_spline
        if ((skips(s)<skips(s+1))) && (s < 2 || (skips(s) > skips(s-1)))
            y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) = y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * s_par(s-first_spline+1:s+3-first_spline+1);
            mask(skips(s)+first_obs:skips(s+1)+first_obs-1) = true;
        end
    end
    y_splined(~mask) = nan;
    x_spline = x_spline + x0;
end

% No Regularization - no variances
function [y_splined, x_spline, s_weights] = spliner_v5(x, y, dxs)
    nObs = length(x);

    % size of the intervall to interpolate
    x_span = x(nObs) - x(1);
    x0 = x(1);
    x = x - x0;
    
    % compute the number of splines needed for the interpolation
    n_splines = ceil((x_span+eps(x_span))/dxs) + 3;

    % compute spline centers
    x_spline = zeros(n_splines,1);
    s_center = x(1) - (((n_splines-3)*dxs-x_span)/2) - dxs;
    for i = 1:n_splines
        x_spline(i) = s_center+(i-1)*dxs;
    end
    
    % init A matrix
    A = zeros(nObs, 4);
    
    tau   = round(rem(x',dxs)/dxs*1e13)/1e13;  % 1e13 rounding necessary to avoid numerical problems
    idx   = floor((x')/dxs)+1; 
    A     = cubicSpline4Col(tau);
    A_idx = [idx(:) idx(:)+1 idx(:)+2 idx(:)+3];
    n_par = max(A_idx(:,4));
    n_obs = numel(x);
    rows  = repmat((1:n_obs)',1,4);
    
    A = sparse(rows, A_idx, A, n_obs, n_par);
    
    idx_null = sum(A~=0) == 0;
    A(:,idx_null) = [];
    
    N = A'*A;
    B = A'*y;
    
    x = N\B;
    
    s_weights = nan(numel(idx_null),1);
    s_weights(~idx_null) = x;
    y_splined = A*x;

    x_spline = x_spline + x0;
end

% Regularization - no variances
function [y_splined, x_spline, s_weights] = spliner_v5R(x, y, dxs, reg_factor)
    nObs = length(x);

    % size of the intervall to interpolate
    x_span = x(nObs) - x(1);
    x0 = x(1);
    x = x - x0;
    
    % compute the number of splines needed for the interpolation
    n_splines = ceil((x_span+eps(x_span))/dxs) + 3;

    % compute spline centers
    x_spline = zeros(n_splines,1);
    s_weights = [];
    s_center = x(1) - (((n_splines-3)*dxs-x_span)/2) - dxs;
    for i = 1:n_splines
        x_spline(i) = s_center+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(n_splines,1);
    N = sparse(n_splines,n_splines);
    TN = zeros(n_splines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    n_skip = 0;              % number of spline to "skip" because ain't intersecting an observation
    used_obs = 0;            % number of observation used in building the N matrix
    y_splined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i) - x_spline(cur_spline))/dxs * 1e13) / 1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            used_obs = used_obs+1;
            n_skip = 0;
            i = i+1;
        else
            cur_local_spline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            n_skip = n_skip +1;
            if (n_skip < 4)
                if (n_skip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) =   sparse(N(cur_local_spline:cur_local_spline+3,cur_local_spline:cur_local_spline+3) + A2'*A2);

                    % Computing TN
                    TN(cur_local_spline:cur_local_spline+3) = TN(cur_local_spline:cur_local_spline+3) + A2' * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % find the next spline whose domain intersect the next observation
                tau = (x(i)-x_spline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-x_spline(cur_spline))/dxs;
                end
                skips(cur_spline) = skips(cur_spline-1);
                used_obs = -1e10;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (n_skip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*A2;

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    s_par = [];
    if (size(N,2)>2)
        %fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = (speye(size(N,2), size(N,2)) - spdiags(ones(size(N,2), 1), 1, size(N,2), size(N,2)) - spdiags(ones(size(N,2), 1), -1, size(N,2), size(N,2)) + spdiags([0; ones(size(N,2) - 2, 1); 0], 0, size(N,2), size(N,2))) * reg_factor;
        s_par = (N+R)\TN;
    else
        if (used_obs < size(N,2))
            %fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
            R = speye(size(N,2)) * reg_factor;
            s_par = (N+R)\TN;
        else
            s_par = N\TN;
        end
    end
    s_weights = [s_weights; s_par];

    % for s = first_spline:cur_spline
    %     if (skips(s) == 0 && s > 1)
    %         skips(s) = skips(s-1);
    %     end
    %     if ((skips(s)<skips(s+1)))
    %         y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) = y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * s_par(s-first_spline+1:s+3-first_spline+1);
    %     end
    % end

    mask = false(size(y_splined));
    for s = first_spline:cur_spline
        if ((skips(s)<skips(s+1))) && (s < 2 || (skips(s) > skips(s-1)))
            y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) = y_splined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * s_par(s-first_spline+1:s+3-first_spline+1);
            mask(skips(s)+first_obs:skips(s+1)+first_obs-1) = true;
        end
    end
    y_splined(~mask) = nan;
    
    x_spline = x_spline + x0;
end

% SYNTAX:
%   [y] = cubicSpline(t)
%
% EXAMPLE:
%   [y] = cubicSpline(1)
%
% INPUT:
%   t = normalized value from [-2 to 2] of the distance between the
%       observation and the center of the cubic spline
%
% OUTPUT:
%   y = value of the spline in the given point
%
% DESCRIPTION:
%   Get the value of the cubic spline with a base of 4 at the given t normalized
%   value
%
% USEFUL VALUES:
%   2/3 = cubicSpline(0);
%   1/6 = cubicSpline(-1);
%   1/6 = cubicSpline(1);
%
% by Andrea Gatti
%
function [y] = cubicSpline(t)

% FAST VECTORIAL IMPLEMENTATION
y = zeros(size(t));
pos = (t > -2) + (t > -1) - (t > 1) - (t > 2);

% pos = 0    => spline = 0;
% pos = 1    => spline = (-abs(t)+2)^3/6;
% pos = 2    => spline = ((-abs(t)+2)^3 - 4*(-abs(t)+1)^3)/6;

p1 = pos==1;
p2 = pos==2;
t = abs(t);
y(p1) = (2-t(p1)).^3/6;
y(p2) = ((2-t(p2)).^3 - 4*(1-t(p2)).^3)/6;

%  PLAIN IMPLEMENTATION
%
% 	if ((t < -2) || (t > 2))
% 		y = 0;
%     else
%         if (t < 0)
%             if (t <= -1)
%                 tmp2=(t+2);
%                 y = (tmp2*tmp2*tmp2 / 6);
%             else
%                 tmp1=(t+1);
%                 tmp2=(t+2);
%                 y = ((tmp2*tmp2*tmp2 - 4*tmp1*tmp1*tmp1) / 6);
%             end
%         else
%             if (t < 1)
%                 tmp1=(1-t);
%                 tmp2=(2-t);
%                 y = ((tmp2*tmp2*tmp2 - 4*tmp1*tmp1*tmp1) / 6);
%             else
%                 tmp2=(2-t);
%                 y = (tmp2*tmp2*tmp2 / 6);
%             end
%         end
%     end
end

function [val] = cubicSpline4Col(t)
    % Compute matrix entry for cubic spline
    %
    % INPUT
    %   t -> 0 : 1
    %   order -> 1,3
    %
    % SYNTAX:
    %  Core_Utils.cubicSplic(t)
    val = zeros(numel(t),4);
    val(:,1) = (1 - t).^3/6;
    val(:,2) = ((2-t).^3 - 4*(1-t).^3)/6;
    val(:,3) = ((1+t).^3 - 4*(t).^3)/6;
    val(:,4) = (t).^3/6;
end
