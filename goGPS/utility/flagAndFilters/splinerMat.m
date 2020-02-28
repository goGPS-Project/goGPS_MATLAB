function [ySplined, xSpline, sWeights, ySplined_ext] = splinerMat(x,y,dxs,regFactor,x_ext)
% SYNTAX:
%   [ySplined, xSpline, sWeights ySplined_ext] = splinerMat(x,y,dxs,regFactor,<x_ext>)
%
% EXAMPLE:
%   [ySplined, xSpline, sWeights, ySplined_ext] = splinerMat(x,y,4,0,x_ext);
%   [ySplined, xSpline, sWeights, ySplined_ext] = splinerMat(x,[y yvar],4,0,x_ext);
%
% INPUT:
%   x            [n x 1] observation time
%   y            [n x 1] observation value
%                [n x 2] observation value and variances (in this case the variances of the data are taken into account)
%   dxs          [1 x 1] spline base size
%   regFactor    [1 x 1] regularization factor on the first derivative
%   x_ext        [m x 1] points in which to compute the interpolation
%
% OUTPUT:
%   ySplined     [n x 1] interpolated observation (on the observation epochs)
%   xSpline      [o x 1] center of the (o) splines
%   sWeights     [o x 1] weights of the splines
%   ySplined_ext [m x 1] spline interpolated in x_ext positions
%
% DESCRIPTION:
%   Interpolate with cubic splines a given dataset
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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
        regFactor = 0;
    end
    if isempty(x)
        x = 1:size(y, 1);
    end
    
    inan = isnan(y);
    if (size(y,2) == 2)
        inan = inan(:,1) | inan(:,2);
    end
    
    x(inan) = [];
    y(inan) = [];
    
    [x, id] = sort(x);
    y = y(id);
    if ~isempty(y)
        if ((nargin == 3) || (regFactor == 0))
            if (size(y,2) == 2)
                [ySplined, xSpline, sWeights] = spliner_v51(x,y(:,1),y(:,2),dxs);
            else
                [ySplined, xSpline, sWeights] = spliner_v5(x,y,dxs);
            end
        else
            if (size(y,2) == 2)
                [ySplined, xSpline, sWeights] = spliner_v51R(x,y(:,1),y(:,2),dxs, regFactor);
            else
                [ySplined, xSpline, sWeights] = spliner_v5R(x,y,dxs, regFactor);
            end
        end
    else
        ySplined = y;
        xSpline = [];
        sWeights = [];
    end

    % Interpolation => using spline to predict in different coordinates
    % (not present in the C version)
    if (nargin == 5)
        if ~isempty(xSpline)
            mask = (isnan(sWeights));
            if (length(mask) > 2)
                mask = mask | [mask(2:end); 0] | [0; mask(1:end-1)];
            end
            sWeights = interp1(xSpline(~mask),sWeights(~mask),xSpline);
            if (size(x_ext,1)==1)
                x_ext = x_ext';
            end
            ySplined_ext = zeros(size(x_ext,1),1);
            for s = 1:length(xSpline)
                tau = round((x_ext-repmat(xSpline(s),length(x_ext),1))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
                ySplined_ext = ySplined_ext + sWeights(s)*cubicSpline(tau);
            end
        else
            ySplined_ext = nan(numel(x_ext),1);
        end
    else
        ySplined_ext = ySplined;
    end
    tmp = nan(numel(inan),1);
    tmp(~inan) = ySplined;
    ySplined = tmp;
end

% No Regularization + variances
function [ySplined, xSpline, sWeights] = spliner_v51(x,y,yvar,dxs)
    nObs = length(x);

    % size of the intervall to interpolate
    xspan = x(nObs) - x(1);

    % compute the number of splines needed for the interpolation
    nSplines = ceil(xspan/dxs) + 3;

    % compute spline centers
    xSpline = zeros(nSplines,1);
    sWeights = [];
    sCenter = x(1) - (((nSplines-3)*dxs-xspan)/2) - dxs;
    for i = 1:nSplines
        xSpline(i) = sCenter+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(nSplines,1);
    N = sparse(nSplines,nSplines);
    TN = zeros(nSplines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    %tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    nSkip = 0;              % number of spline to "skip" because ain't intersecting an observation
    usedObs = 0;            % number of observation used in building the N matrix
    ySplined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i)-xSpline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            usedObs = usedObs+1;
            nSkip = 0;
            i = i+1;
        else
            curLocalSpline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            nSkip = nSkip +1;
            if (nSkip < 4)
                if (nSkip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    iQ = sparse(diag(1./yvar(first_obs+skips(cur_spline):i-1)));
                    N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) = sparse(N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) + A2'*iQ*A2);

                    % Computing TN
                    TN(curLocalSpline:curLocalSpline+3) = TN(curLocalSpline:curLocalSpline+3) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % If I skip more than 3 times the spline solutions are independent,
                % I can start solving my filtering for the first i points

                sPar = [];
                if (usedObs < size(N,2))
                    fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
                    R = sparse(eye(curLocalSpline)*1e-9);
                    sPar = (N(1:curLocalSpline,1:curLocalSpline)+R)\TN(1:curLocalSpline);
                else
                    sPar = (N(1:curLocalSpline,1:curLocalSpline))\TN(1:curLocalSpline);
                end
                sWeights = [sWeights; sPar];

                for s = first_spline:cur_spline-3
                    ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
                end

                % find the next spline whose domain intersect the next observation
                tau = (x(i)-xSpline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-xSpline(cur_spline))/dxs;
                    if (tau > 2),
                        sWeights(cur_spline) = nan;
                    end
                end
                first_spline = cur_spline;
                first_obs = i;
                skips(cur_spline) = 0;

                A = zeros(nObs-i, 4);
                N = sparse(nSplines-cur_spline+1, nSplines-cur_spline+1);
                TN = zeros(nSplines-cur_spline+1,1);
                usedObs = 0;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (nSkip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        iQ = sparse(diag(1./yvar(first_obs+skips(cur_spline):i-1)));
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = sparse(N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*iQ*A2);

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    sPar = [];
    if (usedObs < size(N,2))
        fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = speye(size(N,2), size(N,2))*1e-9;
        sPar = (N+R)\TN;
    else
        sPar = N\TN;
    end
    sWeights = [sWeights; sPar];

    for s = first_spline:cur_spline
        if ((skips(s)<skips(s+1)))
            ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
        end
    end
end

% Regularization + variances
function [ySplined, xSpline, sWeights] = spliner_v51R(x, y, yvar, dxs, regFactor)
    xSpline = [];
    nObs = length(x);

    % size of the intervall to interpolate
    xspan = x(nObs) - x(1);

    % compute the number of splines needed for the interpolation
    nSplines = ceil(xspan/dxs) + 3;

    % compute spline centers
    xSpline = zeros(nSplines,1);
    sWeights = [];
    sCenter = x(1) - (((nSplines-3)*dxs-xspan)/2) - dxs;
    for i = 1:nSplines
        xSpline(i) = sCenter+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(nSplines,1);
    N = sparse(nSplines,nSplines);
    TN = zeros(nSplines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    nSkip = 0;              % number of spline to "skip" because ain't intersecting an observation
    usedObs = 0;            % number of observation used in building the N matrix
    ySplined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i)-xSpline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            usedObs = usedObs+1;
            nSkip = 0;
            i = i+1;
        else
            curLocalSpline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            nSkip = nSkip +1;
            if (nSkip < 4)
                if (nSkip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    iQ = sparse(diag(1./yvar(first_obs + skips(cur_spline):i-1)));
                    N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) =   sparse(N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) + A2'*iQ*A2);

                    % Computing TN
                    TN(curLocalSpline:curLocalSpline+3) = TN(curLocalSpline:curLocalSpline+3) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % find the next spline whose domain intersect the next observation
                tau = (x(i)-xSpline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-xSpline(cur_spline))/dxs;
                end
                skips(cur_spline) = skips(cur_spline-1);
                usedObs = -1e10;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (nSkip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        iQ = sparse(diag(1./yvar(first_obs+skips(cur_spline):i-1)));
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = sparse(N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*iQ*A2);

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * iQ * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    sPar = [];
    if (size(N,2)>2)
        fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = sparse(eye(size(N,2))-diag(ones(size(N,2)-1,1),1)-diag(ones(size(N,2)-1,1),-1) + diag([0; ones(size(N,2)-2,1); 0]));
        R = R*regFactor;
        sPar = (N+R)\TN;
    else
        if (usedObs < size(N,2))
            fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
            R = sparse(eye(size(N,2))*regFactor);
            sPar = (N+R)\TN;
        else
            sPar = N\TN;
        end
    end
    sWeights = [sWeights; sPar];

    for s = first_spline:cur_spline
        if ((skips(s)<skips(s+1)))
            ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
        end
    end
end

% No Regularization - no variances
function [ySplined, xSpline, s_weights] = spliner_v5(x,y,dxs)
    xSpline = [];
    nObs = length(x);

    % size of the intervall to interpolate
    xspan = x(nObs) - x(1);

    % compute the number of splines needed for the interpolation
    nSplines = ceil(xspan/dxs) + 3;

    % compute spline centers
    xSpline = zeros(nSplines,1);
    s_weights = [];
    sCenter = x(1) - (((nSplines-3)*dxs-xspan)/2) - dxs;
    for i = 1:nSplines
        xSpline(i) = sCenter+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(nSplines,1);
    N = sparse(nSplines, nSplines);
    TN = zeros(nSplines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    n_skip = 0;              % number of spline to "skip" because ain't intersecting an observation
    usedObs = 0;            % number of observation used in building the N matrix
    ySplined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i)-xSpline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            usedObs = usedObs+1;
            n_skip = 0;
            i = i+1;
        else
            curLocalSpline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            n_skip = n_skip +1;
            if (n_skip < 4)
                if (n_skip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) =   sparse(N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) + A2'*A2);

                    % Computing TN
                    TN(curLocalSpline:curLocalSpline+3) = TN(curLocalSpline:curLocalSpline+3) + A2' * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % If I skip more than 3 times the spline solutions are independent,
                % I can start solving my filtering for the first i points

                if (usedObs < size(N,2))
                    fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
                    R = sparse(eye(curLocalSpline)*1e-9);
                    sPar = (N(1:curLocalSpline,1:curLocalSpline)+R)\TN(1:curLocalSpline);
                else
                    sPar = (N(1:curLocalSpline,1:curLocalSpline))\TN(1:curLocalSpline);
                end
                s_weights = [s_weights; sPar];

                for s = first_spline:cur_spline-3
                    ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
                end

                % find the next spline whose domain intersect the next observation
                tau = round((x(i)-xSpline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = round((x(i)-xSpline(cur_spline))/dxs *1e13)/1e13; % 1e13 rounding necessary to avoid numerical problems
                    if (tau > 2)
                        s_weights(cur_spline) = nan;
                    end
                end
                first_spline = cur_spline;
                first_obs = i;
                skips(cur_spline) = 0;

                A = zeros(nObs-i, 4);
                N = sparse(nSplines-cur_spline+1,nSplines-cur_spline+1);
                TN = zeros(nSplines-cur_spline+1,1);
                usedObs = 0;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (n_skip == 0)
        last_spline = cur_spline+3;
        A2 = A((skips(cur_spline)+1):i-first_obs,:);
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = sparse(N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*A2);

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    if (usedObs < size(N,2))
        fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = sparse(eye(size(N,2))*1e-9);
        sPar = (N+R)\TN;
    else
        sPar = N\TN;
    end
    s_weights = [s_weights; sPar];

    for s = first_spline:cur_spline
        if (skips(s) == 0 && s > 1)
            skips(s) = skips(s-1);
        end
        if ((skips(s)<skips(s+1)))
            ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
        end
    end
end

% Regularization - no variances
function [ySplined, xSpline, sWeights] = spliner_v5R(x,y,dxs, regFactor)
    nObs = length(x);

    % size of the intervall to interpolate
    xspan = x(nObs) - x(1);

    % compute the number of splines needed for the interpolation
    nSplines = ceil(xspan/dxs) + 3;

    % compute spline centers
    xSpline = zeros(nSplines,1);
    sWeights = [];
    sCenter = x(1) - (((nSplines-3)*dxs-xspan)/2) - dxs;
    for i = 1:nSplines
        xSpline(i) = sCenter+(i-1)*dxs;
    end

    % init A matrix
    A = zeros(nObs, 4);
    skips = zeros(nSplines,1);
    N = sparse(nSplines,nSplines);
    TN = zeros(nSplines, 1);

    cur_spline = 1;          % first spline whose domain intersect the observation
    tau = 0;                % normalized distance between the observation and the center of the cur_spline
    i = 1;                  % index of the first observation
    first_obs = i;          % first observation used in the current A matrix
    first_spline = 1;       % first spline used in the current A matrix
    nSkip = 0;              % number of spline to "skip" because ain't intersecting an observation
    usedObs = 0;            % number of observation used in building the N matrix
    ySplined = zeros(length(y),1);           % output
    skips(1) = 0;
    while (i <= nObs)
        % Compute the distance between the current observation and the current spline
        tau = round((x(i) - xSpline(cur_spline))/dxs * 1e13) / 1e13; % 1e13 rounding necessary to avoid numerical problems
        if (tau <= 2)
            % fill the design matrix
            A(i-first_obs+1,:) = cubicSpline([tau (tau-1) (tau-2) (tau-3)]);
            usedObs = usedObs+1;
            nSkip = 0;
            i = i+1;
        else
            curLocalSpline = cur_spline-first_spline+1;
            skips(cur_spline+1) = i-first_obs;
            % This block of the A matrix is completed
            % Computing N
            nSkip = nSkip +1;
            if (nSkip < 4)
                if (nSkip == 1)
                    A2 = A((skips(cur_spline)+1):i-first_obs,:);
                    N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) =   sparse(N(curLocalSpline:curLocalSpline+3,curLocalSpline:curLocalSpline+3) + A2'*A2);

                    % Computing TN
                    TN(curLocalSpline:curLocalSpline+3) = TN(curLocalSpline:curLocalSpline+3) + A2' * y(first_obs+skips(cur_spline):i-1);
                end
                cur_spline = cur_spline +1;
            else
                % find the next spline whose domain intersect the next observation
                tau = (x(i)-xSpline(cur_spline))/dxs;
                while (tau > 2)
                    cur_spline = cur_spline+1;
                    tau = (x(i)-xSpline(cur_spline))/dxs;
                end
                skips(cur_spline) = skips(cur_spline-1);
                usedObs = -1e10;
            end
        end
    end

    skips(cur_spline+1) = i-first_obs;
    if (nSkip == 0)
        A2 = A(skips(cur_spline)+1:i-first_obs,:);
        N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) = N(cur_spline-first_spline+1:cur_spline+3-first_spline+1,cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2'*A2;

        % Computing TN
        TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) = TN(cur_spline-first_spline+1:cur_spline+3-first_spline+1) + A2' * y(first_obs+skips(cur_spline):i-1);
    end

    % find the interpolation for the last subset of observations
    sPar = [];
    if (size(N,2)>2)
        %fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
        R = (speye(size(N,2), size(N,2)) - spdiags(ones(size(N,2), 1), 1, size(N,2), size(N,2)) - spdiags(ones(size(N,2), 1), -1, size(N,2), size(N,2)) + spdiags([0; ones(size(N,2) - 2, 1); 0], 0, size(N,2), size(N,2))) * regFactor;
        sPar = (N+R)\TN;
    else
        if (usedObs < size(N,2))
            %fprintf('WARNING: Regularization is needed observations are less than splines.\n         Adding 1e-9 on the normal matrix diagonal\n');
            R = speye(size(N,2)) * regFactor;
            sPar = (N+R)\TN;
        else
            sPar = N\TN;
        end
    end
    sWeights = [sWeights; sPar];

    for s = first_spline:cur_spline
        if (skips(s) == 0 && s > 1)
            skips(s) = skips(s-1);
        end
        if ((skips(s)<skips(s+1)))
            ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) = ySplined(skips(s)+first_obs:skips(s+1)+first_obs-1) + A((skips(s):skips(s+1)-1)+1,:) * sPar(s-first_spline+1:s+3-first_spline+1);
        end
    end
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
