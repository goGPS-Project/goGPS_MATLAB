% Special wrapper to regular plot
% Works on regularly sampled data
% When there is a gap of data, it insert a nan value
% to avoid linear interpolation of the data
%
% INPUT
%   t           column array of epochs
%   data        columns of data (could be a matrix)
%   varagin     add other useful parameters of the plot
%
% SYNTAX
%   Core_Utils.plotSep(t, data, varagin);
%
% SEE ALSO
%   plot

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
        
function lh = plotSep(t, data, varargin)
    if not(isempty(data))
        if isa(t, 'handle')
            ax = t;
            t = data;
            data = varargin{1};
            varargin = varargin{2:end};
        else
            ax = gca;
        end
        try
            if numel(t) ~= numel(data) && numel(t) ~= size(data, 1)
                % this means that there is only one data parameter and no time
                varargin = [{data} varargin];
                data = t;
                if size(data, 1) == 1
                    % I want the data to be columnwise
                    data = data';
                end
                t = 1 : size(data, 1);
            else
                if size(data, 1) == 1
                    % I want the data to be columnwise
                    data = data';
                end
            end
        catch
            % probably data is undefined
            data = t;
            if size(data, 1) == 1
                % I want the data to be columnwise
                data = data';
            end
            t = (1 : size(data, 1))';
        end
        for c = 1 : size(data, 2)
            if numel(data) == numel(t)
                [t_col, data_col] = insertNan4Plots(t(:,c), data(:,c));
            else
                [t_col, data_col] = insertNan4Plots(t, data(:,c));
            end
            if isempty(varargin)
                lh = plot(ax, t_col, data_col);
            else
                lh = plot(ax, t_col, data_col, varargin{:});
            end
            hold on;
        end
    end
end

function [t, data_set] = insertNan4Plots(t, data_set)
    % Insert a Nan in a regularly sampled dataset to make
    % plots interrupt continuous lines
    %
    % INPUT
    %   t      epoch of the data [matrix of column arrays]
    %   data   epoch of the data [matrix of column arrays]
    %
    % SYNTAX
    %   [t, data] = Core_Utils.insertNan4Plots(t, data)
    
    t = t(:);
    if size(t, 1) ~= size(data_set, 1)
        % data should be a column array
        data_set = data_set';
    end
    n_set = size(data_set, 2);
    dt = diff(t);
    rate = median(dt);
    id_in = find(dt > 1.5 * rate);
    for x = numel(id_in) : -1 : 1
        t = [t(1 : id_in(x)); (t(id_in(x)) + 1.5 * rate); t((id_in(x)+1) : end)];
        data_set = [data_set(1 : id_in(x), :); nan(1, n_set); data_set((id_in(x)+1) : end, :)];
    end
end
