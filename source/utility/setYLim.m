function ax = setYLim(fh, new_lim, axis_id)
% Set the Y limits to all the axis of a figure
% 
% INPUT
%   fh       figurehandle
%   new_lim  new limits for X axis
%   axis_id  number of the axis (1 is the first inserted)
%
% SYNTAX
%   setXLim(fh, new_line, axis_id)
%   setXLim(new_line, axis_id)
%   setXLim(fh, new_line)
%   setXLim(new_line)

%--------------------------------------------------------------------------
%   ______   ______ _______ _    _ _______ (TM)
%   |_____] |_____/ |______  \  /  |_____|
%   |_____] |    \_ |______   \/   |     |
%                                           v0.9999.1
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      ...
%
%  This software is based on GReD's goGPS software
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%-------------------------------------------------------------------------------
%
%    GReD Geomatics Research & Development - Proprietary Use License
%    Version 1.0, - July, 2023
%
%    This software and its associated code, functions, and tools are the property
%    of GReD (Geomatics Research & Development). All rights reserved.
%
%    Redistribution and use of this software in source and binary forms, with or
%    without modification, are strictly prohibited without explicit permission
%    from GReD, except as permitted by law.
%
%    You may not distribute, sublicense, sell, or otherwise transfer any part of
%    this software, code, functions, or tools, nor grant usage rights to any
%    third party without explicit written permission from GReD.
%
%    Unauthorized use, reproduction, or distribution of any part of this software
%    will be considered a breach of the GReD Geomatics Research & Development
%    Proprietary Use License and may result in legal consequences.
%
%    For licensing inquiries or permissions, please contact GReD at:
%    info@g-red.eu
%
%    For the avoidance of doubt, any part of the software not explicitly
%    mentioned under this license remains under the terms of its
%    respective licence
%
%-------------------------------------------------------------------------------
% 1000010 1010010 1000101 1010110 1000001
%--------------------------------------------------------------------------

    if nargin <= 2
        axis_id = [];
    end
    if not(isa(fh, 'matlab.ui.Figure'))
        if nargin == 2
            axis_id = new_lim;
        end
        new_lim = fh;
        fh = gcf;
    end
    set(0, 'CurrentFigure', fh);
    if isempty(axis_id)
        for a = 1 : numel(fh.Children)
            if isa(fh.Children(a), 'matlab.graphics.axis.Axes')
                subplot(fh.Children(a));
                ylim(new_lim);
            end
        end
    else
        try
            ax = fh.Children(max(1, end + 1 - axis_id));
            subplot(ax);
            ylim(new_lim);
        catch
            ax = axes();
        end
    end
end