% setLastLinesWidth get the last lines handler and change them to a <width> value
%
% SINTAX:
%   setLastLinesWidth(<h>,width)
%   setLastLinesWidth(width)
%
% EXAMPLE:
%   setLastLinesWidth(gcf,2);
%   setLastLinesWidth(2);
%
% INPUT:
%   h       = handler to the figure to modify           <optional argument>
%   width   = requested width of the last lines
%
% DEFAULT VALUES:
%   h       = gcf
%
function setLastLinesWidth(h,width)
if nargin < 2
    width = h;
    h = gcf;
end
hline = findobj(h, 'type', 'line'); 
set(hline(1), 'LineWidth', width);
