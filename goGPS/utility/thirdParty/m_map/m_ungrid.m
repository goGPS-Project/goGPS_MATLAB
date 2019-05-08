function m_ungrid(goptn)
% M_UNGRID Removes a grid;
%          M_UNGRID deletes a map grid, but leaves any plotted
%          data.
%
%          M_UNGRID XXX
%          or
%          M_UNGRID('XXX')
%
%          can be used to remove other parameters plotted using an
%          M_XXX command (e.g. M_UNGRID COAST).

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/Apr/97 
%
% 14/11/98 - Added possible option to remove other tagged items.
% Nov/2017 - removed a loop around deletion.

%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

if nargin==0
    mstr='m_grid_';
else
    mstr=['m_' lower(goptn)];
    if strncmpi('utm',goptn,3) && strncmp('m_utm_grid',get(gca,'tag'),5)
         set(gca,'visible','off');
         return
    end
end


hh=get(gca,'children');

things=get(hh,'tag');
if length(hh)==1, things={things}; end

htags_del=strncmp(mstr,things,length(mstr));
if ~isempty(htags_del)
    delete(hh(htags_del));
end

if strncmp('m_grid',mstr,6)
  set(gca,'visible','on');
end


