function m_ungrid(goptn);
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
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

if nargin==0,
  mstr='m_grid_';
else
  mstr=['m_' lower(goptn)];
end;

hh=get(gca,'children');

things=get(hh,'tag');
if length(hh)==1, things={things}; end;

for i=1:length(hh),
  if ~isempty(things{i}) & strmatch(mstr,things{i}),
   delete(hh(i));
 end;
end;

if strmatch('m_grid',mstr),
  set(gca,'visible','on');
end;


