function wysiwyg
%WYSIWYG -- this function is called with no args and merely
%       changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.
%
%       Particularly useful after calling ORIENT to set a papersize.

%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
%    
% Changes:  Dec/2017 - updated for post-2014b matlab
   
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
set(gcf,'position',pos,'PaperPositionMode','manual');
set(gcf,'units',unis);
  
