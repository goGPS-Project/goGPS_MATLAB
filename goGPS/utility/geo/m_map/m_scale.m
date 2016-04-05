function scale_factor=m_scale(scale_factor),
% M_SCALE Draws the map at a specified scale
%         After M_GRID has been called, the map is sized to fit within
%         the plot region of the current axes. If you want to force the
%         map to appear at a specific scale (e.g. 1:250000), call
%
%         M_SCALE(scale_factor)
%
%         where the map scale is 1:scale_factor. The map will be drawn
%         with its origin at bottom left within the figure limits (which
%         can be set with calls to ORIENT or by setting the 'paperposition'
%         property of the figure).
%
%         M_SCALE('auto') returns to autoscaling.
%
%         SC=M_SCALE returns the current scale (this is useful if you want to
%         know roughly wht scale you currently have)
%
%         see also M_PROJ, M_GRID.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 13/Nov/1998
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

%     1/Feb/99 - added scale calculations for return value.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 4/DEc/11 - isstr to ischar
%

% Need earth radius, in centimeters.
erad=637813700; %cm (from WGS-84)

if nargin==0,

 au=get(gca,'units');
 set(gca,'units','centimeter');
 a_pos=get(gca,'position');
 set(gca,'units',au);

 map_x=diff(get(gca,'xlim'));
 map_y=diff(get(gca,'ylim'));


 if map_x>10, % we are probably using meters (i.e. UTM coords)
   scale_factor=max(map_x/(a_pos(3)/100),map_y/(a_pos(4)/100));
 else
   scale_factor=max(map_x*erad/(a_pos(3)),map_y*erad/(a_pos(4)));
 end;

else
  if ischar(scale_factor),
   gca_pos=get(gca,'userdata');
   if ~isempty(gca_pos),
     set(gca,'position',gca_pos);
   end;
  else

   fu=get(gcf,'units');
   set(gcf,'units','centimeter');
   f_pos=get(gcf,'position');
   set(gcf,'units',fu);

   map_x=diff(get(gca,'xlim'));
   map_y=diff(get(gca,'ylim'));

   if map_x>10, % we are probably using meters (i.e. UTM coords)
     map_x=map_x/scale_factor*100;
     map_y=map_y/scale_factor*100;
   else
     map_x=map_x*erad/scale_factor;
     map_y=map_y*erad/scale_factor;
   end;


   if map_x> f_pos(3) | map_y>f_pos(4),
     disp('Warning - map larger than current window at this scale');
   elseif map_x< f_pos(3)/2 & map_y<f_pos(4)/2,
     disp('Warning - map much smaller than current window at this scale');
   end;

   if isempty(get(gca,'userdata')),
    set(gca,'userdata',get(gca,'position'));
   end;

   gca_u=get(gca,'unit');
   set(gca,'unit','centi','position',[1 1 map_x map_y],'unit',gca_u)

  end;
end;

if nargout==0,
  clear scale_factor
end;



