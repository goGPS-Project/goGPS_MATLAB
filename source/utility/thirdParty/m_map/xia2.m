
%% Pic2:
figure(1)
subplot(121)
m_proj('ortho','lat',50','long',120');
m_coast('patch','r');
m_grid('linest','-','xticklabels',[],'yticklabels',[]);

patch(.55*[-1 1 1 -1],.25*[-1 -1 1 1]-.55,'w'); 
text(0,-.55,'ÖÐ¹úµØÍ¼','fontsize',25,'color','b',...
'verticalalignment','middle','horizontalalignment','center');



%% Pic1:
subplot(122)

m_proj('lambert','long',[80 180],'lat',[30 80]);
m_coast('patch',[1 .85 .7]);
[CS,CH]=m_elev('contourf',[500:500:4000]);
%   m_elev('pcolor');
m_grid('box','fancy','tickdir','in');
colormap(flipud(copper));
xlabel('Conic Projection of North America with elevations','visible','on');
m_contfbar([0 .3],.9,CS,CH);
