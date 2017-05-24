%% plot residual wrt elevation/azimuth
% fprintf('                 - Residuals skyplot ...');
f = figure('Name','Phase residuals skyplot','NumberTitle','Off','MenuBar','None','Position',[26 79 967 603],'Visible','off');
title('Phase residuals skyplot','FontName','Verdana','FontSize',14,'FontWeight','Bold','Color',[0 0 1]);
hold on
scatter(azR(:),elR(:),6,(abs(RES_PHASE1_FIXED(:))));
plot(azR(outliers_PHASE1==1),elR(outliers_PHASE1==1),'ok');
set(gca,'xlim',[0 360]);
set(gca,'ylim',[0 90]);
set(gca,'XTick',0:30:360)
grid on
xlabel('Azimuth','FontName','Verdana','FontSize',10,'FontWeight','Bold');
ylabel('Elevation','FontName','Verdana','FontSize',10,'FontWeight','Bold');
caxis([0 max_phase_residual]);
cb = colorbar;
cblabel=get(cb,'YTickLabel');
cblabel(end)=cellstr(['>' num2str(max_phase_residual,'%.2f') 'm']);
set(cb,'YTickLabel', cblabel);
colormap(goGNSS.RES_COLORMAP);
annotation(f,'textbox',[0.80 0.001 0.2381 0.04638],...
    'String',{['ÿgoGPS, ',datestr(now)]},...
    'HorizontalAlignment','center',...
    'FontSize',6,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 0 0]);

%print PNG
print(f, '-dpng', [filerootOUT '_PHASE_residuals_skyplot']);
%print PDF
set(f,'PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4');
paperSize = get(f,'PaperSize');
set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
print(f, '-dpdf', [filerootOUT '_PHASE_residuals_skyplot']);
%save FIG
set(f,'visible','on'); set(f,'MenuBar','figure');
saveas(f, [filerootOUT '_PHASE_residuals_skyplot'], 'fig');
% close(f);


f = figure('Name','Code residuals skyplot','NumberTitle','Off','MenuBar','None','Position',[26 79 967 603],'Visible','off');
title('Code residuals skyplot','FontName','Verdana','FontSize',14,'FontWeight','Bold','Color',[0 0 1]);
hold on
scatter(azR(:),elR(:),6,(abs(RES_CODE1_FIXED(:))));
plot(azR(outliers_CODE1==1),elR(outliers_CODE1==1),'ok');
set(gca,'xlim',[0 360]);
set(gca,'ylim',[0 90]);
set(gca,'XTick',0:30:360)
grid on
xlabel('Azimuth','FontName','Verdana','FontSize',10,'FontWeight','Bold');
ylabel('Elevation','FontName','Verdana','FontSize',10,'FontWeight','Bold');
caxis([0 max_code_residual]);
cb = colorbar;
cblabel=get(cb,'YTickLabel');
cblabel(end)=cellstr(['>' num2str(max_code_residual,'%d') 'm']);
set(cb,'YTickLabel', cblabel);
colormap(goGNSS.RES_COLORMAP);
annotation(f,'textbox',[0.80 0.001 0.2381 0.04638],...
    'String',{['ÿgoGPS, ',datestr(now)]},...
    'HorizontalAlignment','center',...
    'FontSize',6,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 0 0]);
%print PNG
print(f, '-dpng', [filerootOUT '_CODE_residuals_skyplot']);
%print PDF
set(f,'PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4');
paperSize = get(f,'PaperSize');
set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
print(f, '-dpdf', [filerootOUT '_CODE_residuals_skyplot']);
%save FIG
set(f,'visible','on'); set(f,'MenuBar','figure');
saveas(f, [filerootOUT '_CODE_residuals_skyplot'], 'fig');
% close(f);

save([filerootOUT '_residuals'], 'azR' , 'elR', 'RES_CODE1_FIXED', 'RES_PHASE1_FIXED', 'outliers_PHASE1', 'outliers_CODE1');

% fprintf(' done\n');
