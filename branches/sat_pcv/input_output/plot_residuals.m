function plot_residuals(constellations, RES_PHASE, RES_CODE, outliers_PHASE, outliers_CODE, filerootOUT)


if length(find(isfinite(outliers_PHASE))) > 1
    plot_phase=1;
else
    plot_phase=0;
end

if length(find(isfinite(outliers_CODE))) > 1
    plot_code=1;
else
    plot_code=0;
end


%% GPS
if constellations.GPS.enabled == 1    
    
    %PHASE GRAPHS
    if plot_phase == 1
        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f,'PaperSize');
        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
        for i=1:constellations.GPS.numSat
            
            figure(f);
            subplot(7,5,i);
            hold on;
            title(sprintf('G%02d',constellations.GPS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_PHASE(constellations.GPS.indexes(i),:)*1000,'.b');
            hold on
            plot(find(outliers_PHASE(constellations.GPS.indexes(i),:) == 1),RES_PHASE(constellations.GPS.indexes(i), outliers_PHASE(constellations.GPS.indexes(i),:) == 1)*1000,'.r');
        end
        
        print(f , '-dpdf', [filerootOUT '_GPS_PHASE_residuals']);
        %remove figure
        close(f);
    end
    
    if plot_code == 1
        %CODE GRAPHS        
        f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f1,'PaperSize');
        set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
        
        for i=1:constellations.GPS.numSat
            
            figure(f1);
            subplot(7,5,i);
            hold on;
            title(sprintf('G%02d',constellations.GPS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_CODE(constellations.GPS.indexes(i),:),'.b');
            hold on
            plot(find(outliers_CODE(constellations.GPS.indexes(i),:) == 1),RES_CODE(constellations.GPS.indexes(i), outliers_CODE(constellations.GPS.indexes(i),:) == 1),'.r');
        end
        print(f1, '-dpdf', [filerootOUT '_GPS_CODE_residuals']);
        
        close(f1);
    end
end





%% GLONASS
if constellations.GLONASS.enabled == 1
    
    
    %PHASE GRAPHS
    if plot_phase == 1
        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f,'PaperSize');
        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
        
        for i=1:constellations.GLONASS.numSat
            
            figure(f);
            subplot(6,4,i);
            hold on;
            title(sprintf('R%02d',constellations.GLONASS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_PHASE(constellations.GLONASS.indexes(i),:)*1000,'.b');
            hold on
            plot(find(outliers_PHASE(constellations.GLONASS.indexes(i),:) == 1),RES_PHASE(constellations.GLONASS.indexes(i), outliers_PHASE(constellations.GLONASS.indexes(i),:) == 1)*1000,'.r');
        end
        
        print(f , '-dpdf', [filerootOUT '_GLONASS_PHASE_residuals']);
        %remove figure
        close(f);
    end
    
    
    if plot_code == 1
        %CODE GRAPHS
        f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f1,'PaperSize');
        set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
        
        for i=1:constellations.GLONASS.numSat
            figure(f1);
            subplot(6,4,i);
            hold on;
            title(sprintf('R%02d',constellations.GLONASS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_CODE(constellations.GLONASS.indexes(i),:),'.b');
            hold on
            plot(find(outliers_CODE(constellations.GLONASS.indexes(i),:) == 1),RES_CODE(constellations.GLONASS.indexes(i), outliers_CODE(constellations.GLONASS.indexes(i),:) == 1),'.r');
        end
        print(f1, '-dpdf', [filerootOUT '_GLONASS_CODE_residuals']);
        
        close(f1);
    end
end




%% Galileo
if constellations.Galileo.enabled == 1
    

    %PHASE GRAPHS
    if plot_phase == 1
        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f,'PaperSize');
        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
        
        for i=1:constellations.Galileo.numSat
            figure(f);
            subplot(6,5,i);
            hold on;
            title(sprintf('E%02d',constellations.Galileo.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_PHASE(constellations.Galileo.indexes(i),:)*1000,'.b');
            hold on
            plot(find(outliers_PHASE(constellations.Galileo.indexes(i),:) == 1),RES_PHASE(constellations.Galileo.indexes(i), outliers_PHASE(constellations.Galileo.indexes(i),:) == 1)*1000,'.r');
        end
        
        print(f , '-dpdf', [filerootOUT '_Galileo_PHASE_residuals']);
        %remove figure
        close(f);
    end
    
    
    %CODE GRAPHS
    if plot_code == 1
        f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f1,'PaperSize');
        set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
        
        for i=1:constellations.Galileo.numSat
            %CODE GRAPHS
            figure(f1);
            subplot(6,5,i);
            hold on;
            title(sprintf('E%02d',constellations.Galileo.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_CODE(constellations.Galileo.indexes(i),:),'.b');
            hold on
            plot(find(outliers_CODE(constellations.Galileo.indexes(i),:) == 1),RES_CODE(constellations.Galileo.indexes(i), outliers_CODE(constellations.Galileo.indexes(i),:) == 1),'.r');
        end
        print(f1, '-dpdf', [filerootOUT '_Galileo_CODE_residuals']);
        
        close(f1);
    end
    
end



%% Galileo
if constellations.BeiDou.enabled == 1
    
    
    %PHASE GRAPHS
    if plot_phase == 1
        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f,'PaperSize');
        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
        
        for i=1:constellations.BeiDou.numSat
            figure(f);
            subplot(7,6,i);
            hold on;
            title(sprintf('C%02d',constellations.BeiDou.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_PHASE(constellations.BeiDou.indexes(i),:)*1000,'.b');
            hold on
            plot(find(outliers_PHASE(constellations.BeiDou.indexes(i),:) == 1),RES_PHASE(constellations.BeiDou.indexes(i), outliers_PHASE(constellations.BeiDou.indexes(i),:) == 1)*1000,'.r');
        end
        
        print(f , '-dpdf', [filerootOUT '_BeiDou_PHASE_residuals']);
        %remove figure
        close(f);
    end
    
    %CODE GRAPHS
    if plot_code == 1
        f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f1,'PaperSize');
        set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
        
        for i=1:constellations.BeiDou.numSat
            
            figure(f1);
            subplot(7,6,i);
            hold on;
            title(sprintf('C%02d',constellations.BeiDou.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_CODE(constellations.BeiDou.indexes(i),:),'.b');
            hold on
            plot(find(outliers_CODE(constellations.BeiDou.indexes(i),:) == 1),RES_CODE(constellations.BeiDou.indexes(i), outliers_CODE(constellations.BeiDou.indexes(i),:) == 1),'.r');
        end
        print(f1, '-dpdf', [filerootOUT '_BeiDou_CODE_residuals']);
        close(f1)
    end;

end



%% QZSS
if constellations.QZSS.enabled == 1
    
    
    %PHASE GRAPHS
    if plot_phase == 1
        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f,'PaperSize');
        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
        
        for i=1:constellations.QZSS.numSat
            %PHASE GRAPHS
            figure(f);
            subplot(2,2,i);
            hold on;
            title(sprintf('J%02d',constellations.QZSS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_PHASE(constellations.QZSS.indexes(i),:)*1000,'.b');
            hold on
            plot(find(outliers_PHASE(constellations.QZSS.indexes(i),:) == 1),RES_PHASE(constellations.QZSS.indexes(i), outliers_PHASE(constellations.QZSS.indexes(i),:) == 1)*1000,'.r');
        end
        
        print(f , '-dpdf', [filerootOUT '_QZSS_PHASE_residuals']);
        %remove figure
        close(f);
    end
    
    %CODE GRAPHS
    if plot_code == 1
        f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
        paperSize = get(f1,'PaperSize');
        set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
        
        for i=1:constellations.QZSS.numSat
            figure(f1);
            subplot(2,2,i);
            hold on;
            title(sprintf('C%02d',constellations.QZSS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
            grid on;
            set(gca,'FontName','Verdana');
            set(gca,'FontSize',7);
            ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
            plot(RES_CODE(constellations.QZSS.indexes(i),:),'.b');
            hold on
            plot(find(outliers_CODE(constellations.QZSS.indexes(i),:) == 1),RES_CODE(constellations.QZSS.indexes(i), outliers_CODE(constellations.QZSS.indexes(i),:) == 1),'.r');
        end
        print(f1, '-dpdf', [filerootOUT '_QZSS_CODE_residuals']);
        
        close(f1);
    end

end






