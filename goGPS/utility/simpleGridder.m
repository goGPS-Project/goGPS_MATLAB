% Standalone methods:
%   [gData, wGrid] = simpleGridder.go(phi, lambda , data, gridStep);
%   [gData, wGrid, row, col] = simpleGridder.goCollector(phi, lambda , data, gridStep)
%   [gData, wGrid] = simpleGridder.goXYZ(X, Y, Z, data, gridStep);
classdef simpleGridder < handle
    
    % =========================================================================
    %    PROPERTIES
    % =========================================================================
    
    properties (GetAccess = 'public', SetAccess = 'public')
        data = [];
        flag = [];
        theta = [];
        lambda = [];
        r = [];
    end
    
    % =========================================================================
    %    PUBLIC METHODS
    % =========================================================================
    
    methods
        % -------------------------------------------------------------------------
        % Creator
        % -------------------------------------------------------------------------
        function obj = simpleGridder(data, flag, phi, lambda, r)
            obj.data = data;
            if (nargin == 5)
                obj.flag = flag;
                obj.theta = (pi/2-phi).*180/pi;
                obj.lambda = lambda.*180/pi;
                obj.r = r;
                %if min(obj.lambda) > 0;
                %    obj.lambda = obj.lambda - 180;
                %end
            else
                obj.flag = false(size(obj.data));
                obj.theta = (pi/2-flag).*180/pi;
                obj.lambda = phi.*180/pi;
                obj.r = lambda;
                %if min(obj.lambda) > 0;
                %    obj.lambda = obj.lambda - 180;
                %end
            end
        end
        
        % -------------------------------------------------------------------------
        % Distructor
        % -------------------------------------------------------------------------
        function delete(obj)
        end
    end
    
    methods (Access = 'public')
        function [gData, wGrid, phiGrid, lambdaGrid, row, col] = grid(obj, intL, intT, recursion, showFig)
            if nargin < 3
                % Step of the grid
                intL = 0.5;
                intT = 0.5;
            end
            if nargin < 4
                recursion = false;
            end
            if nargin < 5
                showFig = true;
            end
            
            if showFig
                figure(1); clf;
            end
            rSat = mean(obj.r);
            X = []; Y = []; gData = []; pos=1;
            
            step = 1;
            allIntL = sort(360./divisor(360/intL));
            allIntT = sort(180./divisor(180/intT));
            while (sum(pos) > 0)
                intL = allIntL(min(step,length(allIntL)));
                intT = allIntT(min(step,length(allIntT)));
                
                % coordinates of the grid knots
                lambdaGrid = [(-180+(intL/2)):intL:(180-intL/2)]';
                thetaGrid = [(intT/2):intT:180-(intT/2)]';
                %if (isempty(X))
                %    [X,Y] = meshgrid(lambdaGrid', thetaGrid);
                %else
                %    [x,y] = meshgrid(lambdaGrid', thetaGrid);
                %end
                
                % Find in witch cell the observation falls
                col = max(1,min(floor((obj.lambda+180)/intL)+1,length(lambdaGrid)));
                row = max(1,min(floor(obj.theta/intT)+1,length(thetaGrid)));
                
                g = zeros(length(thetaGrid),length(lambdaGrid));
                wGrid = (zeros(length(thetaGrid),length(lambdaGrid)));
                
                %h = waitbar(0,['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid']);
                fprintf(['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid\n']);
                fprintf('Obs: %8d\\%8d',0,length(obj.data));
                %g(row+(col-1)*size(g,1)) = g(row+(col-1)*size(g,1)) + obj.data;
                %wGrid(row+(col-1)*size(g,1)) = wGrid(row+(col-1)*size(g,1)) + 1; 
                for iExt = 1:100000:length(obj.data)
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%9d\\%9d',iExt,length(obj.data));
                    %waitbar(iExt/length(obj.data),h);
                    for i = iExt:min((iExt+100000), length(obj.data))
                        g(row(i),col(i)) = g(row(i),col(i)) + obj.data(i);
                        wGrid(row(i),col(i)) = wGrid(row(i),col(i)) + 1;
                    end
                end
                g = g./(wGrid+1e-60);
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%9d\\%9d\n',length(obj.data),length(obj.data));
                %close(h);
                if showFig
                    imagesc(g);
                end

                if (isempty(gData))
                    gData = g;
                    lowRes = 0;
                else
                    %lowRes = interp2(x,y,g,X,Y,'linear');
                    %lowRes(isnan(lowRes))=0;
                    fprintf(['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid (%dx%d)\n'],size(g,1),size(g,2));
                    
                    rY = round(size(gData,1)/size(g,1));
                    rX = round(size(gData,2)/size(g,2));
                    lowRes = zeros(size(gData));
                    for dx=1:rX
                        for dy=1:rY
                            lowRes(dy:rY:end,dx:rX:end) = g;
                        end
                    end
                    pos = gData == 0;
                    gData(pos) = lowRes(pos);
                    if showFig
                        fprintf('Showing fig\n');
                        colorLimits = [mean(gData(gData~=0))-3*std(gData(gData~=0)) mean(gData(gData~=0))+3*std(gData(gData~=0))];
                        figure(1); subplot(1,2,1); imagesc(gData); caxis(colorLimits); subplot(1,2,2); imagesc(lowRes); caxis(colorLimits); drawnow;
                    end
                end
                step = step + 1;
                % Force no recurision if asked
                if ~recursion
                    pos = 0;
                end
                
            end
            lambdaGrid = [(-180+(allIntL(1)/2)):allIntL(1):(180-allIntL(1)/2)]';
            thetaGrid = [(allIntT(1)/2):allIntT(1):180-(allIntT(1)/2)]';
            phiGrid = 90-thetaGrid;
        end
        
        function [gData, wGrid, row, col] = gridCollector(obj, intL, intT, recursion, showFig)
            [gData, wGrid, phiGrid, lambdaGrid, row, col] = grid(obj, intL, intT, recursion, showFig);
        end
        
        function [gData, phiGrid, lambdaGrid] = gridD(obj, intL, intT, recursion, showFig)
            if nargin < 3
                % Step of the grid
                intL = 0.5;
                intT = 0.5;
            end
            if nargin < 4
                recursion = true;
            end
            if nargin < 5
                showFig = true;
            end
            
            if showFig
                figure(1); clf;
            end
            rSat = mean(obj.r);
            X = []; Y = []; gData = []; pos=1;
            
            step = 1;
            allIntL = sort(360./divisor(360/intL));
            allIntT = sort(180./divisor(180/intT));
            while (sum(pos) > 0)
                intL = allIntL(min(step,length(allIntL)));
                intT = allIntT(min(step,length(allIntT)));
                
                % coordinates of the grid knots
                lambdaGrid = [(-180+(intL/2)):intL:(180-intL/2)]';
                thetaGrid = [(intT/2):intT:180-(intT/2)]';
                %if (isempty(X))
                %    [X,Y] = meshgrid(lambdaGrid', thetaGrid);
                %else
                %    [x,y] = meshgrid(lambdaGrid', thetaGrid);
                %end
                
                % Find in witch cell the observation falls
                col = max(1,min(floor((obj.lambda+180)/intL)+1,length(lambdaGrid)));
                row = max(1,min(floor(obj.theta/intT)+1,length(thetaGrid)));
                
                % Compute weights
                % Planar distance
                Dx = abs(obj.lambda - (lambdaGrid(col)))/(intL/2);
                Dy = abs(obj.theta - (thetaGrid(row)))/(intT/2);
                d = Dx.^2+Dy.^2;
                
                % Spherical distance at different altitudes
                %d = (obj.r.^2 + rSat^2 - 2*obj.r*rSat .* (sin(obj.theta).*sin(thetaGrid(row)).*cos(obj.lambda-lambdaGrid(col)) + cos(obj.theta).*cos(thetaGrid(row))));
                % Spherical distance at rMean
                %d = (rSat.^2 + rSat^2 - 2*rSat*rSat .* (sin(obj.theta).*sin(thetaGrid(row)).*cos(obj.lambda-lambdaGrid(col)) + cos(obj.theta).*cos(thetaGrid(row))));
                
                W = 1./(d+1e-100);
                %W = (1-(Dx)).*(1-Dy);
                
                g = zeros(length(thetaGrid),length(lambdaGrid));
                wGrid = zeros(length(thetaGrid),length(lambdaGrid));
                
                h = waitbar(0,['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid']);
                fprintf(['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid\n']);
                fprintf('Obs: %8d\\%8d',0,length(obj.data));
                for iExt = 1:100000:length(obj.data)
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8d\\%8d',iExt,length(obj.data));
                    waitbar(iExt/length(obj.data),h);
                    for i = iExt:min((iExt+100000), length(obj.data))
                        g(row(i),col(i)) = g(row(i),col(i)) + obj.data(i)*W(i);
                        wGrid(row(i),col(i)) = wGrid(row(i),col(i)) + W(i);
                    end
                end
                g = g./(wGrid+1e-60);
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8d\\%8d\n',length(obj.data),length(obj.data));
                close(h);
                if (isempty(gData))
                    gData = g;
                    lowRes = 0;
                else
                    %lowRes = interp2(x,y,g,X,Y,'linear');
                    %lowRes(isnan(lowRes))=0;
                    fprintf(['Interpolating on a ' num2str(intT), 'x' num2str(intL), ' grid (%dx%d)\n'],size(g,1),size(g,2));
                    
                    rY = round(size(gData,1)/size(g,1));
                    rX = round(size(gData,2)/size(g,2));
                    lowRes = zeros(size(gData));
                    for dx=1:rX
                        for dy=1:rY
                            lowRes(dy:rY:end,dx:rX:end) = g;
                        end
                    end
                    pos = gData == 0;
                    gData(pos) = lowRes(pos);
                    if showFig
                        colorLimits = [mean(gData(gData~=0))-3*std(gData(gData~=0)) mean(gData(gData~=0))+3*std(gData(gData~=0))];
                        figure(1); subplot(1,2,1); imagesc(gData); caxis(colorLimits); subplot(1,2,2); imagesc(lowRes); caxis(colorLimits); drawnow;
                    end
                end
                step = step + 1;
                % Force no recurision if asked
                if ~recursion
                    pos = 0;
                end
                
            end
            lambdaGrid = [(-180+(allIntL(1)/2)):allIntL(1):(180-allIntL(1)/2)]';
            thetaGrid = [(allIntT(1)/2):allIntT(1):180-(allIntT(1)/2)]';
            phiGrid = 90-thetaGrid;
        end
        
        function correctT(obj)
            % Prepearing data for sinthesis
            load FGOCE_TIM_R5
            gm = 3.986004415e+5;
            rTer = 6.3781363e+3;
            hSat = mean(obj.r)-rTer;
            
            intL = 0.5;
            intT = 0.5;
            
            % coordinates of the grid knots
            lambdaGrid = [(-180+(intL/2)):intL:(180-intL/2)]';
            thetaGrid = [(intT/2):intT:180-(intT/2)]';
            % Find in witch cell the observation falls
            col = max(1,min(floor((obj.lambda+180)/intL)+1,length(lambdaGrid)));
            row = max(1,min(floor(obj.theta/intT)+1,length(thetaGrid)));
            
            fprintf('Sinthesys of Tr\n')
            %[Tr] = sintesiGrid(((90-thetaGrid)./180).*pi, ((lambdaGrid+180)./180).*pi, cnm, snm, 0, 250, 0, 250, gm, rTer, hSat, 1);
            Y = manipulatorInterface(thetaGrid*(pi/180), lambdaGrid*(pi/180), ones(size(thetaGrid))*mean(obj.r)*1e3, 0, 1, 'FGOCE_TIM_R5', 0, 280, 0, 280, 3, 3, '');
            Trrr = Y.Trrr; clear Y;
            fprintf('Removing Trr effect\n')
            
            h = waitbar(0, 'Correcting for altitude');
            for iExt = 1:100000:length(obj.data)
                waitbar(iExt/length(obj.data),h);
                for i = iExt:min((iExt+100000), length(obj.data))
                    obj.data(i) = obj.data(i) - Trrr(row(i),col(i))*(obj.r(i)-(rTer+hSat));
                end
            end
            close(h);
        end
        function correctTrr(obj)
            % Prepearing data for sinthesis
            load FGOCE_TIM_R4
            gm = 3.986004415e+5;
            rTer = 6.3781363e+3;
            hSat = mean(obj.r)-rTer;
            
            intL = 0.5;
            intT = 0.5;
            
            % coordinates of the grid knots
            lambdaGrid = [(-180+(intL/2)):intL:(180-intL/2)]';
            thetaGrid = [(intT/2):intT:180-(intT/2)]';
            % Find in witch cell the observation falls
            col = max(1,min(floor((obj.lambda+180)/intL)+1,length(lambdaGrid)));
            row = max(1,min(floor(obj.theta/intT)+1,length(thetaGrid)));
            
            fprintf('Sinthesys of Tr\n')
            %[Tr] = sintesiGrid(((90-thetaGrid)./180).*pi, ((lambdaGrid+180)./180).*pi, cnm, snm, 0, 250, 0, 250, gm, rTer, hSat, 1);
            Y = manipulatorInterface(thetaGrid*(pi/180), lambdaGrid*(pi/180), ones(size(thetaGrid))*mean(obj.r)*1e3, 0, 1, 'FGOCE_TIM_R4', 0, 250, 0, 250, 1, 0, '');
            Tr = Y.Tr; clear Y;
            fprintf('Removing Tr effect\n')
            
            h = waitbar(0, 'Correcting for altitude');
            for iExt = 1:100000:length(obj.data)
                waitbar(iExt/length(obj.data),h);
                for i = iExt:min((iExt+100000), length(obj.data))
                    obj.data(i) = obj.data(i) - Tr(row(i),col(i))*(obj.r(i)-(rTer+hSat));
                end
            end
            close(h);
        end
        
        function rmEllipsoidT(obj)
            % Prepearing data for sinthesis
            cnm = zeros(11);
            snm = zeros(11);
            
            cnm(1,1) = 1;
            cnm(3,1) = -1.0826298213099999e-03/sqrt(5);
            cnm(5,1) = +2.3709112005300000e-06/sqrt(9);
            cnm(7,1) = -6.0834649888199997e-09/sqrt(13);
            cnm(9,1) = +1.4268108792000000e-11/sqrt(17);
            cnm(11,1) = -1.2143927588200000e-14/sqrt(21);
            save model0_.mat cnm snm
            V = manipulatorInterface(obj.theta*(pi/180), (obj.lambda+180)*(pi/180), obj.r*1e3, 0, 0, 'model0_', 0, 10, 0, 0, 0, 0, '-f');
            delete model0_.mat;
            obj.data = obj.data-V.T;
        end
        function rmEllipsoidTrr(obj)
            % Prepearing data for sinthesis
            cnm = zeros(11);
            snm = zeros(11);
            
            cnm(1,1) = 1;
            cnm(3,1) = -1.0826298213099999e-03/sqrt(5);
            cnm(5,1) = +2.3709112005300000e-06/sqrt(9);
            cnm(7,1) = -6.0834649888199997e-09/sqrt(13);
            cnm(9,1) = +1.4268108792000000e-11/sqrt(17);
            cnm(11,1) = -1.2143927588200000e-14/sqrt(21);
            save model0_.mat cnm snm
            V = manipulatorInterface(obj.theta*(pi/180), (obj.lambda+180)*(pi/180), obj.r*1e3, 0, 0, 'model0_', 0, 10, 0, 0, 2, 3, '-f');
            delete model0_.mat;
            obj.data = obj.data-V.Trr;
        end
    end
    
    methods (Static, Access = 'public')
        function [gData, wGrid, pG, lG] = go(phi, lambda , data, gridStep)
            sg = simpleGridder(data, phi, lambda, phi*0);
            [gData, wGrid, pG, lG] = sg.grid(gridStep(1), gridStep(end), false, false);
            gData(wGrid<1)=nan;
            wGrid(wGrid<1)=nan;
            clear sg;
        end

        function [gData, wGrid, row, col] = goCollector(phi, lambda , data, gridStep)
            sg = simpleGridder(data, phi, lambda, phi*0);
            [gData, wGrid, row, col] = sg.gridCollector(gridStep, gridStep, false, false);
            wGrid(wGrid<1)=nan;
            clear sg;
        end

        function [gData, wGrid, row, col] = goCollectorXYZ(X, Y , Z, data, gridStep)
            [phi, la, r] = cart2polar(X, Y, Z);
            sg = simpleGridder(data, phi, la, phi*0);
            [gData, wGrid, row, col] = sg.gridCollector(gridStep, gridStep, false, false);
            wGrid(wGrid<1)=nan;
            clear sg;
        end
        
        function [gData, wGrid] = goXYZ(X, Y, Z, data, gridStep)
            [phi, la, r] = cart2polar(X, Y, Z);
            sg = simpleGridder(data, phi, la, r);
            [gData, wGrid, pG, lG] = sg.grid(gridStep, gridStep, false, false);
            clear sg;
            gData(wGrid<1)=nan;
            wGrid(wGrid<1)=nan;
        end
    end
    
    
end
