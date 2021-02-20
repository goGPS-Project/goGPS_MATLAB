% setColorMapGat(limitsC, perc, thrL, thrH);
function setColorMapGat(limitsC, perc, thrL, thrH)
if nargin < 2
    perc = 0.8;
end
if nargin < 3
    thrL=[];
    %thrL = [((limitsC(1)-mean(limitsC))*0.05)+mean(limitsC) ((limitsC(2)-mean(limitsC))*0.05)+mean(limitsC)];
end
if nargin < 4
    thrH = [];
end

%%
set(gcf,'renderer','OpenGL');
cDepth = 2^11;

minC = min(limitsC);
maxC = max(limitsC);

cmap = colormap(gat(cDepth,1));
cax = caxis();
rangeC = cax(2)-cax(1);
rangeCnew = maxC-minC;

step = rangeC/cDepth;
caxVal = [cax(1)+step/2:step:cax(2)-step/2]';

[val posMax] = min(abs(caxVal-maxC));
[val posMin] = min(abs(caxVal-minC));


% Set central colormap to perc size
posDiff = posMax-posMin+1;
cmapCenter = gat(round(posDiff./perc),1);
cmapCenter = cmapCenter(round((length(cmapCenter)-posDiff)/2):round((length(cmapCenter)-posDiff)/2)+posDiff-1,:);

% Set lateral colormaps to the remaining perc size
cmapLow = gat(round((posMin)/((1-perc)/2)),1);
cmapLow = cmapLow(1:posMin-1,:);
cmapHigh = gat(round((length(caxVal)-posMax)/((1-perc)/2)),1);
cmapHigh = cmapHigh(end-(length(caxVal)-posMax-1):end,:);
cmap = [cmapLow; cmapCenter; cmapHigh];

% Central threshold
if length(thrL) == 2
    [val posMax] = min(abs(caxVal-thrL(2)));
    [val posMin] = min(abs(caxVal-thrL(1)));
    f = 4;
    %g = gray(f*(posMax-posMin+1)); if mod(length(g)/f,2)==1, w=[1 1 1]; else, w = []; end,   g = [g(ceil(end*(f*2-1)/(f*2))+1:end,:); w; flipud(g(ceil(end*(f*2-1)/(f*2))+1:end,:))];
    g = gray(2*f*(posMax-posMin+1)); if mod(length(g)/f,2)==1, w=[1 1 1]; else, w = []; end,   g = [g(ceil(end*(f*2-1)/(f*2))+1:end,:); w;];
    cmap(posMin:posMax,:)=g; % set to white
end

% High threshold
if length(thrH) == 2
    [val posMax] = min(abs(caxVal-thrH(2))); posMax = posMax + 1;
    [val posMin] = min(abs(caxVal-thrH(1))); posMin = posMin + 1;
    cmap(1:posMin,:)=cmap(1:posMin,:)*0+1; % set to white
    cmap(posMax:end,:)=cmap(posMax:end,:)*0+1; % set to white
end

colormap(cmap);
end
