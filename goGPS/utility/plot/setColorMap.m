function setColorMap(cmap_name, c_limits, perc, thr_low, thr_high)
if nargin < 3
    perc = 0.8;
end
if nargin < 4
    thr_low = [];
    %thrL = [((limitsC(1)-mean(limitsC))*0.05)+mean(limitsC) ((limitsC(2)-mean(limitsC))*0.05)+mean(limitsC)];
end
if nargin < 5
    thr_high = [];
end

%%
set(gcf,'renderer','OpenGL');
cDepth = 2^11;

minC = min(c_limits);
maxC = max(c_limits);

cmap = colormap(Cmap.get(cmap_name, cDepth));
cax = caxis();
rangeC = cax(2)-cax(1);
rangeCnew = maxC-minC;

step = rangeC/cDepth;
caxVal = [cax(1)+step/2:step:cax(2)-step/2]';

[val posMax] = min(abs(caxVal-maxC));
[val posMin] = min(abs(caxVal-minC));


% Set central colormap to perc size
posDiff = posMax-posMin+1;
cmapCenter = Cmap.get(cmap_name, round(posDiff./perc));
cmapCenter = cmapCenter(round((length(cmapCenter)-posDiff)/2):round((length(cmapCenter)-posDiff)/2)+posDiff-1,:);

% Set lateral colormaps to the remaining perc size
cmapLow = Cmap.get(cmap_name, round((posMin)/((1-perc)/2)));
cmapLow = cmapLow(1:posMin-1,:);
cmapHigh = Cmap.get(cmap_name, round((length(caxVal)-posMax)/((1-perc)/2)));
cmapHigh = cmapHigh(end-(length(caxVal)-posMax-1):end,:);
cmap = [cmapLow; cmapCenter; cmapHigh];

% Central threshold
if length(thr_low) == 2
    [val posMax] = min(abs(caxVal-thr_low(2)));
    [val posMin] = min(abs(caxVal-thr_low(1)));
    f = 4;
    g = gray(f*(posMax-posMin+1)); if mod(length(g)/f,2)==1, w=[1 1 1]; else, w = []; end,   g = [g(ceil(end*(f*2-1)/(f*2))+1:end,:); w; flipud(g(ceil(end*(f*2-1)/(f*2))+1:end,:))];
    cmap(posMin:posMax,:)=g; % set to white
end

% High threshold
if length(thr_high) == 2
    [val posMax] = min(abs(caxVal-thr_high(2))); posMax = posMax + 1;
    [val posMin] = min(abs(caxVal-thr_high(1))); posMin = posMin + 1;
    cmap(1:posMin,:)=cmap(1:posMin,:)*0+1; % set to white
    cmap(posMax:end,:)=cmap(posMax:end,:)*0+1; % set to white
end

colormap(cmap);
end
