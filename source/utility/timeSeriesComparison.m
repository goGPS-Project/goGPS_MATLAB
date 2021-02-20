function [diff_data] = timeSeriesComparison(t1,data1,t2,data2,mode)
% compare two time series the first one with lower spacing between epochs the second one with finer spacing
if strcmp(mode,'aggregate')
    edges = zeros(1,length(t2)+2+100);
    j = 1;
    rate = round(nan_mean(diff(t1))*86400)/86400;
    empty_idx = [];
    for i = 1:(length(t1) -1)
        edges(j) = t1(i) - (rate/2);
        j = j + 1;
        if ((t1(i+1) - t1(i)) > 1.5*rate)
            edges(j) = t1(i) + (rate/2);
            j = j + 1;
            empty_idx = [empty_idx; j];
        end
    end
    edges(j) = t1(end) + (rate/2);
    edges(j+1:end) = [];
    
    %edges = (data1.time.first.getMatlabTime-eps()) : 1:(data1.time.last.getMatlabTime+1);
    Y = discretize(t2,edges,'IncludedEdge','right');
    idx = false(size(Y));
    for ei = empty_idx'
        idx = idx | Y == ei;
    end
    idx = idx | isnan(Y);
    Y(idx) = [];
    ztds = zero2nan(data2(~idx));
    avg_data = accumarray(Y,ztds,[],@nan_mean);
    avg_data_time = edges(1:end-1) + (edges(2:end) - edges(1:end-1))/2;
    if length(avg_data_time) > max(Y)
    avg_data_time(max(Y)+1:end) = [];
    end
    avg_data_time = avg_data_time(avg_data~=0);
    avg_data = avg_data(avg_data~=0);
    [LIA,LocB] = ismembertol(avg_data_time, t1,1e-8);
    diff_data = nan(size(data1));
    diff_data(LocB(LocB~=0)) = zero2nan(avg_data(LIA)) - zero2nan(data1(LocB(LocB~=0)));
    
elseif strcmp(mode,'interpolate')
    diff_data = data2 - interp1(t1,data1,t2,'linear');
elseif strcmp(mode,'spline')
     edges = zeros(1,length(t2)+2+100);
    j = 1;
    rate = round(nan_mean(diff(t1))*86400)/86400;
    empty_idx = [];
    for i = 1:(length(t1) -1)
        edges(j) = t1(i) - (rate/2);
        j = j + 1;
        if ((t1(i+1) - t1(i)) > 1.5*rate)
            edges(j) = t1(i) + (rate/2);
            j = j + 1;
            empty_idx = [empty_idx; j];
        end
    end
    edges(j) = t1(end) + (rate/2);
    edges(j+1:end) = [];
    
    %edges = (data1.time.first.getMatlabTime-eps()) : 1:(data1.time.last.getMatlabTime+1);
    Y = discretize(t2,edges,'IncludedEdge','right');
    idx = false(size(Y));
    for ei = empty_idx'
        idx = idx | Y == ei;
    end
    idx = idx | isnan(Y);
    Y(idx) = [];
    avg_data_time = edges(1:end-1) + (edges(2:end) - edges(1:end-1))/2;
    if length(avg_data_time) > max(Y)
    avg_data_time(max(Y)+1:end) = [];
    end
    avg_data_time = avg_data_time(unique(Y));
    [LIA,LocB] = ismembertol(avg_data_time, t1,1e-7);
    
    [~,~,~, splined2] = splinerMat(t2, data2,rate,1e-5, avg_data_time);
    diff_data = nan(size(data1));
    diff_data(LocB(LocB~=0)) = zero2nan(splined2(LIA)) - zero2nan(data1(LocB(LocB~=0)));
    
end