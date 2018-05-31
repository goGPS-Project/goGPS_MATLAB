function diff = timeSeriesComparison(t1,data1,t2,data2,mode)
% compare two time series the first one with lower spacing between epochs the second one with finer spacing
if strcmp(mode,'aggregate')
    diff = zeros(size(data1));
    t_t =   [-inf ;serialize(t1); inf];
    for i = 1 : length(t_t)-2
        idx_t = t2 > mean(t_t(i:i+1)) & t2 < mean(t_t(i+1:i+2));
        diff(i) = data1(i) - mean(data2(idx_t));
    end
    
elseif strcmp(mode,'interpolate')
    diff = data2 - interp1(t1,data1,t2,'linear');
end
end