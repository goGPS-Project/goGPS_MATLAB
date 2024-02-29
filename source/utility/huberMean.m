function mu = huberMean(data, hub_thrs)
data = data(~isnan(data));
if nargin< 2
    hub_thrs = 2.5;
end
w = ones(size(data));
cycle_max = 25;
i = 1;
mu = -inf;
mu_old = inf;
std = 0.1;
while abs(mu - mu_old) > 1e-4*std & i < cycle_max
    mu_old  =  mu;
    mu = sum(data.*w)/sum(w);
    res = data - mu;
    if i == 1
        std = 1.5*median(abs(res));
    else
        std = sum(abs(res).*w)/sum(w);
    end
    w = ones(size(data));
    idx_rw = abs(res)>hub_thrs*std;
    w(idx_rw) = (hub_thrs*std)./abs(res(idx_rw));

    i = i + 1;
end
end