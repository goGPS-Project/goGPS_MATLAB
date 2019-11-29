t0 = tic;
core.exec({'EMPTY T1'});
toc(t0)
tic
core.exec({'LOAD T1 @30s'});
toc
%profile off; profile('-timer', 'cpu', '-memory', 'off', '-timestamp', 'off'); 
%profile off; profile('-timer', 'cpu'); profile on
%profile off; profile('-timer', 'real'); profile on
%t0 = tic;
%core.exec({'PREPRO T* -s=G'});
%toc(t0)
%profile off; profile viewer
%
t0 = tic;
core.rec.work.preProcessing()
%%core.exec({'PREPRO T1 -s=G'});
toc(t0)
tic
%%
%rec(3).work.reset2AprioriTropo
%rec(3).work.staticPPP('G', 1:2760);
%rec(3).work.staticPPPNew('G', 1:2760);
%rec(3).work.staticPPPNew('G');
rec(3).work.staticPPP('G');
%core.exec({'PPP T3'});
toc
rec(3).work.showZwd(false);
%setAllLinesStyle('-')
%%
tic
core.exec({'PUSHOUT T1'});
toc
toc(t0)
tic
core.exec({'SHOW ZTD T1'});
toc
%profile off; profile on; tic;  toc; profile off; profile viewer
%%
tic
core.exec({'PPP T1'});
toc
tic
core.exec({'PPP T2 -U'});
toc
%%
profile off; profile('-timer', 'real'); profile on
core.exec({'PPP T1 -U'});
profile off; profile viewer
%%

tmp = zero2nan(this.getDtS(go_id));
for s = 1: numel(go_id)
    tmp0(this.sat.avail_index(:,go_id(s)),s) = zero2nan(this.getDtS(go_id(s)));
end
% PREPRO TEST

res = rec(2).work.sat.res;
sensor = zero2nan(diff(Core_Utils.diffAndPred(movmedian(res, 5), 1, [], 'linear')));
figure; plot(sensor);
% find the worse satellite for each epoch
[~, id_worse] = max(abs(sensor)');
% create a leave one out set
loo = sensor;
loo((1:numel(id_worse)) + (id_worse - 1) .* numel(id_worse)) = nan;
m_loo = mean(zero2nan(res), 2, 'omitnan');
figure; plot(sensor - m_loo)


msensor = median(zero2nan(sensor), 2, 'omitnan');
ssensor = std(zero2nan(sensor)', 'omitnan');
sensor = bsxfun(@minus, sensor, msensor);
figure; plot(msensor);
hold on; plot(ssensor - movmedian(ssensor, 3));

