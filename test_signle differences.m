[ph ,~, idx_ph_l]=rec(2).getPhases;
idx_ph = find(idx_ph_l);
s_ph = rec(2).getSyntPhObs;
d_ph = ph - s_ph;

sing_diff = (d_ph(:,1) -d_ph(:,3))/trg.wl(idx_ph(1));
sing_diff2 = (d_ph(:,1) -d_ph(:,2))/trg.wl(idx_ph(1));
figure; plot(sing_diff)
series = sing_diff(1:1500);
series2 = sing_diff2(1:1500);
jmp_idx = 1084;
jmp_median = estimateJump(series2,1084,1)
