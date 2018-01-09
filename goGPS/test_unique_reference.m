%test unique reference
clear
close all
idx_i = [0 0 0 0  0 0 0 0 0 1 1 0 0 0   0 0 0 0 0 0 0 0 zeros(1,100) ones(1,100) 0 ones(1,100) zeros(1,100) 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 zeros(1,100) ];
idx = flagShrink(idx_i',4)
idx = idx';

idxe = idx;
filt = exp(-(-2:0.5:2).^2);
filt = filt./sum(filt);
filt_dx= filt;
filt_dx(1:10) = 0;
filt_dx = filt_dx ./ sum(filt_dx);
filt_sx= filt;
filt_sx(11:end) = 0;
filt_sx = filt_sx ./ sum(filt_sx);
idx_lf = conv(idxe,filt_sx,'same');
idx_rg = conv(idxe,filt_dx,'same');
idx_f = conv(idxe,filt,'same');
figure
idx1 = min(circshift(idx_f',4),circshift(idx_f',-4));
idx1(~idx) = 0;
plot(conv(idx1,filt,'same'),'r')
plot(idx1,'b')
hold on
%  plot(circshift(idx_f',11),'b')
%  plot(circshift(idx_f',-11),'g')

plot(idx_i,'k')