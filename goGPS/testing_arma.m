n = 500; A = randn(n + 1, n); Qahat = A'*A;  ahat = randn(n,1);
fprintf('------------------------------------------------------------\n');
tic; [Qzhat1, Z1, L1, D1, zhat1, iZt1] = decorrelMex(Qahat, ahat); toc
profile off;
profile on;
tic; [Qzhat2, Z2, L2, D2, zhat2, iZt2] = decorrel_v2(Qahat, ahat); toc;
profile off;
profile viewer;

tic; [Qzhat, Z, L, D, zhat, iZt] = decorrel(Qahat, ahat); toc;
d = (Qzhat1 - Qzhat2)./(Qzhat1 + eps);
d = log10(abs(d));
figure(100); imagesc(d); colorbar; colormap(flipud(Cmap.get('hot')));

