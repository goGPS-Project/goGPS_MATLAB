%----------------------------------------------------------------------------------------------
% REPRESENTATION OF IONOSPHERIC DELAY
%----------------------------------------------------------------------------------------------

corr = zeros(size(elR));
for t = 1 : length(time_GPS)
    index = find(abs(conf_sat(:,t)) == 1)';
    corr(index,t) = iono_error_correction(phi_KAL(t,1), lam_KAL(t,1), azR(index,t), elR(index,t), time_GPS(t), iono, []);
end

coltab = jet(2*nSatTot);
figure; hold on; grid on; title('Ionospheric delay')
xlabel(['epoch [' num2str(interval) ' sec]'])
ylabel('[m]')
for s = 1 : nSatTot
    index = find(abs(conf_sat(s,:)) == 1)';
    h = plot(index,corr(s,index),'r.');
    set(h,'Color',coltab(2*s-1,:));
end