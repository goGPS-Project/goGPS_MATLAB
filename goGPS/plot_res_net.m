figure; 
[~,wl,idx_ph] = rec(1).work.getPhases(); idx_ph = find(idx_ph);

o_code = 'L1C'
subplot(6,1,1); plot(rec(1).work.sat.res_ph_by_ph(:,rec(1).work.system(idx_ph)' == 'G' & strLineMatch(rec(1).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec1 G' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15); 
o_code = 'L2L'
subplot(6,1,2); plot(rec(1).work.sat.res_ph_by_ph(:,rec(1).work.system(idx_ph)' == 'G' & strLineMatch(rec(1).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec1 G' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15); 
o_code = 'L1C'
[~,wl,idx_ph] = rec(2).work.getPhases(); idx_ph = find(idx_ph);

subplot(6,1,3); plot(rec(2).work.sat.res_ph_by_ph(:,rec(2).work.system(idx_ph)' == 'E' & strLineMatch(rec(2).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec2 E' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'L7Q'
[~,wl,idx_ph] = rec(3).work.getPhases(); idx_ph = find(idx_ph);

subplot(6,1,4); plot(rec(3).work.sat.res_ph_by_ph(:,rec(3).work.system(idx_ph)' == 'E' & strLineMatch(rec(3).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec3 E' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'L1C'
[~,wl,idx_ph] = rec(4).work.getPhases(); idx_ph = find(idx_ph);
subplot(6,1,5); plot(rec(4).work.sat.res_ph_by_ph(:,rec(4).work.system(idx_ph)' == 'R' & strLineMatch(rec(4).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec4 R' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'L2I'
[~,wl,idx_ph] = rec(7).work.getPhases(); idx_ph = find(idx_ph);
subplot(6,1,6); plot(rec(7).work.sat.res_ph_by_ph(:,rec(7).work.system(idx_ph)' == 'C' & strLineMatch(rec(7).work.obs_code(idx_ph,:),o_code)),'.'); ylim([-0.02 0.02]); xlim([0 720]); title(['Rec7 C' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);

figure; 
o_code = 'C1C'
[~,idx_pr]= rec(3).work.getPseudoRanges();idx_pr = find(idx_pr);

subplot(6,1,1); plot(rec(3).work.sat.res_pr_by_pr(:,rec(3).work.system(idx_pr)' == 'G' & strLineMatch(rec(3).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec3 G' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15); 
o_code = 'C2L'
subplot(6,1,2); plot(rec(3).work.sat.res_pr_by_pr(:,rec(3).work.system(idx_pr)' == 'G' & strLineMatch(rec(3).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec3 G' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15); 
o_code = 'C1C'
[~,idx_pr] = rec(5).work.getPseudoRanges(); idx_pr = find(idx_pr);

subplot(6,1,3); plot(rec(5).work.sat.res_pr_by_pr(:,rec(5).work.system(idx_pr)' == 'E' & strLineMatch(rec(5).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec5 E' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'C7Q'
subplot(6,1,4); plot(rec(5).work.sat.res_pr_by_pr(:,rec(5).work.system(idx_pr)' == 'E' & strLineMatch(rec(5).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec5 E' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'C1C'
[~,idx_pr] = rec(6).work.getPseudoRanges(); idx_pr = find(idx_pr);
subplot(6,1,5); plot(rec(6).work.sat.res_pr_by_pr(:,rec(6).work.system(idx_pr)' == 'R' & strLineMatch(rec(6).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec6 R' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);
o_code = 'C2I'
subplot(6,1,6); plot(rec(6).work.sat.res_pr_by_pr(:,rec(6).work.system(idx_pr)' == 'C' & strLineMatch(rec(6).work.obs_code(idx_pr,:),o_code)),'.'); ylim([-2 2]); xlim([0 720]); title(['Rec6 C' o_code]); ylabel('res [m]'); set(gca,'fontweight','bold','fontsize',15);