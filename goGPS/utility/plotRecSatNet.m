function plotRecSatNet(r_id,s_id)
% show the connection between satellite and receiver at that epoch
u_r = unique(r_id);
u_s = unique(s_id);
if length(u_r) > length(u_s)
    s_r = 1;
    s_s = length(u_r)/length(u_s);
else
    s_r = length(u_s)/length(u_r);
    s_s = 1;
end
figure;
r_pos = s_r:s_r:(length(u_r)*s_r);
scatter((s_r:s_r:(length(u_r)*s_r))',zeros(size(u_r)),'^');
hold on
cc = Core.getConstellationCollector();
a_id = cc.getAntennaId(u_s);
p_id = 0;
s_pos = [];
for s = unique(a_id(:,1))'
    idx_s = a_id(:,1) == s;
    if s == 'G'
        symbl = 'o';
    elseif s == 'R'
        symbl = '+';
    elseif s == 'E'
        symbl = '*';
    elseif s == 'C'
        symbl = 'x';
    elseif s == 'J'
        symbl = 's';
    elseif s == 'I'
        symbl = 'd';
    elseif s == 'S'
        symbl = 'v';
    end
    s_pos = [s_pos (p_id +(s_s:s_s:(sum(idx_s)*s_s)))];
    scatter((p_id +(s_s:s_s:(sum(idx_s)*s_s))'),ones(sum(idx_s),1),symbl);  
    p_id = p_id + sum(idx_s)*s_s;
end
sat2pos = zeros(max(u_s),1);
for s = 1:size(a_id,1)
    text(s_pos(s),1.1,a_id(s,:),'HorizontalAlignment','center');
end
sat2pos(u_s) = s_pos;
rec2pos = zeros(max(u_r),1);
rec2pos(u_r) = r_pos;
for o = 1:length(r_id)
    line([rec2pos(r_id(o)) sat2pos(s_id(o))], [0 1]);
end
ylim([0 1.2])



end