function [  fix_L2 ] = fix_jump_L2( til_L2,PRN,threshold )
%FIX_JUMP_L2 この関数の概要をここに記述
%   詳細説明をここに記述
fix_L2=til_L2(PRN,:);

    notnan_L2=fix_L2(~isnan(til_L2(PRN,:)));
    fix_idx=find(abs(diff(notnan_L2))>threshold);
    
    L2_fixed=notnan_L2;
    for i=1:length(fix_idx)
    L2_fixed(fix_idx(i)+1:end)=L2_fixed(fix_idx(i)+1:end)-(L2_fixed(fix_idx(i)+1)-L2_fixed(fix_idx(i)));
    end
    
    fix_L2(~isnan(til_L2(PRN,:)))=L2_fixed;


end

