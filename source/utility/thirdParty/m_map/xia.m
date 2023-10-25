clc;
clear all;
%% T1:
load time
len=length(time);
for i=2:len
    dnb = datevec(time(i));   % datetime获取当天时间
    dna = datevec(time(i-1));   % 获取用户输入时间,格式如左边
    t=fix(etime(dnb,dna)); % 时间差,单位是秒
    if t~= 60*60
        str=datestr(time(i-1)+datenum('10000-00-00 01:00:00'),'yyyy-mm-dd HH:MM:SS');
       fprintf('缺失时间下标: %d, 时间为: %s\n',i,str);
    end
end

