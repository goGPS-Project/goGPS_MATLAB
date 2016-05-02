function [ants,apc1,apc2,pcv1,pcv2,valid]=readantex2(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read antex antenna parameters file
% [func]   : read antex antenna parameters file
% [argin]  : file   = file path
% [argout] : ants   = receiver antenna type/satellite name
%            apc1   = L1 receiver anttena phase center offset  [n,e,u] (m)
%                     L1 satellite anttena phase center offset [x,y,z] (m)
%            apc2   = L2 receiver anttena phase center offset  [n,e,u] (m)
%                     L2 satellite anttena phase center offset [x,y,z] (m)
%            pcv1   = L1 receiver anttena pcv (el=90:-5:0deg) (m)
%                     L1 satellite anttena pcv (nadir=0:1:15deg) (m)
%            pcv2   = L2 receiver anttena pcv (el=90:-5:0deg) (m)
%                     L2 satellite anttena pcv (nadir=0:1:15deg) (m)
%            valid  = valid perod [start,end] (mjd)
% [note]   : no support of antenna azimuth angle dependency
% [version]: $Revision: 2 $ $Date: 06/07/08 14:19 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/04/22   0.1  new
%-------------------------------------------------------------------------------
ants={}; apc1=[]; apc2=[]; pcv1=[]; pcv2=[]; valid=[];

fd=fopen(file,'rt'); if fd<0, return, end

% read antex header
type=ReadHead(fd);

% read antex body
[ants,apc1,apc2,pcv1,pcv2,valid]=ReadData(fd);

fclose(fd);

% read antex header ------------------------------------------------------------
function type=ReadHead(fd)
type='';
while 1
    str=fgetl(fd); if ~isstr(str), break; end
    label=SubStr(str,61,20);
    if findstr(label,'ANTEX VERSION / SYST')
        type=SubStr(str,21,1);
    elseif findstr(label,'END OF HEADER'), break, end
end

% read rinex clock body --------------------------------------------------------
function [ants,apc1,apc2,pcv1,pcv2,valid]=ReadData(fd)
ants={}; azs=[]; els=[]; apc=[]; pcv=[]; valid=[]; n=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    label=SubStr(str,61,20);
    if findstr(label,'START OF ANTENNA')
        n=n+1;
        types{n}='';
        apc1(n,:)=zeros(1,3);
        apc2(n,:)=zeros(1,3);
        pcv1(n,:)=zeros(1,19);
        pcv2(n,:)=zeros(1,19);
        valid(n,:)=[0,inf];
        freq=0;
    elseif findstr(label,'TYPE / SERIAL NO')
        ants{n}=deblank(strrep(SubStr(str,1,20),'NONE','')); % omit radome-none
        if findstr(ants{n},'BLOCK')
            ants{n}=['GPS',deblank(SubStr(str,22,19))]; % satellite antenna
        end
%     elseif findstr(label,'VALID FROM')
%         valid(n,1)=caltomjd(StrToEpoch(str,1,43)');
%     elseif findstr(label,'VALID UNTIL')
%         valid(n,2)=caltomjd(StrToEpoch(str,1,43)');
    elseif findstr(label,'START OF FREQUENCY')
        freq=sscanf(SubStr(str,4,3),'G%d');
    elseif findstr(label,'END OF FREQUENCY')
        freq=0;
    elseif freq==1
        if findstr(label,'NORTH / EAST / UP')
            for m=1:3, apc1(n,m)=StrToNum(str,(m-1)*10+1,10)*1E-3; end
        elseif findstr(SubStr(str,1,8),'NOAZI')
            for m=1:19, pcv1(n,m)=StrToNum(str,(m-1)*8+9,8)*1E-3; end
        end
    elseif freq==2
        if findstr(label,'NORTH / EAST / UP')
            for m=1:3, apc2(n,m)=StrToNum(str,(m-1)*10+1,10)*1E-3; end
        elseif findstr(SubStr(str,1,8),'NOAZI')
            for m=1:19, pcv2(n,m)=StrToNum(str,(m-1)*8+9,8)*1E-3; end
        end
    end
end

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string to number -------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=str2num(str(pos:min(pos+ns-1,length(str)))); if isempty(num), num=0; end

% string to epoch --------------------------------------------------------------
function epoch=StrToEpoch(str,pos,ns)
epoch=sscanf(SubStr(str,pos,ns),'%d %d %d %d %d %f');
