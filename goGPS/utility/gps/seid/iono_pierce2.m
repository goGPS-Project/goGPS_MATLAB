function [posp,z]=iono_pierce2(azel,posr,H)
RE=6370000; 
zr=pi/2-azel(2); z=asin(RE*sin(zr)/(RE+H));
posp(1)=asin(cos(zr-z)*sin(posr(1))+sin(zr-z)*cos(posr(1))*cos(azel(1)));
posp(2)=posr(2)+asin(sin(zr-z)*sin(azel(1)/cos(posp(1))));

%  for s = 1 : 32
% pos = find(azR(s,:)>180);
% azR(s,pos) = azR(s,pos) - 360;
% azel(s,:,1)=azR(s,:);
% azel(s,:,2)=elR(s,:);
% end
% for s = 1 : 32
% for i = 1 : size(azel,2)
% [posp(s,i,:)]=iono_pierce(azel(s,i,:).*pi/180,posr(1,:));
% end
% end
% figure
% hold on
% for s =  1 : 32
% plot(posp(s,:,2)*180/pi,posp(s,:,1)*180/pi,'.');
% end