[WMO lat lon alt] = textread('radiosonde.txt', '%f %f %f %f');
[IGS]=textread('igs_list.txt', '%s');
igs_list=[IGS]
crd_list=[lat lon alt WMO];
%coord_list = [];
dist=[];
receiver=[];
WMO=[];
lat=[];
lon=[];
alt=[];
Lia=[];
for r = 1 : length(crd_list)
    r
    [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(crd_list(r,1), crd_list(r,2), 20);
     %coord_list = [coord_list; coo.xyz(1,1), coo.xyz(1,2)];
     d_km = deg2km(d);
     for s = 1 : numel(d)
         if d_km(s) < 15
             dist=[dist;d_km(s)];
             receiver=[receiver; marker_list(s)];
             WMO=[WMO; crd_list(r,4)];
             lat=[lat; crd_list(r,1)];
             lon=[lon; crd_list(r,2)];
             alt=[alt; crd_list(r,3)];
             Lia = ismember(receiver,igs_list);
         end
     end
end
out=[receiver num2cell(dist) num2cell(WMO) num2cell(lat) num2cell(lon) num2cell(alt) num2cell(Lia)];