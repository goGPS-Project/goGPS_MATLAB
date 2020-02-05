[Y M D h m s lat lon] = textread('tgf9RAD.txt', '%f %f %f %f %f %f %f %f');
crd_list=[lat lon Y M D h m s];
%coord_list = [];
dist=[];
receiver=[];
year=[];
month=[];
day=[];
hour=[];
minute=[];
second=[];
for r = 1 : length(crd_list)
    r
    [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(crd_list(r,1), crd_list(r,2), 20);
     %coord_list = [coord_list; coo.xyz(1,1), coo.xyz(1,2)];
     d_km = deg2km(d);
     for s = 1 : numel(d)
         if d_km(s) < 45
             dist=[dist;d_km(s)];
             receiver=[receiver; marker_list(s)];
             year=[year; crd_list(r,3)];
             month=[month; crd_list(r,4)];
             day=[day; crd_list(r,5)];
             hour=[hour; crd_list(r,6)];
             minute=[minute; crd_list(r,7)];
             second=[second; crd_list(r,8)];
         end
     end
end
out=[receiver num2cell(dist) num2cell(year) num2cell(month) num2cell(day) num2cell(hour) num2cell(minute) num2cell(second)];