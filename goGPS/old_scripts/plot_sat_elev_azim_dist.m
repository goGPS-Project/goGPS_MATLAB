%----------------------------------------------------------------------------------------------
% REPRESENTATION OF AZIMUTH, ELEVATION AND DISTANCE FOR VISIBILE SATELLITES
%----------------------------------------------------------------------------------------------

% if (mode == goGNSS.MODE_PP_KF_CP_DD)

   coltab = jet(2*nSatTot);

   f1 = figure; hold on; grid on; title('Azimuth')
   f2 = figure; hold on; grid on; title('Elevation')
   f3 = figure; hold on; grid on; title('Distance')
   k = 1;
   for i = 1 : nSatTot
      index = find(abs(conf_sat_OUT(i,:)) == 1)';
      if ~isempty(index)
         %azimuth
         figure(f1)
         h = plot(index,azR(i,index),'b.'); grid on;
         set(h,'Color',coltab(2*i-1,:));
         %elevation
         figure(f2)
         h = plot(index,elR(i,index),'r.');
         set(h,'Color',coltab(2*i-1,:));
         %distance
         figure(f3)
         h = plot(index,distR(i,index)*1e-6,'g.');
         set(h,'Color',coltab(2*i-1,:));
         %legend
         list{k} = num2str(i);
         k = k+1;
      end
   end
   figure(f1); hold off; legend(list)
   figure(f2); hold off; legend(list)
   figure(f3); hold off; legend(list)
   clear f1 f2 f3
   clear i k h
   clear coltab

% end
