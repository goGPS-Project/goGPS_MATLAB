%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE VISIBLE SATELLITES CONFIGURATION
%----------------------------------------------------------------------------------------------

% if (mode == goGNSS.MODE_PP_KF_CP_DD)

   %figure
   %imagesc(abs(conf_sat_OUT)), grid;
   %colormap(1-gray);

   figure
   subplot('position',[0.1 0.35 0.8 0.55]);
   hold on; grid on;
   title('Satellite configuration')
   for i = 1 : nSatTot
      index = find(abs(conf_sat_OUT(i,:)) == 1);
      index_cs = intersect(index, find(conf_cs(i,:) == 1));
      index_pivot = intersect(index, find(pivot_OUT == i));
      if ~isempty(index)
         plot(index,i*ones(size(index)),'b.-');
         plot(index_pivot,i*ones(size(index_pivot)),'r.');
         plot(index_cs,i*ones(size(index_cs)),'g.');
      end
   end
   axis([1 size(conf_sat_OUT,2) 0.5 32.5])
   hold off;
   clear i index index_cs index_pivot

   subplot('position',[0.1 0.1 0.8 0.2]);
   hold on; grid on;
   s1 = sum(abs(conf_sat_OUT));
   plot(s1,'b.-');
   s2 = [0; pivot_OUT(2:end) - pivot_OUT(1:end-1)];
   plot(find(s2>0),s1(s2>0),'r.')
   s3 = sum(conf_cs);
   plot(find(s3>0),s3(s3>0),'g.');
   axis([1 size(conf_sat_OUT,2) 0 max(s1)])
   hold off;
   clear s1 s2 s3

% end
