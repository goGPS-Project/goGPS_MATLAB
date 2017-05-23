%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE S/N RATIO FOR MASTER AND ROVER
%----------------------------------------------------------------------------------------------

if (mode == goGNSS.MODE_PP_KF_CP_DD)

   coltab = jet;
   coltab = [1 1 1; coltab([1 16 32 48 56],:)];

   figure
   subplot(2,1,1);
   imagesc(snr_R.*abs(conf_sat_OUT),[-0.5 60.5]);
   title('Rover S/N ratio');
   axis xy; colormap(coltab);
   h = colorbar; set(h,'YTick',0:10:60);
   subplot(2,1,2);
   imagesc(snr_M.*abs(conf_sat_OUT),[-0.5 60.5]);
   title('Master S/N ratio');
   axis xy; colormap(coltab);
   h = colorbar; set(h,'YTick',0:10:60);

   clear h coltab
end
