%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE 2D TRAJECTORY on Google Maps
%----------------------------------------------------------------------------------------------

%if any positioning was done (either post-processing or real-time, not constrained)
if ((goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(lam_KAL)))
    %2D plot
    figure
    plot(lam_KAL, phi_KAL, '.r');
    xlabel('lon [deg]'); ylabel('lat [deg]');
    plot_google_map('MapType','hybrid');
end
