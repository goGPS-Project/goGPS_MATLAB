if (flag_tropo & flag_tropo_gradient)

    idx = 1 : 10 : length(estim_tropo);

    u = Xhat_t_t_OUT(end-nC,idx);
    v = Xhat_t_t_OUT(end-nC-1,idx);

    quiver(datenum(date_R(idx,:)),zeros(1,length(idx)),u,v)
    grid on

    xlabel('hour-of-day')
    datetick('x','hh')
end