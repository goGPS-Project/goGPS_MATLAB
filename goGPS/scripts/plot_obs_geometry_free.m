if (length(frequencies) > 1)
    %geometry free (rover)
    for i = 1 : nSatTot
        index_pr1 = find(pr1_R(i,:) ~= 0)';
        index_ph1 = find(ph1_R(i,:) ~= 0)';
        index1 = intersect(index_pr1, index_ph1);
        index_pr2 = find(pr2_R(i,:) ~= 0)';
        index_ph2 = find(ph2_R(i,:) ~= 0)';
        index2 = intersect(index_pr2, index_ph2);
        index = intersect(index1, index2);
        if ~isempty(index)
            ph_GF = compute_geometry_free(ph1_R, ph2_R, lambda, zeros(size(ph1_R)));
            figure
            plot(index,ph_GF(i,index),'b.-'); grid on;
            title(['Geometry free for SATELLITE ',num2str(i)]);
            xlabel(['epoch [' num2str(interval) ' sec]']);
            ylabel('[m]');
        end
    end
end
