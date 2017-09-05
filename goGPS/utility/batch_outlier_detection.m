function [outliers] = batch_outlier_detection(time_series, interval)

    outlier_thres = 1e-3;
    batch_size = 3600; %seconds

    num_batches = floor(length(time_series)*interval/batch_size);
    batch_idx = 1 : batch_size/interval : batch_size*num_batches/interval;

    if (num_batches == 0 && ~isempty(time_series))
        num_batches = 1;
        batch_idx = 1;
    end

    outliers = [];
    for b = 1 : num_batches
        start_idx = batch_idx(b);
        if (b == num_batches)
            end_idx = length(time_series);
        else
            end_idx = batch_idx(b+1)-1;
        end
        [~,~,outliers_batch] = deleteoutliers(time_series(start_idx:end_idx), outlier_thres);
        if (~isempty(outliers_batch))
            outliers = [outliers; outliers_batch(:)]; %#ok<AGROW>
        end
    end
end
