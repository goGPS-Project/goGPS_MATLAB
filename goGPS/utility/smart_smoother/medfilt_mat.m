% SYNTAX:
%    data_out = medfilt_mat(data_in, filter_size)
%
% DESCRIPTION:
%    Compute a moving window median to filter the data in input
%
function filtered_data = medfilt_mat(data, filter_size)
    % compute a moving window median to filter the data in input
    data = data(:);
    filtered_data = zeros(size(data));
    filter_size = filter_size + mod(filter_size+1,2); % filter size must be odd
    half_size = (filter_size-1) / 2;

    if (filter_size <= length(data))
        
        % Solution with increasing window at the left border
        win = data(1);
        filtered_data(1) = data(1);
        for i = 2 : half_size + 1
            win = inSortedData(win, data(i * 2 - 1));
            win = inSortedData(win, data((i-1) * 2));
            filtered_data(i) = win((length(win)+1)/2);
        end
        
        % Solution with constant median in the borders
%         % init the win with the first "part" of the dataset sorted
%         win = sort(data_in(1 : filter_size));
%         % get the median value
%         m = win(half_size + 1);
%         % first half window values will be equal to the median of the first window
%         filtered_data(1:half_size + 1) = m;
        
        for i = (half_size + 2) : (length(data) - half_size)
            % fprintf('Center: %d  In: %f(%2d)  Out %f(%2d)\n', i, data_in(i+half_size), i+half_size, data_in(i-half_size-1), i-half_size-1);
            win = inOutSortedData(win, data(i + half_size), data(i - half_size - 1));
            filtered_data(i) = win(half_size + 1);
        end
        
        % Solution with decreasing window at the right border
        for i = length(data) - half_size + 1 : length(data) -1
            win = outSortedData(win, data(-length(data) + 2 * (i - 1)));
            win = outSortedData(win, data(-length(data) + 2 * i - 1));
            filtered_data(i) = win((length(win)+1)/2);
        end
        filtered_data(end) = data(end);
        % Solution with constant median in the borders
%         % get the median value
%         m = win(half_size + 1);
%         % last half window values will be equal to the median of the first window
%         filtered_data(length(data_in) - half_size : end) = m;
    else
        error('The data to be filtered is shorter than the median window dimension')
    end
end

function data = inOutSortedData(data, in_key, out_key)
% insert a new key in data and remove another
% data = inOutSortedData(data, in_key, out_key)
% data must be a column array

    if (in_key ~= out_key)
        data(find(data == out_key, 1, 'first')) = []; % remove out_key
        insert_pos = find(data < in_key, 1, 'last');
        if isempty(insert_pos) % The key need to be inserted in the first position
            data = [in_key; data];
        else
            data = [data(1:insert_pos); in_key; data(insert_pos+1:end)];
        end
    end
end

function data = inSortedData(data, in_key)
% insert a new key in data
% data = inSortedData(data, in_key)
% data must be a column array

    insert_pos = find(data < in_key, 1, 'last');
    if isempty(insert_pos) % The key need to be inserted in the first position
        data = [in_key; data];
    else
        data = [data(1:insert_pos); in_key; data(insert_pos+1:end)];
    end
end

function data = outSortedData(data, out_key)
% remove a key in data
% data = outSortedData(data, out_key)
% data must be a column array
    data(find(data == out_key, 1, 'first')) = []; % remove out_key
end
