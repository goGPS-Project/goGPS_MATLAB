%----------------------------------------------------------------------------------------------
% REPRESENTATION OF PSEUDORANGE DERIVATIVES
%----------------------------------------------------------------------------------------------

%master
for i = 1 : nSatTot
    index = find(pr1_M(i,:) ~= 0)';
    if ~isempty(index)
        pr = pr1_M(i,index);
        j = 1;
        diff_index = 0;
        while (j <= length(index) & diff_index ~= 1)
            diff_index = index(j+1) - index(j);
            interval = time_M(index(j+1)) - time_M(index(j));
            j = j + 1;
        end
        pr_der1 = zeros(length(index)-1,1);
        pr_der2 = zeros(length(index)-2,1);
        for j = 1 : length(index)-1
            if (index(j+1) == index(j) + 1)
                pr_der1(j) = (pr(j+1)-pr(j))/interval;
            else
                if (j > 1)
                    pr_der1(j) = pr_der1(j-1);
                end
            end
            if (j <= length(index)-2)
                if (index(j+2) == index(j) + 2)
                    pr_der2(j) = (pr(j+2)-2*pr(j+1)+pr(j))/interval^2;
                else
                    if (j > 1)
                        pr_der2(j) = pr_der2(j-1);
                    end
                end
            end
        end
        pos = find(pr_der1 == 0);
        pr_der1(pos) = [];
        pr_der2(pos) = [];
        index(pos) = [];
        figure
        plot(index(1:end-1), pr_der1,'b.');
        title(['MASTER: PSEUDORANGE FIRST DERIVATIVE for SATELLITE ',num2str(i)]);
        figure
        plot(index(1:end-2), pr_der2,'b.');
        title(['MASTER: PSEUDORANGE SECOND DERIVATIVE for SATELLITE ',num2str(i)]);
    end
end
