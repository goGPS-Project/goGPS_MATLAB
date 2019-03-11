%----------------------------------------------------------------------------------------------
% REPRESENTATION OF PHASE RANGE DERIVATIVES
%----------------------------------------------------------------------------------------------

%rover
for i = 1 : nSatTot
    index = find(ph1_R(i,:) ~= 0)';
    if ~isempty(index)
        ph = ph1_R(i,index);
        j = 1;
        diff_index = 0;
        while (j <= length(index) & diff_index ~= 1)
            diff_index = index(j+1) - index(j);
            interval = time_R(index(j+1)) - time_R(index(j));
            j = j + 1;
        end
        ph_der1 = zeros(length(index)-1,1);
        ph_der2 = zeros(length(index)-2,1);
        for j = 1 : length(index)-1
            if (index(j+1) == index(j) + 1)
                ph_der1(j) = (ph(j+1)-ph(j))/interval;
            else
                if (j > 1)
                    ph_der1(j) = ph_der1(j-1);
                end
            end
            if (j <= length(index)-2)
                if (index(j+2) == index(j) + 2)
                    ph_der2(j) = (ph(j+2)-2*ph(j+1)+ph(j))/interval^2;
                else
                    if (j > 1)
                        ph_der2(j) = ph_der2(j-1);
                    end
                end
            end
        end
        pos = find(ph_der1 == 0);
        ph_der1(pos) = [];
        ph_der2(pos) = [];
        index(pos) = [];
        figure
        plot(index(1:end-1), ph_der1,'b-');
        title(['ROVER: PHASE FIRST DERIVATIVE for SATELLITE ',num2str(i)]);
        figure
        plot(index(1:end-2), ph_der2,'b-');
        title(['ROVER: PHASE SECOND DERIVATIVE for SATELLITE ',num2str(i)]);
    end
end

%master
for i = 1 : nSatTot
    index = find(ph1_M(i,:) ~= 0)';
    if ~isempty(index)
        ph = ph1_M(i,index);
        j = 1;
        diff_index = 0;
        while (j <= length(index) & diff_index ~= 1)
            diff_index = index(j+1) - index(j);
            interval = time_M(index(j+1)) - time_M(index(j));
            j = j + 1;
        end
        ph_der1 = zeros(length(index)-1,1);
        ph_der2 = zeros(length(index)-2,1);
        for j = 1 : length(index)-1
            if (index(j+1) == index(j) + 1)
                ph_der1(j) = (ph(j+1)-ph(j))/interval;
            else
                if (j > 1)
                    ph_der1(j) = ph_der1(j-1);
                end
            end
            if (j <= length(index)-2)
                if (index(j+2) == index(j) + 2)
                    ph_der2(j) = (ph(j+2)-2*ph(j+1)+ph(j))/interval^2;
                else
                    if (j > 1)
                        ph_der2(j) = ph_der2(j-1);
                    end
                end
            end
        end
        pos = find(ph_der1 == 0);
        ph_der1(pos) = [];
        ph_der2(pos) = [];
        index(pos) = [];
        figure
        plot(index(1:end-1), ph_der1,'b-');
        title(['MASTER: PHASE FIRST DERIVATIVE for SATELLITE ',num2str(i)]);
        figure
        plot(index(1:end-2), ph_der2,'b-');
        title(['MASTER: PHASE SECOND DERIVATIVE for SATELLITE ',num2str(i)]);
    end
end
