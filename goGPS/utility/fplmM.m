% Function for the computation of Legendre functions
% High Definition
% Matlab version
function [P] = fplmM(l, m, theta)
%%
    lMin = l(1);
    lMax = l(end);
    mMin = m(1);
    mMax = m(end);
    
    if (size(l,1)==1)
        l=l';
    end
    if (size(m,1)==1)
        m=m';
    end
    if (size(theta,1)==1)
        theta=theta';
    end
    
%%
    r1 = zeros(mMax,1);

    % computing root 1 (bl)
    r1(1) = sqrt(3);
    i = (2:mMax);
    l = (1:lMax);
    
    r1(2:end) = sqrt((2*i+1)./(2*i));
%%
    % Init P (result matrix)
    
    P = zeros(mMax,length(m),length(theta));

    % Computing Pmm

    % Skip the calculous of the first l(1)-1 lines
    % Go from 0 to lMin-1
    Ptmp0 = ones(length(theta),1);
    for mfix = 1:mMin-1
       Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
       Ptmp0 = Ptmp1;
    end
        
    % for each m
    r2 = zeros(length(l)-1,1);
    r3 = zeros(length(l)-1,1);
  
    for mfix = mMin:mMax
        % Computing Pmm --------------------------------------------------
        Ptmp1 = r1(mfix) * sin(theta).*Ptmp0;
        Ptmp0 = Ptmp1;
        
        % Save in the results matrix the Pmm element
        P(mfix, mfix-mMin+1, :) = Ptmp1(:);
        
        % Computing Plm --------------------------------------------------
       
        % get the row
        r = mfix+1;
        
        
        % computing root 2 (clm)
        r2 = sqrt(((2*l(r:end)+1).*(2*l(r:end)-1))./((l(r:end) + mfix).*(l(r:end) - mfix)));
        
        % computing root 3 (dlm)
        r3 = sqrt(((2*l(r:end)+1).*(l(r:end)+mfix-1).*(l(r:end)-mfix-1))./((2*l(r:end)-3).*(l(r:end)+mfix).*(l(r:end)-mfix)));
        
        Pl1 = Ptmp1;
        Pl2 = zeros(size(Ptmp1));
        for lfix = r:lMax
            tmp = r2(lfix-r+1) * cos(theta).*Pl1 - r3(lfix-r+1).*Pl2;
            P(lfix,mfix-mMin+1, :) = tmp(:);
            Pl2 = Pl1;
            Pl1 = tmp(:);
        end
    end;

    P = P(lMin:lMax,:,:);
 