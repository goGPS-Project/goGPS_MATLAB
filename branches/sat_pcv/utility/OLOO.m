function [imax,xfin,s2fin,ufin,Cxx,uout]=OLOO(A, y, Q)
%Purpose:   perform LS on blocks of correlated observations
%           identify one (block) outlier
%           reject it
%           re-estimate unknowns
%           according to the theory in "L. Biagi and S. Caldera.An efficient leave one block out approach to identify outliers.Journal of Applied Geodesy, Volume 7, Issue 1, pages 11–19, 2013"
%1.0: Stefano Caldera, 22.05.2014
%input:
%   A: design matrix
%   y: observations vector
%   Q: cofactor matrix
%output
%   imax: index of the rejected blocks
%   x_fin: estimated parameters without outlier
%   s2fin: a posteriori sigma without outlier
%   ufin: estimated residuals without outlier
%   Cxx: parameters covariance of the final solution
%   uout: residual of outlier observation

% this version is optimized to manage l.o.o. of 1 observation at time, no
% blocks!


%%
global FTABLE;

n_blocks=length(y);
[m,n]=size(A);           % m: number of observations, n: number of unknowns
uout=NaN;
if m - n > 0
   
    %% compute the global solution
    
    Q=Q./(min(diag(Q)));
    invQ=inv(Q);
    Ninv=inv(A'*invQ*A);
    xcap=Ninv*A'*invQ*y;
    um=y-A*xcap;
    s2cap=um'*invQ*um/(m-n);
    
    Qm=diag(Q);
    
    if m - n > 1
        %% start outliers rejection
        Im=eye(n);
        Bm=Ninv*A';
        Cm=diag(A*Bm);
        Km=Qm-Cm;
        Kminv=Km.^-1;
        wm=Qm.*Kminv.*um;
        s2m=((m-n).*s2cap-um.*Kminv.*um)./(m-n-1);      

%         keyboard
        
%         clear temp
%         tic               
%         temp=[];
%         for i=1:n_blocks
%            temp= blkdiag(temp,Bm(:,i)');             
%         end
%         toc
%         %[Bm',zeros((size(Bm',1)-1)*size(Bm',2),size(Bm',1))'];
%         
%         tic
%         D=reshape([Bm',zeros((size(Bm',1))*size(Bm',2),size(Bm',1))']',[],1)';
%         D=reshape(D(1:end-(size(Bm',1))*size(Bm',2)),(size(Bm',1))*size(Bm',2),size(Bm,2))';
%         Qw1=D*(eye(n*m)+reshape(A',[],1).*reshape(repmat(Kminv,1,4)',[],1)*Bm(:)');
%         toc
%         
%         %Qw2=kron(eye(n_blocks),Bm')*(eye(n*m)+reshape(A',[],1).*reshape(repmat(Kminv,1,4)',[],1)*Bm(:)');
%         %toc
% tic
        Qw=zeros(n_blocks,1);
        for i=1:n_blocks
            Qw(i)=Qm(i)+Bm(:,i)'*(Im+A(i,:)'*Kminv(i)*Bm(:,i)')*A(i,:)';
            %Qw3=Bm(:,i)'*(Im+A(i,:)'*Kminv(i)*Bm(:,i)');
        end
%toc
        Qwinv=Qw.^-1;
        deg2=m-n-1;
        
        F=wm.*Qwinv.*wm./s2m;        
        Flim=FTABLE(deg2);
        
        
        %% apply final solution
        % find maximum F(i)/Flim(i)
        [Fmax,imax]=max(abs(F./Flim));

        if (Fmax<1)
            % no outlier
            imax=0;
            xfin=xcap;
            %yfin=y;
            %Afin=A;
            %Qfin=Q;
            s2fin=s2cap;
            Ninvfin=Ninv;
            ufin=um;            
            Cxx=s2fin*Ninvfin;
            
        else
            % if the maximum ratio exceedes the threshold, the observation is eliminated from the solution
            uout=um(imax);
            xfin=xcap-Bm(:,imax)*Kminv(imax)*um(imax);
            yfin=y;
            yfin(imax)=[];

            Afin=A;
            Afin(imax,:)=[];
            
            %Qfin=Q;
            %Qfin(imax,:)=[];
            %Qfin(:,imax)=[];
            
            s2fin=s2m(imax);
            
            Ninvfin=Ninv+Bm(:,imax)*Kminv(imax)*Bm(:,imax)';
            ufin=yfin-Afin*xfin;
            Cxx=s2fin*Ninvfin;

        end

    else
        imax=0;
        xfin=xcap;
        %yfin=y;
        %Afin=A;
        %Qfin=Q;
        s2fin=s2cap;
        Cxx=s2cap*Ninv;
        ufin=um;
    end
    
else
    % detection is not possibile
    imax=0;
    xfin=[];
    Cxx=[];
    %yfin=[];
    ufin=[];
    s2fin=NaN;
end

