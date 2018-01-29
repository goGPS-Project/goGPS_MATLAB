%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scmidt quasi normalized associated functions discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 15;

% co to colatitude
theta = 0:0.001:1;
P = zeros(n,n,length(theta));
for j = 1:length(theta)
    for i = 0:n-1;
        P(1:i+1,i+1,j) = legendre(i,cos(theta(j)),'sch');
    end
end