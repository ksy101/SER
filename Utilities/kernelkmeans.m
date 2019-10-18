function [mask_t] = kernelkmeans(X,Y,num_k,mask_kt, gam,lambda)

% X = YLt in R{5 * N}
N = size(X,2);
Ke1 = zeros(N);
Ke2 = zeros(N);
mask_t = zeros(size(mask_kt));
gam1 = 0.5/gam;
% gam2 = 0.5/gam(2);
gam2 = gam1;
% lambda = 0.7;
for n = 1:N
    for n2 = 1:N
        Ke1(n,n2) = exp(-gam1*sum(sum((X(:,n)-X(:,n2)).^2)));
        Ke2(n,n2) = exp(-gam2*sum(sum((Y(:,n)-Y(:,n2)).^2)));
    end
end
Ke = lambda*Ke1 + (1-lambda)*Ke2;
%% initialization
converged = 0;
% Assign all objects into one cluster except one
% Kernel K-means is *very* sensitive to initial conditions.  Try altering
% this initialisation to see the effect.

di = zeros(N,num_k);
mask_k1=zeros(num_k,N);

for k = 1:num_k
    mask_k1(k,mask_kt==k)=1;  
end


while ~converged
    
    for k = 1:num_k
%         mask_k1=zeros(size(mask_kt));
%         mask_k1(k,mask_kt==k)=1;  
        Nk=sum(mask_k1(k,:)==1);
        % Compute kernelised distance
        di(:,k) = diag(Ke) - (2/(Nk))*sum(repmat(mask_k1(k,:),N,1).*Ke,2) + ...
            Nk^(-2)*sum(sum((mask_k1(k,:)'*mask_k1(k,:)).*Ke));
    end
    mask_k1_old = mask_k1;
    mask_k1 = (di == repmat(min(di,[],2),1,num_k))';
    mask_k1 = 1.0*mask_k1;
    if sum(sum(mask_k1_old~=mask_k1))==0
        converged = 1;
    end
end
for k = 1:num_k
    mask_t(mask_k1(k,:)==1)=k;  
end


