function [Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws,n,f] = coupled_DL_n(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)
% Based on Semi-Coupled Dictionary Learning <Wang2012>


[dimX, numX]        =       size(Xp);
dimY                =       size(Alphap, 1);
numD                =       size(Dp, 2);
rho                 =       par.rho;
lambda1             =       par.lambda1;
lambda2             =       par.lambda2;
mu                  =       par.mu;
sqrtmu              =       sqrt(mu);
nu                  =       par.nu;
nIter               =       par.nIter;
NUM_P               =       par.NUM_P;
% t0                  =       par.t0;
epsilon             =       par.epsilon;
param.lambda       = lambda1; % not more than 20 non-zeros coefficients
param.lambda2       = lambda2;
param.mode          = 2;       % penalized formulation
param.approx=0;
param.K = par.K;
% param.L = par.L;
f = 0;
n_ran=zeros(dimY,1);
%         Xs_tal = [Xs];    % size [178 655360]
%         Xp_tal = [Xp];    % size [178 655360]
        
for t = 1 : nIter
    % Alphat = mexLasso(Xt,D,param);
    tic
    f_prev = f;
    
    if numX == 65536*NUM_P


        Alphas_sgl = zeros(par.K,65536,NUM_P);
        Alphap_sgl = zeros(par.K,65536,NUM_P);
        for num_P = 1:NUM_P
            Xs_sgl = Xs(:,(num_P-1)*65536+1:num_P*65536);
            Xp_sgl = Xp(:,(num_P-1)*65536+1:num_P*65536);        
%             Alphas_sgl(:,:,num_P) = mexLasso(Xs_sgl,[Ds; sqrtmu * Ws] ,param); % size [128 65536]
            Alphas_sgl(:,:,num_P) = mexLasso(Xs_sgl, Ds, param); % broke the chain
            Alphap_sgl(:,:,num_P) = mexLasso([Xp_sgl;sqrtmu * full(Alphas_sgl(:,:,num_P))], [Dp; sqrtmu * Wp], param);
        end

        Alphas = reshape(Alphas_sgl,par.K,[]);                        % size [128 655360]
        
        Alphap = reshape(Alphap_sgl,par.K,[]);
        
    else     
%         Alphas = mexLasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws],param);
        Alphas = mexLasso(Xs, Ds, param); % broke the chain
%         Xp_D = Dp*Ws*Alphas;
%         Alphas = mexLasso([Xs;sqrtmu*Xp_D], [Ds;sqrtmu*Dp*Ws],param); % broke the chain
        Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);
    end
%  

    Alphas = gpuArray(full(Alphas));
    Alphap = gpuArray(full(Alphap));    
    Xs = gpuArray(Xs);
    Ds = gpuArray(Ds);
    Xp = gpuArray(Xp);
    Dp = gpuArray(Dp);
    Wp = gpuArray(Wp);
    Ws = gpuArray(Ws);
    
    Alphas_sum = sum(Alphas,1);
    idx_n = zeros(size(idx_n),'gpuArray');
    idx_n(Alphas_sum~=0)=1;
    Alphas(:,idx_n)+n;
    
dictSize = par.K;


    % Update D    
    for i=1:dictSize
       ai        =    Alphas(i,:);
       Y         =    Xs-Ds*Alphas+Ds(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Ds(:,i)    =    di;
    end
    
%     Ds = gather(Ds);

    for i=1:dictSize
       ai        =    Alphap(i,:);
       Y         =    Xp-Dp*Alphap+Dp(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Dp(:,i)    =    di;
    end
    
%     Update W
%     Ws = Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;
%     Wp = Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;    
    Ws = (1 - rho) * Ws  + rho * Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;
    Wp = (1 - rho) * Wp  + rho * Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;
    % Alpha = pinv(D' * D + lambda2 * eye(numD)) * D' * X;
    P1 = Xp - Dp * Alphap;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphap, 1);    
    P3 = Alphas - Wp * Alphap;
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Wp, 'fro');
    fp = 1 / 2 * P1 + P2 + mu * (P3 + P4)
    
    P1 = Xs - Ds * Alphas;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphas, 1);    
    P3 = Alphap - Ws * Alphas;
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Ws, 'fro'); 
    fs = 1 / 2 * P1 + P2 + mu * (P3 + P4)
    f = fp + fs;

    
    Alphas = gather(Alphas);
    Xs = gather(Xs);
    Ds = gather(Ds);
    Alphap = gather(Alphap);
    Xp = gather(Xp);
    Dp = gather(Dp);
    Wp = gather(Wp);
    Ws = gather(Ws);
    f = gather(f);

    
    time=toc;
    fprintf('%d-th: Energy: %d Time: %d\n',t,f,time);
    save tempDict_SR_NL Ds Dp Ws Wp par param f time t;
    % fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4), E);
    
    if (abs(f_prev - f) / f < epsilon)
        break;
    end
end




    