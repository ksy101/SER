function [Alphap, Alphas] = coupled_DL_test(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)


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
% t0                  =       par.t0;
epsilon             =       par.epsilon;
param.lambda       = lambda1; % not more than 20 non-zeros coefficients
param.lambda2       = lambda2;
param.mode          = 2;       % penalized formulation
param.approx=0;
param.K = par.K;
% param.L = par.L;
f = 0;


for t = 1 : nIter
    % Alphat = mexLasso(Xt,D,param);
    tic
    f_prev = f;
    
    if numX == 65536*11
        Xs_tal = [Xs;sqrtmu * full(Alphap)];    % size [178 655360]
        Xp_tal = [Xp;sqrtmu * full(Alphas)];    % size [178 655360]
        Alphas_sgl = zeros(par.K,65536,11);
        Alphap_sgl = zeros(par.K,65536,11);
        for num_P = 1:11
            Xs_sgl = Xs_tal(:,(num_P-1)*65536+1:num_P*65536);
            Xp_sgl = Xp_tal(:,(num_P-1)*65536+1:num_P*65536);        
            Alphas_sgl(:,:,num_P) = mexLasso(Xs_sgl,[Ds; sqrtmu * Ws] ,param); % size [128 65536]
            Alphap_sgl(:,:,num_P) = mexLasso(Xp_sgl, [Dp; sqrtmu * Wp],param);
        end
        Alphas = reshape(Alphas_sgl,par.K,[]);                        % size [128 655360]
        Alphap = reshape(Alphap_sgl,par.K,[]);
        
    else     
        Alphas = mexLasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws],param);
        Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);
    end
    
    Xp = Dp*Alphap;




    P1 = Xp - Dp * Alphap;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphap, 1);    
    P3 = Alphas - Wp * Alphap;
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Wp, 'fro');
    fp = 1 / 2 * P1 + P2 + mu * (P3 + P4);
    
    P1 = Xs - Ds * Alphas;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphas, 1);    
    P3 = Alphap - Ws * Alphas;  
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Ws, 'fro'); 
    fs = 1 / 2 * P1 + P2 + mu * (P3 + P4);
    f = fp + fs;
 
    
    time=toc;
    
    fprintf('%d-th: Energy: %d; Time: %d\n',t,f,time);
    
    if (abs(f_prev - f) / f < epsilon)
%         break;
        epsilon = epsilon/10;
        param.lambda = param.lambda/10;
    end

%     save tempDict_SR_NL Ds Dp Ws Wp par param i;
    % fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4), E);
end



    