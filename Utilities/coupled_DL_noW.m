function [Alphap, Alphas, Xp, Xs, Dp, Ds,f] = coupled_DL_noW(Alphap, Alphas, Xp, Xs, Dp, Ds,  par)
% Semi-Coupled Dictionary Learning
% Shenlong Wang
% Reference: S Wang, L Zhang, Y Liang and Q. Pan, "Semi-coupled Dictionary Learning with Applications in Super-resolution and Photo-sketch Synthesis", CVPR 2012

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
    Alphas = mexLasso(Xs, Ds,param);
    Alphap = mexLasso(Xp, Dp,param);
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
    

    % Alpha = pinv(D' * D + lambda2 * eye(numD)) * D' * X;
    P1 = Xp - Dp * Alphap;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphap, 1);    
    P3 = Alphas -  Alphap;
    P3 = P3(:)'*P3(:) / 2;

    fp = 1 / 2 * P1 + P2 + mu * (P3 );
    
    P1 = Xs - Ds * Alphas;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphas, 1);    
    P3 = Alphap - Alphas;
    P3 = P3(:)'*P3(:) / 2;

    fs = 1 / 2 * P1 + P2 + mu * (P3 );
    f = fp + fs;
    time=toc;
    fprintf('%d-th: Energy: %d Time: %d\n',t,f,time);
    if (abs(f_prev - f) / f < epsilon)
        break;
    end

%     save tempDict_SR_NL Ds Dp Ws Wp par param i;
    % fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4), E);
end



    