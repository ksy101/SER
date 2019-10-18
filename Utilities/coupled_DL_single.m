function [Alphap, Alphas, XP, XS, Dp, Ds, Wp, Ws,f] = coupled_DL( Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, siz,par,param,params)
% Semi-Coupled Dictionary Learning
% Shenlong Wang
% Reference: S Wang, L Zhang, Y Liang and Q. Pan, "Semi-coupled Dictionary Learning with Applications in Super-resolution and Photo-sketch Synthesis", CVPR 2012

[dimX, numX]        =       size(Xp);
% dimY                =       size(Alphap, 1);
dictSize                =       size(Dp, 2);
rho                 =       par.rho;
lambda1             =       par.lambda1;
% lambda2             =       par.lambda2;
% lambda3             =       par.lambda3;
mu                  =       par.mu;
sqrtmu              =       sqrt(mu);
nu                  =       par.nu;
nIter               =       par.nIter;
% t0                  =       par.t0;
epsilon             =       par.epsilon;
% param.lambda       = 1e-5; % not more than 20 non-zeros coefficients
% param.lambda2       = lambda2;
% param.mode          = 2;       % penalized formulation
% param.approx=0;
% % param.K = par.K;
% param.L = par.L;
f = 0;
% Xt=[Xp;Xs];
% Dt=[Dp;Ds];
% Alphap0 = Alphap;
% Alphas0 = Alphas;
% Ws1=Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1)));
% Wp1=Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1)));

for t = 1 : nIter

    f_prev = f;
%     Alphas = mexLasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws],params);
    Alphas = mexLasso(Xs, Ds, params);
%     Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);
    Alphap = mexLasso(Xp, Dp, param);
% 
%     Alphas = L1Lasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws]);
%     Alphap = L1Lasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp]);
%     Alphas = L1Lasso(Xs, Ds,lambda3);
%     Alphap = L1Lasso(Xp, Dp,lambda3);  

%     Alphat = mexLasso(Xt, Dt, param);
    

    
    % Update D
    for i=1:dictSize
       ai        =    Alphas(i,:);
       Y         =    Xs-Ds*Alphas+Ds(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Ds(:,i)    =    di;
    end
    for i=1:dictSize
       ai        =    Alphap(i,:);
       Y         =    Xp-Dp*Alphap+Dp(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Dp(:,i)    =    di;
    end
%     
%     for i=1:size(Dt,2)
%        ai        =    Alphat(i,:);
%        Y         =    Xt-Dt*Alphat+Dt(:,i)*ai;
%        di        =    Y*ai';
%        di        =    di./(norm(di,2) + eps);
%        Dt(:,i)    =    di;
%     end

% Alphap=Alphat;
% Alphas=Alphat;
% Dp = Dt(1:50,:);
% Ds = Dt(51:end,:);
    % Update W
%     Ws = Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;
%     Wp = Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;    
    
    Wp = (1 - rho) * Wp  + rho * Alphas * Alphap' * pinv(Alphap * Alphap' + nu * eye(size(Alphap, 1))) ;
    Ws = (1 - rho) * Ws  + rho * Alphap * Alphas' * pinv(Alphas * Alphas' + nu * eye(size(Alphas, 1))) ;
    
%     Ws = l2ls_learn_basis_dual(Alphap, Alphas, 1);
%     Wp = l2ls_learn_basis_dual(Alphas, Alphap, 1);
%     Ws = (1 - rho) * Ws  + rho *  Ws1;
%     Wp = (1 - rho) * Wp  + rho *  Wp1;


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
%     if (abs(f_prev - f) / f < epsilon)
%         break;
%     end

%     fprintf('Energy: %d\n',f);
%     save tempDict_SR_NL Ds Dp Ws Wp par param i;
%  YH=reshape((Dp * Alphap)',50,50,60);
%  YL=reshape((Ds * Alphas)',50,50,60);
 
%  GT = reshape(Xp',siz);
 XP = reshape((Dp * Alphap)',siz);
 XS = reshape((Dp * Ws * Alphas)',siz);
%  E=error_dic(GT,XS);
E=norm(Ws*Alphas-Alphap);
Ep=norm(Dp*Alphap-Xp);
Es=norm(Ds*Alphas-Xs);
Sp=full(max(sum(Alphap~=0,1)));
Ss=full(max(sum(Alphas~=0,1)));

%     fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4),E);
    fprintf('Iter: %d, Eh : %d, El : %d, E : %d\n', t, Ep, Es,E);
    fprintf('Iter: %d, Sh : %d, Sl : %d\n', t, Sp,Ss);
%         if E < 10
%         break;
%         end
end
    