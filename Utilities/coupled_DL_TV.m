function [Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws,f] = coupled_DL_TV(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)
% Based on Semi-Coupled Dictionary Learning <Wang2012>


[dimX, numX]        =       size(Xs);
height              =       sqrt(numX);
width              =       sqrt(numX);
dimY                =       size(Alphap, 1);
numD                =       size(Dp, 2);
rho                 =       par.rho;
lambda1             =       par.lambda1;
lambda2             =       par.lambda2;
lambda3             =       par.lambda3;
mu                  =       par.mu;
sqrtmu              =       par.sqrtmu;
nu                  =       par.nu;
nIter               =       par.nIter;
NUM_P               =       par.NUM_P;
nIter_P             =       par.nIter_P;
Tw             =       par.Tw;
% t0                  =       par.t0;
epsilon             =       par.epsilon;
param.lambda       = lambda1; % not more than 20 non-zeros coefficients
param.lambda2       = lambda2;
param.mode          = 2;       % penalized formulation
param.approx=0;
param.K = par.K;
% param.L = par.L;
f = 0;
% nn=0;
% beta = 1;
% z = zeros(size(Alphas));
% d = zeros(size(z));
paramp = param;
paramp.lambda =  lambda3;
nIter_W = nIter_P;
dictSize = par.K;

% % Dp_pre
paramD.K = par.K;
paramD.iter=300; 
paramD.mode=2;
paramD.lambda  = lambda3; % 1e-1    
Dp_pre = mexTrainDL(Xp, paramD);

% % TV
  L2=8.0;
  tau=par.tau;
  theta = par.theta;
  sigma=1.0/(L2*tau);
  p=zeros(height, width, 2*dimX,'gpuArray');

for t = 1 : nIter
    
    
    tic
    f_prev = f;
%% update sparsity coefficients   
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
     
        % % initial
%         Alphas = mexLasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws],param);


        % % broke the chain 
        % % v1: || Alphas - Wp * Alphap ||
%         Alphas = mexLasso(Xs, Ds, param);
%         Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);


        % % v2: || Alphap - Ws * Alphas ||
        Alphas = mexLasso(Xs, Ds, param);
        Alphas(Alphas<1e-3)=0;
%         sum_AlphasIsZero = sum(Alphas,1);
%         idx_AlphasIsZero = find(sum_AlphasIsZero==0);
        Alphap = mexLasso([Xp;sqrtmu * full(Ws*Alphas)], [Dp; sqrtmu * eye(size(Alphap,1))],paramp);     
        Alphap(Alphap<5e-3)=0;
       
    end
%  


    Alphas = gpuArray(full(Alphas));
    Xs = gpuArray(Xs);
    Ds = gpuArray(Ds);
    Alphap = gpuArray(full(Alphap));
    Xp = gpuArray(Xp);
    Dp = gpuArray(Dp);
    Wp = gpuArray(Wp);
    Ws = gpuArray(Ws);
 
%% Update D    
    for i=1:dictSize
       ai        =    Alphas(i,:);
       Y         =    Xs-Ds*Alphas+Ds(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Ds(:,i)    =    di;
    end

    for ii = 1:nIter_P
        for i=1:dictSize
           ai        =    Alphap(i,:);
           Y         =    Xp-Dp*Alphap+Dp(:,i)*ai;
           di        =    Y*ai';
           di        =    di./(norm(di,2) + eps);
           if sum(di)~=0
               Dp(:,i)    =    di;
           else
               Dp(:,i)    =    Dp_pre(:,i);
           end
        end     
    end    

%%     Update W   
    Wp = (1 - rho) * Wp  + rho * Alphas * Alphap' * pinv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;
    
    % % initial
    for ii = 1:nIter_W
        Ws = (1 - rho) * Ws  + rho * Alphap * Alphas' * pinv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;        
%         Ws = soft(Ws,Tw);
    end
%% TV

    Xs = (1 - rho)*Xs + rho*Ds*Alphas;      
    u = reshape(Xs', height, width, dimX);
    ux=u(:, [2:width, width],:) - u;   
    uy=u([2:height, height], :, :) - u;
    p=p + sigma*cat(3, ux, uy);
    normep=max(1, sqrt(sum(p.^2,3))); 
    p = p./repmat(normep,1,1,size(p,3));

    %----- soft shrinkage---------------
    % compute divergence in div
    div1=[p([1:height-1], :, [dimX+1:end]); zeros(1, width, dimX,'gpuArray')] - [zeros(1, width,dimX,'gpuArray'); p([1:height-1], :, [dimX+1:end])];
    div2=[p(:, [1:width-1], 1:dimX)  zeros(height, 1, dimX,'gpuArray')] - [zeros(height, 1, dimX,'gpuArray')  p(:, [1:width-1], 1:dimX)] ;
    div = sum(cat(3,div1,div2),3);
    v = soft(u+tau*div,theta);
    u=v + 1*(v-u);
    Xs = reshape(u,[],dimX)';    
    
    
    
    

    P1p = Xp - Dp * Alphap;
    P1p = P1p(:)'*P1p(:) / 2;
    P2p = lambda1 *  norm(Alphap, 1);    
    P3p = Alphas - Wp * Alphap;
    P3p = P3p(:)'*P3p(:) / 2;
    P4p = nu * norm(Wp, 'fro');
    fp = 1 / 2 * P1p + P2p + mu * (  P4p);
    
    
    P1 = Xs - Ds * Alphas;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphas, 1);    
    P3 = Alphap - Ws * Alphas;
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Ws, 'fro'); 
    fs = 1 / 2 * P1 + P2 + mu * (  P4);
    f = fp + fs;
    fm =  norm(Xp - Dp*Ws*Alphas,2);

    Alphas = gather(Alphas);
    Xs = gather(Xs);
    Ds = gather(Ds);
    Alphap = gather(Alphap);
    Xp = gather(Xp);
    Dp = gather(Dp);
    Wp = gather(Wp);
    Ws = gather(Ws);
    fp  = gather(fp);
    fs = gather(fs);
    f  = gather(f);
    fm  = gather(fm);
    P1p  = gather(P1p);
    P2p  = gather(P2p);
    P3  = gather(P3);
    P4p  = gather(P4p);
    
    time=toc;
    fprintf('%d-th: Energy: %d Time: %d\n',t,fm,time);
    rec_path = sprintf('%stempDict_SR_NL_MU%s_nIter%s.mat','./Results/',num2str(sqrtmu),num2str(t)); 
%     save tempDict_SR_NL Ds Dp Ws Wp par param f time t;
    save(rec_path,'Ds' ,'Dp','Alphas', 'Alphap','Ws', 'Wp' ,'par', 'param', 'fs','fp','P1p','P2p','P3','P4p','fm', 'time', 't')
    % fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4), E);
    
    if (abs(f_prev - f) / f < epsilon) && P3 < 3e2
        break;
%         epsilon = epsilon/10;
%         param.lambda = param.lambda/10;
    end
end




    