function [Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws,f] = coupled_DL(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)
% Based on Semi-Coupled Dictionary Learning <Wang2012>


[dimX, numX]        =       size(Xp);
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
nIter_W             =       par.nIter_W;
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
% paramp.lambda2 =  lambda2*1e-2;
% nIter_W = nIter_P;

% % Dp_pre
% paramD.K = par.K;
% paramD.iter=300; 
% paramD.mode=2;
% paramD.lambda  = lambda3; % 1e-1    
% Dp_pre = mexTrainDL(Xp, paramD);
fs = 0; fp = 0;
P1p = 0;P2p = 0;P3 = 0; P4p = 0; fm = 0;
time = 0;
for t = 1 : nIter
    
    
    tic
    f_prev = f;
%% update sparsity coefficients   
        num_pix = 320*320;
%     if numX == (320*320)*NUM_P


        Alphas_sgl = zeros(par.K,num_pix,NUM_P);
        Alphap_sgl = zeros(par.K,num_pix,NUM_P);
        for num_P = 1:NUM_P
            Xs_sgl = Xs(:,(num_P-1)*num_pix+1:num_P*num_pix);
            Xp_sgl = Xp(:,(num_P-1)*num_pix+1:num_P*num_pix);        
%             Alphas_sgl(:,:,num_P) = mexLasso(Xs_sgl,[Ds; sqrtmu * Ws] ,param); % size [128 65536]
            Alphas_sgl(:,:,num_P) = mexLasso(Xs_sgl, Ds, param); % broke the chain
            Alphap_sgl(:,:,num_P) = mexLasso([Xp_sgl;sqrtmu * full(Ws*Alphas_sgl(:,:,num_P))], [Dp; sqrtmu * eye(size(Alphap_sgl(:,:,num_P),1))], paramp);
        end

        Alphas = reshape(Alphas_sgl,par.K,[]);                        % size [128 655360]
        Alphap = reshape(Alphap_sgl,par.K,[]);
 

        
        
%     else     
     
        % % initial
%         Alphas = mexLasso([Xs;sqrtmu * full(Alphap)], [Ds; sqrtmu * Ws],param);


        % % broke the chain 
        % % v1: || Alphas - Wp * Alphap ||
%         Alphas = mexLasso(Xs, Ds, param);
%         Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);


%         % % v2: || Alphap - Ws * Alphas ||
%         Alphas = mexLasso(Xs, Ds, param);
% %         Alphas(Alphas<1e-3)=0;
% 
% %         sum_AlphasIsZero = sum(Alphas,1);
% %         idx_AlphasIsZero = find(sum_AlphasIsZero==0);
% 
%         Alphap = mexLasso([Xp;sqrtmu * full(Ws*Alphas)], [Dp; sqrtmu * eye(size(Alphap,1))],paramp);     
% %         Alphap = 0.5*mexLasso(Xp, Dp, paramp) + 0.5*Ws*Alphas; 
% %         Alphap = 0.5*Alphap + 0.5*mexLasso(full(Ws*Alphas),  eye(size(Alphap,1)),paramp); 
% %         Alphap = mexLasso([Xp;sqrtmu * Xs], [Dp; sqrtmu * Ds],paramp);     
% 
% %         Alphap(Alphap<5e-3)=0;
% %         Alphap(Alphas==0)=0;
% 
%         
%         % % v3
% %         Alphas = mexLasso(Xs, Ds, param); 
% %         Xp_D = Dp*Ws*Alphas;
% %         Alphap = mexLasso([Xp_D;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);
% %         Alphas = mexLasso([Xs;sqrtmu*full(Alphap)], [Ds;sqrtmu*Ws],param); 
% %         Alphap = mexLasso([Xp;sqrtmu * full(Alphas)], [Dp; sqrtmu * Wp],param);
%         
% 
%         % % v4: || WAlpha ||_1
% %         Alphas = mexLasso([Xs;beta*(z+d)], [Ds;beta*Ws], param);
% %         Alphap = mexLasso([Xp;sqrtmu * full(Ws*Alphas)], [Dp; sqrtmu * eye(size(Alphap,1))],param);
% %         z = mexLasso(full(Ws*Alphas)-d,eye(size(z,1)),param);
% %         d = d - (Ws*Alphas-z);
%         
% 
%         % % v5: Alphap = Ws * Alphas 
% %         Alphas = mexLasso(Xs, Ds, param);
% %         Alphap = Alphas;     
%         
% 
%     end
%  


%     Dp_pre1 = mexTrainDL([Xp;sqrtmu * full(Ws*Alphas)], paramD);
%     Dp_pre = Dp_pre1(1:size(Dp,1),:);





    Alphas = gpuArray(full(Alphas));
    Xs = gpuArray(Xs);
    Ds = gpuArray(Ds);
    Alphap = gpuArray(full(Alphap));
    Xp = gpuArray(Xp);
    Dp = gpuArray(Dp);
%     Wp = gpuArray(Wp);
    Ws = gpuArray(Ws);
    

    
    
dictSize = par.K;


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
%            if sum(di)~=0
               Dp(:,i)    =    di;
%            else
%                Dp(:,i)    =    Dp_pre(:,i);
%            end
        end     
    end    


    
    
%%     Update W
%     Ws = Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;
%     Wp = Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;    
    
%     Wp = (1 - rho) * Wp  + rho * Alphas * Alphap' * pinv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;
    
    % % initial
    for ii = 1:nIter_W
        Ws = (1 - rho) * Ws  + rho * Alphap * Alphas' * pinv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;        
%         Ws = soft(Ws,Tw);
    end

    % % for || Ws * Alpha ||_1
%     Ws = ((sqrtmu*Alphap  +beta*(z+d))* Alphas')* pinv((sqrtmu+beta)*Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;

    % % for || Ws ||_1
%     WsT =  mexLasso(Alphap',Alphas',param_w); 
%     Ws = (1 - rho) * Ws + rho * WsT';
    

    % % Ws = Wp^-1
%     Ws = (1 - rho) * Ws  + rho * (Alphas * Alphap' + pinv(Wp)) * pinv(Alphap * Alphap' + (par.nu + 1) * eye(size(Alphap, 1)));

    % % use CNN
    



    P1p = Xp - Dp * Alphap;
    P1p = P1p(:)'*P1p(:) / 2;
    P2p = lambda1 *  norm(Alphap, 1);    
%         size(full(Alphas))
%         size(full(Wp ))
%     size(full(Wp * Alphap))
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
%     Wp = gather(Wp);
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

        rec_path = sprintf('%stempDict_SER_challenge2569.2571_9slice_random_NoM_Lambda1_%s_Lambda3_%s_MU%s_RHO%s_nIter%s.mat','./Results/Challenge/',num2str(lambda1),num2str(lambda3),num2str(sqrtmu),num2str(rho),num2str(t)); 
        save(rec_path,'Ds' ,'Dp','Alphas', 'Alphap','Ws', 'Wp' ,'par', 'param', 'fs','fp','P1p','P2p','P3','P4p','fm', 'time', 't')
    
    if (abs(f_prev - f) / f < epsilon) 
        break;
%         epsilon = epsilon/10;
%         param.lambda = param.lambda/10;
    end
end




    