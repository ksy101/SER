function [X_hat_f,X_hat_ngf,num,mask,mask_k,mask_gmm,YL_gmm]=decompose_kmeans_gf(YL,AL,T,options)
%  guarded image filter with population thresholding
%  INPUT:     option:
%             lambda1  L1-norm sparse
%             lambda2  guarded image filter
%             pt       population thresholding
%             dt       density thresholding

lambda1 = options.lambda1;
lambda2 = options.lambda2;
pt = options.pt;
dt = options.dt;
num_k = options.num_k;
gam = options.gam;
lambda = options.lambda;
B = size(YL,3);

%% filter by groundtruth
YLt=reshape(YL,[],B)';
YL_m = mean(YL,3);

% if any(YL_m(:))
%     num_k = 4;
% else 
%     num_k = 3;
% end

 %% filter by kmeans

 % segmentation via k-means into 2 sub-areas
 mask_kt=kmeans(YLt(1:end,:)',num_k,'MaxIter',300);
 mask_k = reshape(mask_kt,size(YL,1),size(YL,2));
%   figure;imagesc(mask_k)
 %% GMM
 YL_gmm_t = zeros(size(YLt));
 YL_gmm = zeros(size(YL));
 mask_gmm = zeros(size(mask_k));
 if sum(sum(YL(:)))
    for i =1:B
        YLi(:,:,i) = TVL1denoise(YL(:,:,i), 9e-1, 100,1);     
        Yi= reshape(YLi(:,:,i),[],1);
        gm = fitgmdist(Yi,num_k, 'Regularize',1e-6);
        [idx_c(:,i),nlogL(i)] = cluster(gm,Yi);
    end
    [nlogL_max, nlogL_idx] = min(nlogL);
    mask_gmm = reshape(idx_c(:,nlogL_idx),size(mask_k));
    for ii = 1:num_k
        idx_ii = find(idx_c(:,nlogL_idx)==ii);
        for iii = 1:size(idx_ii)
            YL_gmm_t(:,idx_ii(iii)) = median(YLt(:,idx_ii),2);
        end
    end
    YL_gmm=reshape(YL_gmm_t',size(YL));
 end
 %% kernelkmeans
 YL_gmm_t =reshape(YL_gmm,[],B)';
 [mask_t] = kernelkmeans(YLt,YL_gmm_t,num_k,mask_kt,gam,lambda);
 mask = reshape(mask_t,size(YL,1),size(YL,2));
 
 
 
 
  %%
%   close all
 X_hat_mf = zeros(size(YL,1),size(YL,2),size(YL,3));
 X_hat_ngf = zeros(size(YL,1),size(YL,2),size(AL,2),num_k);
%  X_hat_nf_t= zeros(size(AL,2),size(YLt,2));
 X_hat_f_t= zeros(size(AL,2),size(YLt,2));
 
%  mask_tv_total=zeros(size(mask_k,1),size(mask_k,2),2);
%  X_hat_totall = zeros(size(YL,1),size(YL,2),size(AL,2));
 
for k = 1:num_k
    % segmented areas filtered by TV
     mask_k1=zeros(size(mask));
%      mask_tv=zeros(size(mask_k));
     mask_k1(mask==k)=1;   
%      mask_tv1=TVL1denoise(mask_k1, 5e-2, 100);     
%      mask_tv(mask_tv1>=0.5)=1;
%      mask_tv_total(:,:,k)=mask_tv;
     num_musk=sum(mask_k1(:)==1);
     
    % corse decomposition in dif areas
     YLk = YL.*mask_k1;
     YLk_t=reshape(YLk,[],size(YL,3))';
    % l1 method
     [X_hat_t] = reg_L1(AL,YLk_t,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
    % guided filter
%      [X_hat_l] = reg_gf_L1L2(AL,YLk_t,T,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
     X_hat=reshape(X_hat_t',size(YL,1),size(YL,2),size(AL,2)); 
%      X_hat_totall = X_hat + X_hat_totall;
     
    % filter decomposed materials by threshold: Relative population thresholding
     for i=1:size(AL,2) 
         Xf=X_hat(:,:,i);
    %      Xf=X_hat_l(i,:);
    %      num(i,k)=sum(Xf(:)~=0)/sum(musk_k1(:)~=0);
         num(i,k)=sum(Xf(:)>=dt)/num_musk;
         if num(i,k)>=pt
            Xf_nz=Xf(Xf~=0);
            Xf_m=mean(Xf_nz);
            % all the area share the same value
            X_hat_mf(:,:,i)=mask_k1*Xf_m+X_hat_mf(:,:,i);
            % keep the original value
            X_hat_ngf(:,:,i,k)=Xf;
            else
            X_hat_mf(:,:,i)=mask_k1*0+X_hat_mf(:,:,i);
            X_hat_ngf(:,:,i,k)=Xf*0;
         end
     end

%      X_hat_nf_k = squeeze(X_hat_ngf(:,:,:,k));  % materials : exists and water, PMMA
%      X_hat_nf_k_t = reshape(X_hat_nf_k,[],size(AL,2))';
%      X_hat_nf_t =X_hat_nf_k_t+X_hat_nf_t;
    
    %  AL_b=AL(:,3:end);

 
     % without delete background
     y_hat_nf_k_b_t = YLk_t;
     % delete non-exist materials
     num_AL = find(squeeze(num(1:end,k))>pt);
     % normilized for compariation of flier image and decomposed materials (without considering water and PMMA, dueto their bigger magnitude and )
     YL_mk = YL_m.*mask_k1;
     norm_ym = sqrt(mean(mean(YL_mk.^2)));
     if  size(num_AL,1)>1 & ismember(1,num_AL)
         norm_xm = sqrt(mean(mean(X_hat_t(num_AL(2:end),:).^2)));
     else
         norm_xm = sqrt(mean(mean(X_hat_t(num_AL,:).^2)));
     end
     norm_rho = norm_ym/norm_xm;
     YL_mk = YL_mk/norm_rho;
     YL_m_t = reshape(YL_mk,1,[]);

    AL_bf=AL(:,num_AL);
    X_hat_f_k_t = zeros(size(AL,2),size(YLt,2));
    
    if size(AL_bf,2)~=0
          [X_hat_gf_k_t] = reg_gf_L1L2(AL_bf,y_hat_nf_k_b_t,T,YL_m_t,num_AL,'POSITIVITY','yes','lambda1', lambda1,'lambda2',lambda2,'AL_ITERS',300, 'TOL', 1e-6); 
          
          
          X_hat_f_k_t(num_AL,:) =  X_hat_gf_k_t;   
%           xx=reshape(X_hat_f_k_t',size(YL,1),size(YL,2),size(AL,2));
%           figure;imagesc(xx(:,:,1));figure;imagesc(xx(:,:,2));
          X_hat_f_t = X_hat_f_k_t+X_hat_f_t;
%           X_hat_f = reshape(X_hat_f_t',size(YL,1),size(YL,2),size(AL,2));
%           figure;imagesc(X_hat_f(:,:,1));figure;imagesc(X_hat_f(:,:,2));
     end
end
    X_hat_ngf = sum(X_hat_ngf,4);
%     X_hat_f_t(1,:) = X_hat_ngf(1,:);
    X_hat_f = reshape(X_hat_f_t',size(YL,1),size(YL,2),size(AL,2));
     

     
     
     
     
     