function [YH_hat_f_t,YH_hat_nf_t,X_hat_f,X_hat_nf,num,mask_k,mask_tv_total]=interL2H_kmeans_gf(YL,AL,AH,T,lambda1,lambda2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ----------INPUT-------------
% YL - reconstructed sCT images
% AL - decomposition matrix for bins
% AH - decomposition matrix from ser
% T - transform matrix
% lambda1 and lambda2 for fine decomposition: L1 and denoising
% %-------------OUTPUT------------
% YH_hat_f_t || X_hat_f - ROI-wised method - removing basis materials and
%                                            denoising
% YH_hat_nf_t ||  X_hat_nf - without  basis materials optimization and
%                            denosing
% num - population percentage
% mask_k - result of kmeans
% mask_tv_total - kmeans with tv
% -------------------------------
%% Initialization
% lambda1=1e-4;

YL_tv = zeros(size(YL));
for ii = 1:size(YL,3)
    YL_tv(:,:,ii) = TVL1denoise(YL(:,:,ii), 5e1, 100);
end
% YL=YL_tv;

YLt=reshape(YL,[],size(YL,3))';
% YL_m_t = mean(YLt,1);
YL_m = mean(YL,3);
YL_m=TVL1denoise(YL_m, 5e-1, 100);
 X_hat_mf = zeros(size(YL,1),size(YL,2),size(YL,3));
 X_hat_nf = zeros(size(YL,1),size(YL,2),size(YL,3),size(YL,3));
 X_hat_nf_t= zeros(size(AL,2),size(YLt,2));
 X_hat_f_t= zeros(size(AL,2),size(YLt,2));
 % segmentation via k-means into 2 sub-areas: mask is the matrix which the
 % only the selected area is '1', otherwise '0'
 mask_kt=kmeans(YLt',2,'MaxIter',300);
 mask_k = reshape(mask_kt,size(YL,1),size(YL,2));

 mask_tv_total=zeros(size(mask_k,1),size(mask_k,2),2);
 X_hat_totall = zeros(size(YL,1),size(YL,2),size(AH,2));

 %% main LOOP
for k = 1:2
% ROI filtered by TV
     mask_k1=zeros(size(mask_k));
     mask_tv=zeros(size(mask_k));
     mask_k1(mask_k==k)=1;   
     mask_tv1=TVL1denoise(mask_k1, 5e-2, 100);     
     mask_tv(mask_tv1>=0.5)=1;
     mask_tv_total(:,:,k)=mask_tv;
     num_musk=sum(mask_tv(:)==1);
     
% corse decomposition in dif ROI
     YLk = YL.*mask_tv;
     YLk_t=reshape(YLk,[],size(YL,3))';
     % l1 method
     [X_hat_l] = reg_L1(AL,YLk_t,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
% guided filter
%      [X_hat_l] = reg_gf_L1L2(AL,YLk_t,T,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
     X_hat=reshape(X_hat_l',size(YL,1),size(YL,2),size(AH,2)); 
     X_hat_totall = X_hat + X_hat_totall;
     
% filter decomposed materials by threshold: Relative population thresholding
     for i=1:5  
         Xf=X_hat(:,:,i);
    %      Xf=X_hat_l(i,:);
    %      num(i,k)=sum(Xf(:)~=0)/sum(musk_k1(:)~=0);
         num(i,k)=sum(Xf(:)>=1e-3)/num_musk;
         if num(i,k)>=0.6
            Xf_nz=Xf(Xf~=0);
            Xf_m=mean(Xf_nz);
            % all the area share the same value
            X_hat_mf(:,:,i,k)=mask_tv*Xf_m+X_hat_mf(:,:,i);
            % keep the original value
            X_hat_nf(:,:,i,k)=Xf;
            else
            X_hat_mf(:,:,i,k)=mask_tv*0+X_hat_mf(:,:,i);
            X_hat_nf(:,:,i,k)=Xf*0;
         end
     end

     X_hat_nf_k = squeeze(X_hat_nf(:,:,:,k));  % materials : exists and water, PMMA
     X_hat_nf_k_t = reshape(X_hat_nf_k,[],size(AL,2))';
     % 'coarse' decomosition results
     X_hat_nf_t =X_hat_nf_k_t+X_hat_nf_t;
     X_hat_f_k_t = zeros(size(X_hat_nf_k_t));
    %  AL_b=AL(:,3:end);

    % delete background
    %  Xwp = zeros(2,size(X_hat_nf_k_t,2)); 
    %  if num(1,k) > 0.8
    %  Xwp(1,:)=1;
    %  else
    %  Xwp(2,:)=1.18;
    %  end
    %  y_hat_nf_k_nb_t=YLk_t - AL(:,1:2)*Xwp;
    %  y_hat_nf_k_nb_t(y_hat_nf_k_nb_t<0)=0;

 % without removing background materials (water+PMMA) for paper
     y_hat_nf_k_b_t = YLk_t;
 % remove non-exist materials
     num_AL = find(squeeze(num(1:end,k))>0.6);
      % normilized for compariation of flier image and decomposed materials (without considering water and PMMA, dueto their bigger magnitude and )
     YL_mk = YL_m.*mask_tv;
     norm_ym = sqrt(mean(mean(YL_mk.^2)));
     norm_xm = sqrt(mean(mean(X_hat_l(3:end,:).^2)));
     norm_rho = norm_ym/norm_xm;
     YL_mk = YL_mk/norm_rho;
     YL_m_t = reshape(YL_mk,1,[]);
    %  tic
    AL_bf=AL(:,num_AL(1:end));
% FINE decomposition    
    tic
    if size(AL_bf,2)~=0
     [X_hat_gf_k_t] = reg_gf_L1L2(AL_bf,y_hat_nf_k_b_t,T,YL_m_t,num_AL,'POSITIVITY','yes','lambda1', lambda1,'lambda2',lambda2,'AL_ITERS',300, 'TOL', 1e-6); 
%       X_hat_gf_k_t  = reg_TV2(AL_bf,y_hat_nf_k_b_t,T,'POSITIVITY','yes','lambda', 4e-2,'rho', 1,'AL_ITERS',300, 'TOL', 1e-6); 
      
    %  [X_hat_l_k_t] = reg_L1(AL_f,y_hat_nf_k_b_t,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 

    X_hat_f_k_t(num_AL(1:end),:) =  X_hat_gf_k_t;
    % X_hat_f_k_t(1:2,:)=Xwp;
    end
     toc

     X_hat_f_t = X_hat_f_k_t+X_hat_f_t;
end

%  X_hat_nf_t = reshape(X_hat_nf,[],size(AL,2))';
%  X_hat_f_t = X_hat_nf_t;
%  y_hat_nf_k_t=AL*X_hat_nf_t;
%  num_1 = num;
%  num_1(num_1<=0.6)=0;
%  num_AL = find(sum(num_1,2)>0);
%  AL_f = AL(:,num_AL);
%  tic
%  [X_hat_gf_t] = reg_gf_L1(AL_f,y_hat_nf_k_t,T,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
%  toc
%  X_hat_f_t(num_AL,:) =  X_hat_gf_t;
 
     X_hat_f = reshape(X_hat_f_t',size(YL,1),size(YL,2),size(AH,2));
     YH_hat_f_t=AH(:,1:5)*X_hat_f_t;
     YH_hat_nf_t=AH(:,1:5)*X_hat_nf_t;
     X_hat_nf = sum(X_hat_nf,4);
%      
%      figure;
%      subplot(1,2,1);imagesc(X_hat_f(:,:,3));colorbar
% %      subplot(1,3,2);imagesc(X_hat_nf(:,:,3)-X_hat_f(:,:,3));colorbar
%      subplot(1,2,2);plot(X_hat_totall(50,:,3));hold on;plot(X_hat_f(50,:,3));
     
%      subplot(1,3,3);imagesc(mean(YL,3));colorbar
     
     
     
     
     