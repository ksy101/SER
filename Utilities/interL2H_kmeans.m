function [YH_hat_filter_t,YH_hat_nofilter_t,X_hat_filter,X_hat_nofilter,num,musk_k,musk_tv_total]=interL2H_kmeans(YL,AL,AH,lambda1)
%% filter by groundtruth
% filter_one=Xgt;
YLt=reshape(YL,[],size(YL,3))';
% XLt=reshape(XL,[],size(XL,3))';
% if ~exist('lambda2')    
%     
%     [X_hat_l] = reg_L1(AL,YLt,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
%     
%     X_hat=reshape(X_hat_l',size(YL,1),size(YL,2),size(AH,2));
% else
%     
%     [X_hat_l1] = reg_L1(AL,YLt,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
%     [X_hat_l2] = reg_L1(AL,YLt,'POSITIVITY','yes','lambda', lambda2,'AL_ITERS',300, 'TOL', 1e-6);
%     X_hat_l=[X_hat_l2(1:2,:);X_hat_l1(3:end,:)];
%     
%     X_hat=reshape(X_hat_l',size(YL,1),size(YL,2),size(AH,2));
%     
% end
% YH_hat_t=AH(:,1:5)*X_hat_l;
% 
% 
% X_hat_filter=X_hat;
% for i=3:5  
% Xf=X_hat(:,:,i).*filter_one;
% num(i)=sum(Xf(:)~=0);
% if sum(Xf(:)~=0)>=2e3
%  Xf_nz=Xf(Xf~=0);
% Xf_m=mean(Xf_nz);
% X_hat_filter(:,:,i)=filter_one*Xf_m;
% else
%     X_hat_filter(:,:,i)=filter_one*0;
% end
%  Xfct=reshape(X_hat_filter,[],size(AL,2))';
%  YH_hat_filter_t=AH(:,1:5)*Xfct;
% end

 %% filter by kmeans
 
 musk_kt=kmeans(YLt',2,'MaxIter',5000);
 musk_k = reshape(musk_kt,size(YL,1),size(YL,2));
 X_hat_filter = zeros(size(YL,1),size(YL,2),size(YL,3));
 X_hat_nofilter = zeros(size(YL,1),size(YL,2),size(YL,3));
 musk_tv_total=zeros(size(musk_k,1),size(musk_k,2),2);
 
 
for k = 1:2
     musk_k1=zeros(size(musk_k));
     musk_tv=zeros(size(musk_k));
     musk_k1(musk_k==k)=1;   
     musk_tv1=TVL1denoise(musk_k1, 5e-2, 100);     
     musk_tv(musk_tv1>=0.5)=1;
     musk_tv_total(:,:,k)=musk_tv;
     num_musk=sum(musk_tv(:)==1);
     YLk = YL.*musk_tv;
     YLk_t=reshape(YLk,[],size(YL,3))';
     [X_hat_l] = reg_L1(AL,YLk_t,'POSITIVITY','yes','lambda', lambda1,'AL_ITERS',300, 'TOL', 1e-6); 
     X_hat=reshape(X_hat_l',size(YL,1),size(YL,2),size(AH,2));

 for i=1:5  
     Xf=X_hat(:,:,i);
%      num(i,k)=sum(Xf(:)~=0)/sum(musk_k1(:)~=0);
     num(i,k)=sum(Xf(:)>=1e-3)/num_musk;
     if num(i,k)>=0.6
        Xf_nz=Xf(Xf~=0);
        Xf_m=mean(Xf_nz);
        X_hat_filter(:,:,i)=musk_tv*Xf_m+X_hat_filter(:,:,i);
        X_hat_nofilter(:,:,i)=Xf+X_hat_nofilter(:,:,i);
        else
        X_hat_filter(:,:,i)=musk_tv*0+X_hat_filter(:,:,i);
        X_hat_nofilter(:,:,i)=Xf*0+X_hat_nofilter(:,:,i);
     end
 end
 
end
%      X_hat_filter(:,:,1:2)=XL(:,:,1:2);
     Xfct=reshape(X_hat_filter,[],size(AL,2))';
     YH_hat_filter_t=AH(:,1:5)*Xfct;
     Xnfct=reshape(X_hat_nofilter,[],size(AL,2))';
     YH_hat_nofilter_t=AH(:,1:5)*Xnfct;
        
     
%      index=sum(YLt,1)==0;
%      index1=sum(YLt,1)~=0;
%      index_t=reshape(index,size(YL,1),size(YL,2));
%      Xfct1=Xfct;
%      Xfct1(1,index)=mean(XLt(1,index));
%      Xfct1(1,index1)=mean(XLt(1,index1));
%      Xfct1(2,index)=mean(XLt(2,index));
%      Xfct1(2,index1)=mean(XLt(2,index1));     
%      Xfct1_w=Xfct1(1,:);
%      Xfct1_w(Xfct1_w>1)=1;    
%      Xfct1(1,:)=Xfct1_w;
%      Xfct1_p=Xfct1(2,:);
%      Xfct1_p(Xfct1_p>1.18)=1.18;    
%      Xfct1(2,:)=Xfct1_p;     
%      Xfct1(1:2,:)=XLt(1:2,:);
%      X_hat_filter1=reshape(Xfct1',size(YL,1),size(YL,2),size(AH,2));
%      YH_hat_filter_t1=AH(:,1:5)*Xfct1;
     
%      YH_hat_t=AH(:,1:5)*XLt;    
%      YH_hat_filter_t1=YH_hat_filter_t;
%      YH_hat_filter_t1(:,sum(YLt,1)==0)=YH_hat_t(:,sum(YLt,1)==0);
     


% YH_hat_filter_t=YH_hat_filter_t+AH(:,1:5)*X_hat_l;

%for the ferformance of decomposition
% X_hat_l=reshape(X_hat_l',721,780,5);
% figure;imagesc(X_hat_l(:,:,1));
% figure;imagesc(X_hat_l(:,:,2));
% figure;imagesc(X_hat_l(:,:,3));
% figure;imagesc(X_hat_l(:,:,4));
% figure;imagesc(X_hat_l(:,:,5));





% 
% figure;
% plot([20:1:79],squeeze(xh(50,50,:)))
% hold on;plot([25:10:75],squeeze(xl(50,50,:)),'*')
% 
% xlabel('energy')
% ylabel('\mu')
% legend('60 bins', '6 bins', 'SER')

% 
% figure;
% plot([20:1:79],squeeze(xh(50,50,:)))
% hold on;plot([20:1:79],squeeze(Yll(50,50,:)))

% for ii=1:60
% imll(:,:,ii)=xl(:,:,ceil(ii/10));
% end
% 
% 
% figure;
% plot([20:1:79],squeeze(xh(50,50,:)))
% hold on;plot([25:10:75],squeeze(xl(50,50,:)),'*')
% hold on;plot([20:1:79],squeeze(imll(50,50,:)),'.')
% hold on;plot([20:1:79],squeeze(iml(50,50,:)),'-.')
% xlabel('energy')
% ylabel('\mu')
% hl=legend('60 bins', '6 bins', 'flat','inter decomposition')
% set(hl,'Box','off');
