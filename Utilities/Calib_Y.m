

%% calibration
% % ground truth
load('X_gt_Jan2017_precise.mat');
load('A_Jan2017.mat');
load('id_gd.mat');
load('id_water.mat');
YL_t = A_Jan2017*reshape(X_gt,[],4)';
yl_water_gt = YL_t(:,squeeze(id_water(:))==1);
mean_yl_water_gt = mean(yl_water_gt,2);

% % Jan2017 1e-5
load('SART_Jan172017_Tau1e-05.mat');
img_rec1 = reshape(SART_TV,320,320,5);
img_rec1 = permute(img_rec1,[2,1,3]);
img_rec1 = img_rec1/0.044;
img_rec1(img_rec1>2)=0;
figure_single(img_rec1)
Y1_t = reshape(img_rec1,[],size(img_rec1,3))';
yl_water1 = Y1_t(:,squeeze(id_water(:))==1);
mean_yl_water1 = mean(yl_water1,2);
figure;plot(mean_yl_water_gt);
hold on; plot(mean_yl_water1,'-.')
factor_calib1 = mean_yl_water1./mean_yl_water_gt;
img_rec1 = img_rec1./reshape(repmat(squeeze(factor_calib1),1,320*320)',320,320,5);
figure_single(img_rec1)
SART_TV= reshape(img_rec1,[],5);
% save SART_Jan2017_TAU1e-5_Calib  SART_TV

% % Jan2017 5e-5
load('SART_Jan172017_Tau5e-05.mat');
img_rec2 = reshape(SART_TV,320,320,5);
img_rec2 = permute(img_rec2,[2,1,3]);
img_rec2 = img_rec2/0.044;
img_rec2(img_rec2>2)=0;
Y2_t = reshape(img_rec2,[],size(img_rec2,3))';
yl_water2 = Y2_t(:,squeeze(id_water(:))==1);
mean_yl_water2 = mean(yl_water2,2);
figure;plot(mean_yl_water_gt);
hold on; plot(mean_yl_water2)
factor_calib2 = mean_yl_water2./mean_yl_water_gt;
img_rec2 = img_rec2./reshape(repmat(squeeze(factor_calib2),1,320*320)',320,320,5);
SART_TV= reshape(img_rec2,[],5);
% save SART_Jan2017_TAU5e-5_Calib  SART_TV

% % Jan2017 1e-4
load('SART_Jan2017_TAU1e-4.mat');
img_rec3 = reshape(SART_TV,320,320,5);
img_rec3 = permute(img_rec3,[2,1,3]);
img_rec3 = img_rec3/0.044;
img_rec3(img_rec3>2)=0;
Y3_t = reshape(img_rec3,[],size(img_rec3,3))';
yl_water3 = Y3_t(:,squeeze(id_water(:))==1);
mean_yl_water3 = mean(yl_water3,2);
figure;plot(mean_yl_water_gt);
hold on; plot(mean_yl_water3)
factor_calib3 = mean_yl_water3./mean_yl_water_gt;
img_rec3 = img_rec3./reshape(repmat(squeeze(factor_calib3),1,320*320)',320,320,5);
SART_TV= reshape(img_rec3,[],5);
% save SART_Jan2017_TAU1e-4_Calib  SART_TV
% % Jan2017 1e-6
load('SART_Jan2017_TAU1e-6.mat');
img_rec4 = reshape(SART_TV,320,320,5);
img_rec4 = permute(img_rec4,[2,1,3]);
img_rec4 = img_rec4/0.044;
img_rec4(img_rec4>2)=0;
Y4_t = reshape(img_rec4,[],size(img_rec4,3))';
yl_water4 = Y4_t(:,squeeze(id_water(:))==1);
mean_yl_water4 = mean(yl_water4,2);
figure;plot(mean_yl_water_gt);
hold on; plot(mean_yl_water4)
factor_calib4 = mean_yl_water4./mean_yl_water_gt;
img_rec4 = img_rec4./reshape(repmat(squeeze(factor_calib4),1,320*320)',320,320,5);
SART_TV= reshape(img_rec4,[],5);
% save SART_Jan172017_TAU1e-6_Calib  SART_TV


figure;plot(mean_yl_water_gt,'-.');
hold on; plot(mean_yl_water3);plot(mean_yl_water2)
plot(mean_yl_water1);plot(mean_yl_water4)

figure;plot(factor_calib3);
hold on;plot(factor_calib2);plot(factor_calib1);plot(factor_calib4);

figure_single(img_rec3,img_rec2,img_rec1,img_rec4)
%% calculation of Mum transmation matrix -- M
% clear
load('SART_Jan172017_Tau1e-4_Calib.mat');
load('A_Jan2017.mat');
load('X_gt_Jan2017_precise.mat');

AL = A_Jan2017(:,[1,3:4]);
% img_gt = zeros(256,256,3);
% img_gt(:,:,1) = X_gt(33:288,33:288,1);
% img_gt(:,:,2:3) = X_gt(33:288,33:288,3:4);
img_gt = X_gt(:,:,[1,3:4]);

X_t = reshape(img_gt,[],size(img_gt,3))';
img_rec = reshape(SART_TV,320,320,5);
% img_rec = permute(img_rec,[2,1,3]);
% img_rec = img_rec/0.044;
% img_rec(img_rec>2)=0;
% Mul = img_rec(33:288,33:288,:);
Mul = img_rec;
Y_t = reshape(Mul,[],size(Mul,3))';
% YX = Y_t*X_t'*pinv(X_t*X_t');

error = 0;
% for i = 1:10
    AX = (AL*(X_t-error));
    [Mt] = reg_minusI_L1(AX',Y_t','POSITIVITY','yes','lambda', 9e-1,'AL_ITERS',300, 'TOL', 1e-8,'ALPHA',3,'MUU',1e0);
%     error1 = norm(Y_t - AL*X_t,2)
    error_rec = norm(Y_t - Mt'*AL*X_t,2);
    A = Mt'*AL;
    [X_hat_n_t] = reg_L1(A,Y_t,'POSITIVITY','yes','lambda', 1e-5,'AL_ITERS',300, 'TOL', 1e-6); % no-positive for simu
    error = X_hat_n_t - X_t;
    error_dec = norm(error,2);
    fprintf('%d-th: Error_reconstruction: %d Error_decomposition: %d\n',i,error_rec,error_dec);
% end

X_hat_l1_n=reshape(X_hat_n_t',size(Mul,1),size(Mul,2),size(AL,2));
X_hat_l1_n(:,:,1) = X_hat_l1_n(:,:,1)/max(max(X_hat_l1_n(:,:,1)));

[X_hat_t] = reg_L1(AL,Y_t,'POSITIVITY','yes','lambda', 1e-5,'AL_ITERS',300, 'TOL', 1e-6); % no-positive for simu
X_hat_l1=reshape(X_hat_t',size(Mul,1),size(Mul,2),size(AL,2));
X_hat_l1(:,:,1) = X_hat_l1(:,:,1)/max(max(X_hat_l1(:,:,1)));

figure;plot(AL(:,2))
hold on;plot(A(:,2))
figure_single(X_hat_l1,X_hat_l1_n)
M=Mt'

%% calibration---Feb2018

% % 5e-5
load('SART_Feb2018_TAU5e-5.mat');
img_rec2 = reshape(SART_TV,320,320,5);
img_rec2 = permute(img_rec2,[2,1,3]);
img_rec2 = img_rec2/0.044;
img_rec2(img_rec2>2)=0;
% Y2_t = reshape(img_rec2,[],size(img_rec2,3))';
% yl_water2 = Y2_t(:,squeeze(id_water(:))==1);
% mean_yl_water2 = mean(yl_water2,2);
% figure;plot(mean_yl_water_gt);
% hold on; plot(mean_yl_water2)
% factor_calib2 = mean_yl_water2./mean_yl_water_gt;
img_rec2 = img_rec2./reshape(repmat(squeeze(factor_calib2),1,320*320)',320,320,5);
SART_TV= reshape(img_rec2,[],5);
save SART_Feb2018_TAU5e-5_Calib  SART_TV

% % 1e-4
load('SART_Feb2018_TAU1e-4.mat');
img_rec3 = reshape(SART_TV,320,320,5);
img_rec3 = permute(img_rec3,[2,1,3]);
img_rec3 = img_rec3/0.044;
img_rec3(img_rec3>2)=0;
% Y3_t = reshape(img_rec3,[],size(img_rec3,3))';
% yl_water3 = Y3_t(:,squeeze(id_water(:))==1);
% mean_yl_water3 = mean(yl_water3,2);
% figure;plot(mean_yl_water_gt);
% hold on; plot(mean_yl_water3)
% factor_calib3 = mean_yl_water3./mean_yl_water_gt;
img_rec3 = img_rec3./reshape(repmat(squeeze(factor_calib3),1,320*320)',320,320,5);
SART_TV= reshape(img_rec3,[],5);
save SART_Feb2018_TAU1e-4_Calib  SART_TV


