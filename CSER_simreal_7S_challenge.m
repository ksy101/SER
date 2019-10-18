clear all
addpath(genpath('Data'));
addpath('lasso');
addpath('Utilities');
addpath(genpath('SPAM'));
% addpath(genpath('SPAM_linux'));
addpath('decomposition');
addpath('images_plot');
% addpath('./Data/Train/SER_simu_10P4CDL/');
% % ----------LOAD DATA--------------------
% % ------SIMU------
load('AL_3.mat');
load('AL_4.mat');
load('AH_4.mat');
load('AH_3.mat');
% load('GdIwater_5.mat');
% load('GdIwater_50.mat');
% % ------REAL------
load('Ar.mat');
% load('ro_sart_tv.mat');
% %-------10P-------
% load('YH_t.mat');
% load('YL_t.mat');
% %  with D
% load('CDL_train_simreal_1P_sizeD128_breakChain.mat');
% %----water+PMMA-----
% load('waterPMMA.mat');
% load('waterPMMA_11_50.mat');
% load('waterPMMA_tv.mat');
% 
% Mul = Mul_tv;

load('id_gd.mat');
load('id_water.mat');
%%  -------Challenge2571-----
load('2569.2571_9slice_random.mat');


YL_com_t = TRAIN(:,2:12)';
YH_com_t = TRAIN(:,13:end)';





%% ----------------CDL-------Training-----------
% rec_path = sprintf('%stempDict_SER_challenge2569.2571_NoM_Lambda1_%s_Lambda3_%s_MU%s_RHO%s_nIter%s.mat','./Results/Challenge/',num2str( 0.001 ), num2str( 0.01 ),num2str( 0.05 ),...
%     num2str(  0.1  ),num2str( 8  )); 
% load(rec_path);

rec_path = sprintf('%stempDict_SER_challenge2569.2571_9slice_NoM_NoN_Lambda1_%s_Lambda3_%s_MU%s_RHO%s_nIter%s.mat','./Results/Challenge/',...
   num2str( 0.001 ), num2str( 0.01 ),num2str( 0.05 ),num2str(  0.1  ),num2str( 132  )); 
load(rec_path);


tic
% [Alphah, Alphal, YH_D_t, YL_D_t, Dh, Dl, Wh, Wl, f] = cser(YH_com_t,YL_com_t ,1, Dp, Ds, Wp, Ws);
[Alphah, Alphal, YH_D_t, YL_D_t, Dh, Dl, Wh, Wl, f] = cser(YH_com_t,YL_com_t ,1);
time = toc;


% save CDL_train_simreal_1P_sizeD128_5to50 Alphah Alphal YH_D_t YL_D_t Dh Dl Wh Wl f
save CDL_train_simreal_2571 Alphah Alphal YH_D_t YL_D_t Dh Dl Wh Wl f time


















 