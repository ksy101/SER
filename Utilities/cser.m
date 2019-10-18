function [Alphah, Alphal, YH_D_t, YL_D_t, Dh, Dl, Wh, Wl, f] = cser(varargin)
%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------

YH_t = varargin{1};
YL_t = varargin{2};   
Train = varargin{3}; % 1 for train with W ; 0 for test; 2 for train without W
% par.id_patch =  varargin{4};
    
    [B,N]=size(YL_t);
    [E,N]=size(YH_t);

    par.lambda1         =       1e-3;            % Alphas:  1e-3 for size256; 8e-3 for 300
    par.lambda2         =       1e-5;
    par.lambda3         =       1e-2;            % Alphap:   1.2e0 for TYPE3
    par.sqrtmu          =       5e-2;            % effect of AlphaL on ALphaH, sqrtmu should be improved with lambda3; small enough to make sure the reconstructin of H image: 8e-2 for TYPE1; 1e-1 for type3
    par.mu              =       1e-2;            % 

    par.rho             =       1e-1;              % the effect of ancent W--5e-2; 
    par.nu              =       1e-2;              % efficient of identity of W--1e-1
    par.K               =       1024;            % Dic size   
    par.nIter_W         =       1;              % iterations of traing W       
    par.Tw              =       1e-3;
    
%     par.lambda1         =       4e-5;           %1e-3 for size256; 8e-3 for 300
%     par.lambda2         =       1e-5;
%     par.lambda3         =       4e-1;            %paramp
%     par.Tw              =       5e-3;
    par.nIter           =       1000;
    par.epsilon         =       1e-7;    
    par.K               =       1024;            % Dic size
    par.NUM_P           =       N/(320*320);              % num of phantoms for  training
    par.nIter_P         =       1;              % iterations of traing Dic_p
    par.tau             =       1e-4;
    par.theta           =       5e-3;
    

%%     

    param.K = par.K;
    param.iter=300; 
    param.mode=2;
    param.lambda  = par.lambda1; % 1e-1
% if nargin>=4
%     D = varargin{4};
% else
%     D = mexTrainDL([YH_t;YL_t], param);
% end    
% %     Dh = D(1:E,:);
% %     Dl = D(E+1:end,:);

%     Dh = mexTrainDL(YH_t, param);
%     Dl = mexTrainDL(YL_t, param);
if nargin>=4
    Dh = varargin{4};
    Dl = varargin{5};
else
    Dh = mexTrainDL(YH_t, param);
    Dl = mexTrainDL(YL_t, param);    
end
    
    

     
    param.lambda2       = par.lambda2;
    param.approx=0;
     

    
    
    Alphah = mexLasso(YH_t, Dh,param);
    
    Alphal = mexLasso(YL_t, Dl,param);
%     Alphal = varargin{5};
    
if Train==1

    if nargin>=6
        Wh = varargin{6};
        Wl = varargin{7};
    else
        Wl = eye(size(Dh, 2));
        Wh = eye(size(Dl, 2));
    end 
%     size(Wh)
%     size(Wl)
    [Alphah, Alphal, ~, ~, Dh, Dl, Wh, Wl, f] = coupled_DL(Alphah, Alphal, YH_t, YL_t, Dh, Dl, Wh, Wl, par);
    YH_D_t = Dh*Alphah;
    YL_D_t = Dl*Alphal;
elseif Train==2


    [Alphah, Alphal, ~, ~, Dh, Dl, f] = coupled_DL_noW(Alphah, Alphal, YH_t, YL_t, Dh, Dl, par);
    YH_D_t = Dh*Alphah;
    YL_D_t = Dl*Alphal;
    
elseif Train==3

    Wl = eye(size(Dh, 2));
    Wh = eye(size(Dl, 2));
    [Alphah, Alphal, ~, ~, Dh, Dl, Wh, Wl, f] = coupled_DL_TV(Alphah, Alphal, YH_t, YL_t, Dh, Dl, Wh, Wl, par);
    YH_D_t = Dh*Alphah;
    YL_D_t = Dl*Alphal;
else 
    if nargin>=5
        Wh = varargin{5};
        Wl = varargin{6};
    end    
    Alphal = mexLasso([YL_t;par.sqrtmu*(YH_t)], [Dl;par.sqrtmu*Dh*Wl],param);
    Alphah = mexLasso([YH_t;par.sqrtmu * full(Alphal)], [Dh; par.sqrtmu * Wh],param);
    
    [Alphah, Alphal] = coupled_DL_test(Alphah, Alphal, YH_t, YL_t, Dh, Dl, Wh, Wl, par);
end

