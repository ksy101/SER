clear
load('AH_4.mat');
load('A_Jan2017.mat');
load('SART_Jan2017.mat');
% load('SART_Jan172017_Tau5e-05.mat');
load('Mod_A_Jan2017.mat');
img_rec = reshape(SART_TV,320,320,5);
img_rec = permute(img_rec,[2,1,3]);
img_rec = img_rec/0.044;
img_rec(img_rec>2)=0;

%%
[a,b,c] = size(img_rec);
idx_tot = zeros(a,b);
idx_tot0 = zeros(a,b);
idx_tot1 = zeros(a,b);
idx = zeros(a,b);
x = [61.5  92.5    162.5   229.5   255.5   224.5   154.5   87.5    119.5   159.5   196.5   156.5];
y = [163.5 230.5   256.5   225.5   156.5   88.5    62.5    93.5    160.5   198.5   158.5   120.5];
for k = 1:length(x)
    idx = zeros(a,b);
    for i = 1:a
        for j = 1:b
            if sqrt((i-x(k))^2+ (j-y(k))^2)<=16

                idx(i,j) = k;
                
            end
        end
    end
    idx_tot = idx_tot + idx;
end

x0 = [61.5     255.5         159.5      156.5];
y0 = [163.5    156.5         198.5     120.5];
for k = 1:length(x0)
    idx = zeros(a,b);
    for i = 1:a
        for j = 1:b
            if sqrt((i-x0(k))^2+ (j-y0(k))^2)<=16

                idx(i,j) = k;
                
            end
        end
    end
    idx_tot0 = idx_tot0 + idx;
end
idx_tot0(idx_tot0==2)=5;
idx_tot0(idx_tot0==3)=10;
idx_tot0(idx_tot0==4)=12;

x1 = [  92.5    162.5   229.5     224.5   154.5   87.5    119.5     196.5];
y1 = [ 230.5   256.5   225.5      88.5    62.5    93.5    160.5     158.5];
for k = 1:length(x1)
    idx = zeros(a,b);
    for i = 1:a
        for j = 1:b
            if sqrt((i-x1(k))^2+ (j-y1(k))^2)<=12

                idx(i,j) = k;
                
            end
        end
    end
    idx_tot1 = idx_tot1 + idx;
end
idx_tot0(idx_tot1==1)=2;
idx_tot0(idx_tot1==2)=3;
idx_tot0(idx_tot1==3)=4;
idx_tot0(idx_tot1==4)=6;
idx_tot0(idx_tot1==5)=7;
idx_tot0(idx_tot1==6)=8;
idx_tot0(idx_tot1==7)=9;
idx_tot0(idx_tot1==8)=11;

figure;imagesc(idx_tot0)

X_gt = zeros(a,b,4);
X_w = zeros(a,b);
X_P = zeros(a,b);
X_I = zeros(a,b);
X_Gd = zeros(a,b);

X_I(idx_tot0==3) = 0.014;  
X_I(idx_tot0==9) = 0.007;  
X_I(idx_tot0==6) = 0.0035;  
X_I(idx_tot0==8) = 0.001;  
X_Gd(idx_tot0==2) = 0.016;  
X_Gd(idx_tot0==11) = 0.008;  
X_Gd(idx_tot0==4) = 0.004;  
X_Gd(idx_tot0==7) = 0.002; 
X_w(idx_tot~=0) = 1;
X_w(idx_tot==5) = 0;
X_w(idx_tot==10) = 0;


for i = 1:a
    for j = 1:b
        if sqrt((i-158)^2+ (j-160)^2)<=127

            X_P(i,j) = 1.18;

        end
    end
end
X_P(idx_tot~=0) = 0;

figure;imagesc(X_P)

X_gt(:,:,1) = X_w;
X_gt(:,:,2) = X_P;
X_gt(:,:,3) = X_I;
X_gt(:,:,4) = X_Gd;



X_gt_t = reshape(X_gt,[],4)';
Y_gt_t = AH_4*X_gt_t;
Muh = reshape(Y_gt_t',a,b,size(AH_4,1));


