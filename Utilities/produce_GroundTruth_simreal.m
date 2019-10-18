
I = XH(:,:,3);
I4 = zeros(size(I));
I3 = zeros(size(I));
I2 = zeros(size(I));
I1 = zeros(size(I));
I4c = zeros(size(I));
I3c = zeros(size(I));
I2c = zeros(size(I));
I1c = zeros(size(I));
I4c(110:150,200:240) = I(110:150,200:240);
I3c(70:110,110:150) = I(70:110,110:150);
I2c(170:210,45:85) = I(170:210,45:85);
I1c(45:85,45:85) = I(45:85,45:85);

I4(I4c>0.012)=0.014;
I3(I3c>0.006)=0.007;
I2(I2c>0.003)=0.0035;
I1(I1c>0.0008)=0.001;


figure;imagesc(I4)
figure;imagesc(I3)
figure;imagesc(I2)
figure;imagesc(I1)

X_gt(:,:,3)=I1+I2+I3+I4;
%%
Gd = XH(:,:,4);
Gd4 = zeros(size(Gd));
Gd3 = zeros(size(Gd));
Gd2 = zeros(size(Gd));
Gd1 = zeros(size(Gd));
Gd4c = zeros(size(Gd));
Gd3c = zeros(size(Gd));
Gd2c = zeros(size(Gd));
Gd1c = zeros(size(Gd));
Gd4c(45:85,170:210) = Gd(45:85,170:210);
Gd3c(150:190,105:145) = Gd(150:190,105:145);
Gd2c(170:210,170:210) = Gd(170:210,170:210) ;
Gd1c(105:145,20:60) = Gd(105:145,20:60);

Gd4(Gd4c>0.014)=0.016;
Gd3(Gd3c>0.006)=0.008;
Gd2(Gd2c>0.003)=0.004;
Gd1(Gd1c>0.001)=0.002;


figure;imagesc(Gd4)
figure;imagesc(Gd3)
figure;imagesc(Gd2)
figure;imagesc(Gd1)

X_gt(:,:,4)=Gd1+Gd2+Gd3+Gd4;

%%
water = X_gt(:,:,3)+X_gt(:,:,4);
water(water~=0)=1;

waterc1= F_m(20:55,110:145);
waterc = zeros(size(waterc1));
waterc(waterc1<0.18)=1;
water(20:55,110:145) = waterc;

waterc2= F_m(110:145,70:110);
waterc = zeros(size(waterc2));
waterc(waterc2<0.18)=1;
water(110:145,70:110) = waterc;


figure;imagesc(water)

X_gt(:,:,1)=water;

%%

PMMAc = XH(:,:,2);
figure;imagesc(PMMAc)
PMMA = zeros(size(PMMAc));
PMMA(PMMAc>0.1)=1.18;
figure;imagesc(PMMA)
X_gt(:,:,2)=PMMA;
save X_gt_simreal X_gt





