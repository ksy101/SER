
load('A_Jan2017.mat');
id = zeros(size(idx));
id(idx==5)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
muw1 = mean(img_rec1_t(id(:)==1,:),1)/1;
id = zeros(size(idx));
id(idx==6)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
muw2 = mean(img_rec1_t(id(:)==1,:),1)/1;
muw = (muw1+muw2)/2;

id = zeros(size(idx));
id(idx==2)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu10 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.01;

id = zeros(size(idx));
id(idx==3)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu12 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.012;

id = zeros(size(idx));
id(idx==11)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu2 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.002;

id = zeros(size(idx));
id(idx==10)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu1 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.001;

id = zeros(size(idx));
id(idx==9)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu05 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.0005;

id = zeros(size(idx));
id(idx==8)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu02 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.0002;

id = zeros(size(idx));
id(idx==7)=1;
img_rec1_t = reshape(img_rec1.*id,[],5);
mu01 = (mean(img_rec1_t(id(:)==1,:),1)-muw)/0.0001;

Mu_Gd = mean([mu12; mu10; mu2; mu1]);
Mu_Gd2 = mean([mu01; mu02; mu05]);
A_2571 = zeros(5,4);
A_2571(:,1) = muw';
A_2571(:,2) = A_Jan2017(:,2);
A_2571(:,3) = Mu_Gd2';
A_2571(:,4) = Mu_Gd';




