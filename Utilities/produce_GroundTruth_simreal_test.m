

I = X_gt(:,:,1);
I1 = zeros(size(I));
I2 = zeros(size(I));
I3 = zeros(size(I));
I4 = zeros(size(I));
I5 = zeros(size(I));
I6 = zeros(size(I));
I7 = zeros(size(I));
I8 = zeros(size(I));
I9 = zeros(size(I));
I10 = zeros(size(I));

I1c = zeros(size(I));
I2c = zeros(size(I));
I3c = zeros(size(I));
I4c = zeros(size(I));
I5c = zeros(size(I));
I6c = zeros(size(I));
I7c = zeros(size(I));
I8c = zeros(size(I));
I9c = zeros(size(I));
I10c = zeros(size(I));

I1c(Idx_e{1},Idx_e{2}) = I(Idx_e{1},Idx_e{2});
I2c(Idx_e{3},Idx_e{4}) = I(Idx_e{3},Idx_e{4});
I3c(Idx_e{5},Idx_e{6}) = I(Idx_e{5},Idx_e{6});
I4c(Idx_e{7},Idx_e{8}) = I(Idx_e{7},Idx_e{8});
I5c(Idx_e{9},Idx_e{10}) = I(Idx_e{9},Idx_e{10});
I6c(Idx_e{11},Idx_e{12}) = I(Idx_e{11},Idx_e{12});
I7c(Idx_e{13},Idx_e{14}) = I(Idx_e{13},Idx_e{14});
I8c(Idx_e{15},Idx_e{16}) = I(Idx_e{15},Idx_e{16});
I9c(Idx_e{17},Idx_e{18}) = I(Idx_e{17},Idx_e{18});
I10c(Idx_e{19},Idx_e{20}) = I(Idx_e{19},Idx_e{20});

I1(I1c==1)=0.021;
I2(I2c==1)=0.022;
I3(I3c==1)=0.023;
I4(I4c==1)=0.024;
I5(I5c==1)=0.025;
I6(I6c==1)=0.026;
I7(I7c==1)=0.027;
I8(I8c==1)=0.028;
I9(I9c==1)=0.029;
I10(I10c==1)=0.030;

X_gt_test=I1+I2+I3+I4+I5+I6+I7+I8+I9+I10;