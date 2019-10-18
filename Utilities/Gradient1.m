function [Y] = Gradient1(X)
[X_gx, X_gy] = gradient(X);
Y = zeros([size(X),2]);
Y(:,:,1) = X_gx;
Y(:,:,2) = X_gy;
