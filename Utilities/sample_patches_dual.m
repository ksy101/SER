function [HP,LP,id_patch,id_xy] = sample_patches_dual(varargin)
im1         =       varargin{1};
im2         =       varargin{2};
patch_size  =       varargin{3}; 
[nrow, ncol, ~] = size(im1);
if nargin>=4
    id_xy          =       varargin{4}; 
    x = id_xy(1,:);
    y = id_xy(2,:);
else
    x = randperm(nrow-patch_size+1) ;
    y = randperm(ncol-patch_size+1) ;
end


id_xy = [x;y];
[X,Y] = meshgrid(x,y);

xrow = X(:);
ycol = Y(:);
patch_num = length(xrow);


for i = 1: nrow/patch_size
    for j = 1: nrow/patch_size
        id_x(i) = find(X(1,:)==(patch_size*(i-1)+1));
        id_y(j) = find(Y(:,1)==(patch_size*(j-1)+1));
        id_patch(i,j) = (id_x(i)-1) * (ncol-patch_size+1) + id_y(j);
    end
end
        



for ii = 1:patch_num
    row = xrow(ii);
    col = ycol(ii);
    
    Hpatch = im1(row:row+patch_size-1,col:col+patch_size-1,:);
    Lpatch = im2(row:row+patch_size-1,col:col+patch_size-1,:);
 
    HP(:,ii) = Hpatch(:);
    LP(:,ii) = Lpatch(:);

end

