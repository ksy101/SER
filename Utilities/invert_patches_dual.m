function [HP] = invert_patches_dual(im1, patch_size, id)

% [nrow, ncol] = size(im1);


% [~,~,ndim] = size(reshape(im1(:,1),patch_size,patch_size,[]));


% HP = zeros(size(id,1)*patch_size, size(id,1)*patch_size,ndim);
for i = 1: size(id,1)
    for j = 1: size(id,2)
        im_patch = reshape(im1(:,id(i,j)),patch_size,patch_size,[]);

        HP((i-1)*patch_size+1:i*patch_size,(j-1)*patch_size+1:j*patch_size,:) = im_patch;
    end
end
        


