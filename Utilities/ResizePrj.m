function  b_new = ResizePrj(b_binned)

b = reshape(b_binned,643,600,5);


b_midx = reshape(b(1:640,:,:),2,[],5);
b_mid = reshape(b_midx(1,:,:),[],600,5);

b_midy = reshape(permute(b_mid,[2, 1, 3]),2,[],5);
b_mid = permute(reshape(b_midy(1,:,:),[],size(b_mid,1),5),[2,1,3]);
b_new = zeros(320,320,5);
b_new(:,11:310,:) = b_mid;
