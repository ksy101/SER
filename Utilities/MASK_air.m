function out = MASK_air (in,Filter)

M=size(in,3);
out = zeros(size(in));
F_m = mean(Filter,3);
mask = ones(size(F_m));
F_m(F_m<0)=0;
mask = (F_m>0);
for i = 1:M
out(:,:,i) = in(:,:,i).*mask;
end
