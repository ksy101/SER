function out = MASK_density (in)

out = in;
for i = 1:2
ini = in(:,:,i);
ini(ini<0.1)=0;
out(:,:,i) = ini;
end
in4 = in(:,:,4);
in4(in4<0)=0;
out(:,:,4) = in4;
