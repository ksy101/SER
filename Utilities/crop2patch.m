function img_patch = crop2patch(img)
dim = size(img,3);
dec_id      = [291 610; 134 453]; % for phantom
xysize = 320;

img_roi = img(dec_id(1,1):dec_id(1,2),dec_id(2,1):dec_id(2,2),:); 
img_out_t = cat(1,reshape(img(1:dec_id(1,1)-1,:,:),[],dim),reshape(img(dec_id(1,1):dec_id(1,2),[1:dec_id(2,1)-1,dec_id(2,2)+1:end],:),[],dim),...
    reshape(img(dec_id(1,2)+1:end, :,:),[],dim));
img_out_t_0 = [img_out_t; zeros(xysize*xysize-mod(size(img_out_t,1),xysize*xysize),dim)];
img_out = reshape(img_out_t_0,xysize,xysize,[],dim);
img_out = permute(img_out,[1 2 4 3]);
img_patch = zeros(xysize,xysize,dim,size(img_out,4)+1);
img_patch(:,:,:,1) = img_roi;
img_patch(:,:,:,2:size(img_out,4)+1) = img_out;