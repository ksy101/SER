function  img = patch2crop(img_patch)
dec_id      = [291 610; 134 453]; % for phantom
dim = size(img_patch,3);
imgsize = 900;

img = zeros(imgsize,imgsize,dim);

img_roi = img_patch(:,:,:,1);
img_out = img_patch(:,:,:,2:end);
img_out = permute(img_out,[1 2 4 3]);
img_out_t = reshape(img_out,[],dim);
p1 = (dec_id(1,1)-1)*imgsize;
p2 = (dec_id(1,1)-1)*imgsize+(dec_id(1,2)-dec_id(1,1)+1)*((dec_id(2,1)-1)+...
    (imgsize-dec_id(2,2)));
p3 =p2 + (imgsize-dec_id(1,2))*imgsize;

img(dec_id(1,1):dec_id(1,2),dec_id(2,1):dec_id(2,2),:) = img_roi;

img(1:dec_id(1,1)-1,:,:)     =  reshape(img_out_t(1:p1,:),dec_id(1,1)-1,imgsize,dim);

img(dec_id(1,1):dec_id(1,2),[1:dec_id(2,1)-1,dec_id(2,2)+1:end],:)  =  reshape(img_out_t(p1+1:p2,:),...
    (dec_id(1,2)-dec_id(1,1)+1),((dec_id(2,1)-1)+ (imgsize-dec_id(2,2))),dim);

img(dec_id(1,2)+1:end, 1:end,:)     =  reshape(img_out_t(p2+1:p3,:),(imgsize-dec_id(1,2)),imgsize,dim);