function Yt = recoverImage(Ypt,Ysize,psize)


for i = 1:Ysize(1)/psize
    for j = 1:Ysize(2)/psize
        idx = (Ysize(1)-psize+1)*psize*(i-1)+psize*(j-1)+1;
        ypt = Ypt(:,idx);
        yp = reshape(ypt,psize,psize,Ysize(3)); %4 4 50
        
        Y((i-1)*psize+1:i*psize,(j-1)*psize+1:j*psize,:) = yp;

    end
end
Yt = reshape(Y,[],size(Y,3))';