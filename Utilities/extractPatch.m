function Ypt = extractPatch(Yt,Ysize,psize)
Y = reshape(Yt',Ysize);
for i = 1:Ysize(1)-psize+1
    for j = 1:Ysize(2)-psize+1
        yp = Y(i:i+psize-1,j:j+psize-1,:);
        ypt = reshape(yp,[],1);
        if i == 1 && j == 1
            Ypt = ypt;
        else
            Ypt=cat(2,Ypt,ypt);
        end
    end
end
        