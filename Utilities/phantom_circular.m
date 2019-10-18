function Mul=phantom_circular(Mul,R)

Mul1=Mul;
clear Mul
     for i=1:800
         for j= 1:800
            if (i-400)^2+(j-400)^2 > R*R
                continue;
            end
            Mul(i,j, :,:)=Mul1(i,j,:,:);
         end
     end