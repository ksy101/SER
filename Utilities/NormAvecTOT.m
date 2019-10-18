function Mul_norm = NormAvecTOT(Mul,img_tot)

Mul_norm = Mul./repmat((img_tot+eps),1,1,size(Mul,3));
Mul_norm(isnan(Mul_norm))=0;
Mul_norm(isinf(Mul_norm))=0;
Mul_norm(Mul_norm>2)=0;
Mul_max = max(Mul_norm(:));
Mul_norm = Mul_norm/Mul_max;