function ent = BinDec(Ind) 
% // convierte el vector binario ind a entero
two = 2; ent = 0;
n = size(Ind,2) ;

for i=1:n
    ent = ent + Ind(i)*two^(i-1) ;
end
return ;
