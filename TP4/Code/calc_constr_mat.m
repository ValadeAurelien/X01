function PP = calc_constr_mat(Nbpt, Nbtri, Coorneu, Numtri, Nbaretes, ...
                              epsilon)

PP = zeros(Nbpt, 1+Nbaretes/2);
MM = sparse(Nbpt,Nbpt); 
NbaretesNC = Nbaretes - 4;

for i=1:2
    PP(2*i-1, i) = 1;
    PP(2*i,   i) = -1;
end
PP(1, 3) = -1;
PP(4, 3) =  1;

% On fait les bords haut / bas
for i=1:NbaretesNC/4
    PP(4+i, 3+i) = 1;
    PP(4+3*NbaretesNC/4-i+1, 3+i) = -1;
end

% On fait lesNC bordesNC droite / gauche
for i=1:NbaretesNC/4
    PP(4+NbaretesNC/4+i, 3+NbaretesNC/4+i) = 1;
    PP(4+NbaretesNC-i+1, 3+NbaretesNC/4+i) = -1;
end


for l=1:Nbtri
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    Mel=matM_elem(S1, S2, S3);
    for i=1:3
        I = Numtri(l, i);
        for j=1:3
            J = Numtri(l, j);
            MM(I, J) = MM(I, J) + Mel(i, j);
        end
    end 
end

PP = [ 2*diag(MM) PP ];
PP = sparse(PP);
