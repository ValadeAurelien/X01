function Ah = calc_A_homo(namemsh, Atype, eta)

[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemsh);

KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
for l=1:Nbtri
   % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    % calcul des matrices elementaires du triangle l 
    Kel=matK_elem_pbcell(S1, S2, S3, Atype);
    Mel=matM_elem(S1, S2, S3);
    % On fait l'assemmblage de la matrice globale et du second membre
    % A COMPLETER
    for i=1:3
        I = Numtri(l, i);
        for j=1:3
            J = Numtri(l, j);
            MM(I, J) = MM(I, J) + Mel(i, j);
            KK(I, J) = KK(I, J) + Kel(i, j);
        end
    end     
end

NbaretesNS = Nbaretes-4;
PPmod = zeros(Nbpt, 1+NbaretesNS/2);
% On fait les 4 angles sur la même colonne
for i=1:4 
    PPmod(i, 1) = 1/4.;
end
% On fait les bords haut / bas
for i=1:(NbaretesNS)/4
    PPmod(4+i, 1+i) = 1/2.;
    PPmod(4+3*NbaretesNS/4-i+1, 1+i) = 1/2.;
end
% On fait les bordes droite / gauche
for i=1:(NbaretesNS)/4
    PPmod(4+NbaretesNS/4+i, 1+NbaretesNS/4+i) = 1/2.;
    PPmod(4+NbaretesNS-i+1, 1+NbaretesNS/4+i) = 1/2.;
end
% On concatene avec les noeuds du centre 
PPcentre = vertcat(zeros(Nbaretes, Nbpt-Nbaretes), ...
                  eye(Nbpt-Nbaretes));
PP = [PPmod PPcentre];
PP = sparse(PP);

LL_1 = - KK * Coorneu(:, 1);%-(Coorneu(:, 1)' * KK)';
LL_2 = - KK * Coorneu(:, 2);%-(Coorneu(:, 2)' * KK)';

AA = KK+eta*MM;

AAp = PP'*AA*PP;
LLp_1 = PP'*LL_1;
LLp_2 = PP'*LL_2;

Wp_1 = AAp\LLp_1;
Wp_2 = AAp\LLp_2;

W_1 = PP*Wp_1;
W_2 = PP*Wp_2;
Ws = [W_1 W_2];
% $$$ plot(Coorneu(:, 1), W_1, '.');
affiche(W_1, Numtri, Coorneu, 'test');
affiche(W_2, Numtri, Coorneu, 'test 2');
% $$$ disp([max(abs(W_1)), max(abs(W_2))]);

Ah = (Coorneu + Ws)' * KK * (Coorneu + Ws);
% $$$ Ah = zeros(2);
% $$$ for i=1:2
% $$$     for j=1:2
% $$$         Ah(j, i) = ( Coorneu(:,i) + Ws(:,i) )' * KK * ...
% $$$             ( Coorneu(:,j) + Ws(:,j) );
% $$$     end
% $$$ end
