function [UU, MM, KK] = principal_dirichlet_aux(h, namemsh, namemsh_pbcell, ...
                                                Atype, epsilon, eta)
% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemsh);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

if epsilon<=0
    Ah = calc_A_homo(namemsh_pbcell, Atype, eta);
else 
    Ah = eye(2); % on s'en fout...
end
% $$$ disp(sprintf('Ah11 %.10f', Ah(1, 1)));
% $$$ disp(Ah);
% boucle sur les triangles
% ------------------------
for l=1:Nbtri
   % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    % calcul des matrices elementaires du triangle l 
    
    Kel=matK_elem(S1, S2, S3, Ah, Atype, epsilon);
    Mel=matM_elem(S1, S2, S3);
    % On fait l'assemblage de la matrice globale et du second membre
    % A COMPLETER
    for i=1:3
        I = Numtri(l, i);
        for j=1:3
            J = Numtri(l, j);
            MM(I, J) = MM(I, J) + Mel(i, j);
            KK(I, J) = KK(I, J) + Kel(i, j);
        end
    end     
end % for l

% Calcul du second membre L
% -------------------------
% A COMPLETER
% utiliser la routine f.m
FF = f(Coorneu(:,1), Coorneu(:,2), Atype);
LL = MM*FF;
AA = KK;

PP = sparse(horzcat(zeros(Nbpt-Nbaretes, Nbaretes), eye(Nbpt-Nbaretes, Nbpt-Nbaretes)));
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UU0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

