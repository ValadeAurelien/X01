function [UU, MM, KK] = principal_dirichlet_aux(h, namemsh, namemsh_pbcell, ...
                                                visualisation, validation, ...
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

Ah = calc_A_homo(namemsh, Atype, eta)
% boucle sur les triangles
% ------------------------
for l=1:Nbtri
   % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    % calcul des matrices elementaires du triangle l 
    
    Kel=matK_elem(S1, S2, S3, Ah);
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
end % for l

% Calcul du second membre L
% -------------------------
% A COMPLETER
% utiliser la routine f.m
FF = f(Coorneu(:,1), Coorneu(:,2), Atype, epsilon);
LL = MM*FF;
AA = KK;

PP = horzcat(zeros(Nbpt-Nbaretes, Nbaretes), eye(Nbpt-Nbaretes, Nbpt-Nbaretes));
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UU0;
% $$$ UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
% $$$ 
% $$$ AA_exact = zeros(Nbpt, 1);
% $$$ for i=1:Nbpt
% $$$     Amat = A(Coorneu(i,1), Coorneu(i,2), Atype, epsilon);
% $$$     AA_exact(i) = Amat(1);
% $$$ end

% visualisation
% -------------
if (visualisation>0)
% $$$     affiche(abs(UU-UU_exact), Numtri, Coorneu, sprintf(['Dirichlet ' ...
% $$$                         'diff - %s, eta : %f'], namemsh, eta));
    if (visualisation>1)
        affiche(UU, Numtri, Coorneu, sprintf(['Dirichlet UU - %s, ' ...
                            'eta : %f'], namemsh, eta));           
% $$$         affiche(UU_exact, Numtri, Coorneu, sprintf(['Dirichlet ' ...
% $$$                             'UU exact - %s, eta : %f'], namemsh, eta));
% $$$         if (visualisation>2)
% $$$             affiche(AA_exact, Numtri, Coorneu, sprintf(['Dirichlet ' ...
% $$$                                 'AA - %s, eta : %f'], namemsh, eta));           
% $$$             affiche(FF, Numtri, Coorneu, sprintf(['Dirichlet FF - ' ...
% $$$                                 '%s, eta : %f'], namemsh, eta));
% $$$         end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

