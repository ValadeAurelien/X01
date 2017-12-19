function [UU, MM, KK] = principal_dirichlet_aux(h, namemsh, micro_namemsh, ...
                                                macro_Atype, macro_epsilon, ...
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
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,...
 Reftri,Nbaretes,Numaretes,Refaretes] = lecture_msh(namemsh);
[micro_Nbpt,micro_Nbtri,micro_Coorneu,micro_Refneu,micro_Numtri, ...
 micro_Reftri,micro_Nbaretes,micro_Numaretes,micro_Refaretes] = lecture_msh(micro_namemsh);
% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre

%%% calculs à ne faire qu'une fois pour le micro maillage %%%
micro_PP = calc_constr_mat(micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                           micro_Numtri, micro_Nbaretes);




% $$$ disp(sprintf('Ah11 %.10f', Ah(1, 1)));
% $$$ disp(Ah);
% boucle sur les triangles
% ------------------------
for l=1:Nbtri
    disp(sprintf('\r%i /  %i (%.1f %%)             ', l, Nbtri, floor(l/Nbtri*1000)/10));
   % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    % calcul des matrices elementaires du triangle l 
    
    Kel=matK_elem(S1, S2, S3, ...
                  macro_Atype, Atype, macro_epsilon, epsilon, ...
                  micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                  micro_Numtri, micro_Nbaretes, micro_PP);
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

% Calcul du second membre L et AA
% -------------------------------
FF = f(Coorneu(:,1), Coorneu(:,2), macro_Atype);
LL = MM*FF;
AA = KK;

%%% MATRICE DE CHANGEMENT DE BASE %%%
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

