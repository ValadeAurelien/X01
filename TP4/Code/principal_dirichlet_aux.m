function [UU, MM, KK] = principal_dirichlet_aux(h, namemsh, micro_namemsh, ...
                                                macro_Atype, macro_epsilon, ...
                                                Atype, epsilon, fourpointsKquad)
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


%%% Lecture du maillage %%%
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,...
 Reftri,Nbaretes,Numaretes,Refaretes] = lecture_msh(namemsh);
[micro_Nbpt,micro_Nbtri,micro_Coorneu,micro_Refneu,micro_Numtri, ...
 micro_Reftri,micro_Nbaretes,micro_Numaretes,micro_Refaretes] = lecture_msh(micro_namemsh);

%%% Declarations %%%
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre
AllKel = zeros(Nbtri, 9);
AllMel = zeros(Nbtri, 9);

%%% Recentrage du micro-maillage %%%
micro_Coorneu = 1/epsilon * ( micro_Coorneu - ...
                              ( repmat(micro_Coorneu(1, :), micro_Nbpt, 1) + ...
                                repmat(micro_Coorneu(3, :), micro_Nbpt, ...
                                       1) ) / 2 )
;
%%% Calcul de la matrice de contrainte %%%
micro_PP = calc_constr_mat(micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                           micro_Numtri, micro_Nbaretes, epsilon);

%%% Calculs des contributions des cellules  %%%
parfor_progress(Nbtri);
for l=1:Nbtri
 
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    
    
    Kel=matK_elem(S1, S2, S3, ...
                  Atype, epsilon, macro_Atype, macro_epsilon, ...
                  micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                  micro_Numtri, micro_Nbaretes, micro_PP, fourpointsKquad);
    Mel=matM_elem(S1, S2, S3);
    
    %%% Stockage temporaire  %%%
    AllKel(l, :) = [Kel(1, :) Kel(2, :) Kel(2, :)];
    AllMel(l, :) = [Mel(1, :) Mel(2, :) Mel(3, :)];
    
    %%% Affichage avancement %%%
    parfor_progress;
end 

%%% Construction de la matrice de rigidité  %%%
for l=1:Nbtri
    Kel = [ AllKel(l, 1:3) ;
            AllKel(l, 4:6) ;
            AllKel(l, 7:9) ];
    Mel = [ AllMel(l, 1:3) ;
            AllMel(l, 4:6) ;
            AllMel(l, 7:9) ];
    for i=1:3
        I = Numtri(l, i);
        for j=1:3
            J = Numtri(l, j);
            MM(I, J) = MM(I, J) + Mel(i, j);
            KK(I, J) = KK(I, J) + Kel(i, j);
        end
    end     
end

% Calcul du second membre L et AA
% -------------------------------
FF = f(Coorneu(:,1), Coorneu(:,2), macro_Atype);
LL = MM*FF;
AA = KK;

% norm(KK, 'fro')

%%% Matrice de changement de base %%%
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

