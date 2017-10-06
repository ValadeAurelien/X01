function [diff_L2, norm_ex_L2, diff_H1, norm_ex_H1 ] = principal_dirichlet_aux(h, namemsh, visualisation, validation, Acst, alpha)
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

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
   % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l,1), :);
    S2=Coorneu(Numtri(l,2), :);
    S3=Coorneu(Numtri(l,3), :);
    % calcul des matrices elementaires du triangle l 
    
    Kel=matK_elem(S1, S2, S3, Acst, alpha);
    Mel=matM_elem(S1, S2, S3);
    
    % On fait l'assemmblage de la matrice globale et du second membre
    % A COMPLETER
    for i=1:3
        I = Numtri(l,i);
        for j=1:3
            J = Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i, j);
            KK(I,J) = KK(I,J) + Kel(i, j);
        end
    end     
end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1), Coorneu(:,2), Acst, alpha);
LL = MM*FF;

% Projection sur l espace V_0
% ———————————————————
% matrice de projection 
PPtmp = eye(Nbpt, Nbpt);
for a=1:Nbaretes
    PPtmp(Numaretes(a, 1), Numaretes(a, 1)) = 0;
end
c=1;
PP=zeros(Nbaretes, Nbpt);
for i=1:Nbpt    
    if PPtmp(i,i)>0
        PP(c,:)=PPtmp(i,:);
        c=c+1;
    end
end
% $$$ for a=1:Nbaretes
% $$$     V = zeros(Nbpt,1); 
% $$$     V(Numaretes(10,1),1) = 1;
% $$$     disp((PP*V)==zeros(Nbpt-Nbaretes,1));
% $$$ end
AA = MM+KK;
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UU0;
UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
% visualisation
% -------------
if (visualisation>0)
    affiche(abs(UU-UU_exact), Numtri, Coorneu, sprintf(['Dirichlet diff - ' '%s'], namemsh));
    if (visualisation>1)
            affiche(UU, Numtri, Coorneu, sprintf(['Dirichlet UU - %s'], namemsh));           
            affiche(UU_exact, Numtri, Coorneu, sprintf(['Dirichlet ' ...
                                'UU exact - %s'], namemsh));
    end
end

% validation
% ----------
if (validation)
    % Calcul de l erreur L2
    diff = UU_exact - UU;

    diff_L2 = sqrt((diff)'*MM*diff);
    norm_ex_L2 = sqrt((UU_exact)'*MM*UU_exact);
    diff_H1 = sqrt(abs((diff)'*KK*diff));
    norm_ex_H1 = sqrt(abs((UU_exact)'*(KK)*UU_exact));
    % attention de bien changer le terme source (dans FF)
else
    
    diff_L2 = 1;
    norm_ex_L2 = 1;
    diff_H1 = 1;
    norm_ex_H1 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

