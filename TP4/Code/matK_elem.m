function Kel = matK_elem(S1, S2, S3, ...
                         Atype, epsilon, macro_Atype, macro_epsilon, ...
                         micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                         micro_Numtri, micro_Nbaretes, micro_PP, fourpointsKquad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
    disp('l aire d un triangle est nulle!!!');
    error('l aire d un triangle est nulle!!!'); 
end;

%%% Quadrature à quatres points  %%%
if (fourpointsKquad)
    % Transformation F_l du triangle de reference vers le triangle physique
    F_l = @(x1, x2, x3, y1, y2, y3, xref, yref) ...
          [(x2-x1)*xref + (x3-x1)*yref + x1 , ...
           (y2-y1)*xref + (y3-y1)*yref + y1];
    % Jacobienne de F_l
    Jacobian = [ x2-x1 , x3-x1 ; ...
                 y2-y1 , y3-y1];

    % coordonees des points pour la quadrature sur le triangle physique
    X1q = F_l(x1, x2, x3, y1, y2, y3, 1/3,1/3);
    X2q = F_l(x1, x2, x3, y1, y2, y3, 1/5,1/5);
    X3q = F_l(x1, x2, x3, y1, y2, y3, 1/5,3/5);
    X4q = F_l(x1, x2, x3, y1, y2, y3, 3/5,1/5);
    % Poids et points
    w1 = -9/32; w234 = 25/96;

    % Integrale de K_q
    Kel_q = w1*   calc_K_contrib(Atype, epsilon, ...
                                 X1q(1), X1q(2), macro_Atype, macro_epsilon, ...
                                 S1, S2, S3, ...
                                 micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                                 micro_Numtri, micro_Nbaretes, micro_PP) + ...
            w234*(calc_K_contrib(Atype, epsilon, ...
                                 X2q(1), X2q(2), macro_Atype, macro_epsilon, ...
                                 S1, S2, S3, ...
                                 micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                                 micro_Numtri, micro_Nbaretes, micro_PP) + ...
                  calc_K_contrib(Atype, epsilon, ...
                                 X3q(1), X3q(2), macro_Atype, macro_epsilon, ...
                                 S1, S2, S3, ...
                                 micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                                 micro_Numtri, micro_Nbaretes, micro_PP) + ...
                  calc_K_contrib(Atype, epsilon, ...
                                 X4q(1), X4q(2), macro_Atype, macro_epsilon, ...
                                 S1, S2, S3, ...
                                 micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                                 micro_Numtri, micro_Nbaretes, micro_PP) ...
                  );

    %%% Mauvaise formule....  %%%
    Kel = Kel_q * 2 * epsilon * abs(det(Jacobian));

%%% Quadrature à un point  %%%
else
    %%% Quadrature au barycentre  %%%
    Xq = [ (x1 + x2 + x3)  (y1 + y2 + y3)] / 3;
   
    %%% epsilon ou epsilon^2 ? %%%
    %%% De toute manière l'erreur n'est pas contante  %%%
    Kel = epsilon^2 * D/2 * ...
          calc_K_contrib(Atype, epsilon, ...
                         Xq(1), Xq(2), macro_Atype, macro_epsilon, ...
                         S1, S2, S3, ...
                         micro_Nbpt, micro_Nbtri, micro_Coorneu, ...
                         micro_Numtri, micro_Nbaretes, ...
                         micro_PP);
end

%%% Décommenter pour comparer avec la solution du TP1  %%%
% $$$ n = [ y2-y3 x3-x2;
% $$$       y3-y1 x1-x3;
% $$$       y1-y2 x1-x2];
% $$$ Kel2 = zeros(3,3);
% $$$ for i=1:3
% $$$     for j=1:3
% $$$         Kel2(i, j) = n(i, 1)*n(j,1) + ...
% $$$             n(i, 2)*n(j, 2);
% $$$     end
% $$$ end
% $$$ disp(sprintf('%.7f %.7f %.7f', norm(Kel), norm(Kel2), norm(Kel)/norm(Kel2)));
