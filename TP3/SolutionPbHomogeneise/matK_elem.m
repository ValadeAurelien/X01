function [Kel] = matK_elem(S1, S2, S3, Ah, Atype, epsilon)
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

% $$$ % les 3 normales a l'arete opposees (de la longueur de l'arete)
% $$$ norm = zeros(3, 2);
% $$$ norm(1, :) = [y2-y3, x3-x2];
% $$$ norm(2, :) = [y3-y1, x1-x3];
% $$$ norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
    disp('l aire d un triangle est nulle!!!');
    error('l aire d un triangle est nulle!!!'); 
end;

% Transformation F_l du triangle de reference vers le triangle physique
F_l = @(xref,yref) [(x2-x1)*xref + (x3-x1)*yref + x1 , (y2-y1)*xref + (y3-y1)*yref + y1];
% Jacobienne de F_l
Jacobian = [ x2-x1 , x3-x1 ; y2-y1 , y3-y1];

% coordonees des points pour la quadrature sur le triangle physique
X1q = F_l(1/3,1/3);
X2q = F_l(1/5,1/5);
X3q = F_l(1/5,3/5);
X4q = F_l(3/5,1/5);
% Poids et points
w1 = -9/32; w234 = 25/96;

% Integrale de A
if (epsilon>0)
    A_quad = w1*   A(X1q(1), X1q(2), Atype, epsilon) + ...
             w234*(A(X2q(1), X2q(2), Atype, epsilon) + ...
                   A(X3q(1), X3q(2), Atype, epsilon) + ...
                   A(X4q(1), X4q(2), Atype, epsilon));
else
    A_quad = Ah/2;
end
norm_ref = zeros(3, 2);
norm_ref = [-1,-1;1,0;0,1];

% $$$ Aval = A((x1+x2+x3)/3, (y1+y2+y3)/3, Atype, epsilon)/2;
% $$$ disp(sprintf('%f %f %f %f\n%f %f %f %f\n\n', ...
% $$$              Aval(1), Aval(2), Aval(3), Aval(4), ...
% $$$              Aq(1), Aq(2), Aq(3), Aq(4)) );            

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    Kel(i,j) = (     (  A_quad*inv(Jacobian')*(norm_ref(i,:))'  )' * (inv(Jacobian')*norm_ref(j,:)') ) * abs(det(Jacobian));
  end; % j
end; % i

% i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
