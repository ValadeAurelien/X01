function [Kel] = matK_elem_pbcell(S1, S2, S3, Atype)

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Transformation F_l du triangle de reference vers le triangle physique
F_l = @(xref,yref) [(x2-x1)*xref + (x3-x1)*yref + x1 , (y2-y1)*xref + (y3-y1)*yref + y1];
% Jacobienne de F_l
Jacobian = [ x2-x1 , x3-x1 ; y2-y1 , y3-y1];

% coordonees des points pour la quadrature sur le triangle physique
X1q = F_l(1/3,1/3);
X2q = F_l(1/5,1/5);
X3q = F_l(1/5,3/5);
X4q = F_l(3/5,1/5);
w1 = -9/32 ;
w234 = 25/96;

% quadrature du terme source
A_quad = w1*   A(X1q(1), X1q(2), Atype) + ...
         w234*(A(X2q(1), X2q(2), Atype) + ...
               A(X3q(1), X3q(2), Atype) + ...% les 3 normales sur le triangle de reference
               A(X4q(1), X4q(2), Atype));    

norm_ref = zeros(3, 2);
norm_ref = [-1,-1;1,0;0,1];

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    Kel(i,j) = (     (  A_quad*inv(Jacobian')*(norm_ref(i,:))'  )' * (inv(Jacobian')*norm_ref(j,:)') ) * abs(det(Jacobian));
  end; % j
end; % i

end% $$$ 
% $$$ function [Kel] = matK_elem_pbcell(S1, S2, S3, Atype)
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % mat_elem :
% $$$ % calcul la matrices de raideur elementaire en P1 lagrange
% $$$ %
% $$$ % SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
% $$$ %          
% $$$ % INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
% $$$ %                      (vecteurs reels 1x2)
% $$$ %
% $$$ % OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
% $$$ %
% $$$ % NOTE (1) le calcul est exacte (pas de condensation de masse)
% $$$ %      (2) calcul direct a partir des formules donnees par 
% $$$ %          les coordonnees barycentriques 
% $$$ %
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ % preliminaires, pour faciliter la lecture:
% $$$ x1 = S1(1); y1 = S1(2);
% $$$ x2 = S2(1); y2 = S2(2);
% $$$ x3 = S3(1); y3 = S3(2);
% $$$ 
% $$$ % Transformation F_l du triangle de reference vers le triangle physique
% $$$ F_l = @(xref,yref) [(x2-x1)*xref + (x3-x1)*yref + x1 , (y2-y1)*xref + (y3-y1)*yref + y1];
% $$$ % Jacobienne de F_l
% $$$ J_F = [ x2-x1 , x3-x1 ; y2-y1 , y3-y1];
% $$$  Jt_inv = inv(J_F');
% $$$  det_J = det(J_F);
% $$$  
% $$$  % Poids et points
% $$$  w1 = -9/32; w234 = 25/96;
% $$$  
% $$$  X1 = F_l(1/3,1/3);
% $$$  X2 = F_l(1/5,1/5);
% $$$  X3 = F_l(1/5,3/5);
% $$$  X4 = F_l(3/5,1/5);
% $$$  
% $$$  % Integrale de A
% $$$  Aq = w1*   A(X1(1), X1(2), Atype) + ...
% $$$       w234*(A(X2(1), X2(2), Atype) + ...
% $$$             A(X3(1), X3(2), Atype) + ...
% $$$             A(X4(1), X4(2), Atype));
% $$$  
% $$$  norm_ref = zeros(3, 2);
% $$$  norm_ref(1, :) = [-1, -1];
% $$$  norm_ref(2, :) = [1, 0];
% $$$  norm_ref(3, :) = [0, 1];
% $$$  
% $$$  % calcul de la matrice de raideur
% $$$  % -------------------------------
% $$$  Kel = zeros(3,3);
% $$$  for i=1:3
% $$$      for j=1:3
% $$$          Kel(i,j) = ( (Aq*Jt_inv*(norm_ref(i,:))')' * (Jt_inv*norm_ref(j,:)') ) * ...
% $$$              abs(det_J);
% $$$      end; % j
% $$$  end; % i
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %                                                        fin de la routine
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
