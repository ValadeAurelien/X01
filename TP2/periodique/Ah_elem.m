function Ahel = Ah_elem(Sn1, Sn2, Sn3, Atype, W_1, W_2)

S1=Coorneu(Sn1);
S2=Coorneu(Sn2);
S3=Coorneu(Sn3);

WS1_1 = W_1(Sn1); WS1_2 = W_2(Sn1); 
WS2_1 = W_1(Sn2); WS2_2 = W_2(Sn2);
WS3_1 = W_1(Sn3); WS3_2 = W_2(Sn3);
% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% Transfo depuis le triangle unité
F_l = @(x,y) [(x3-x1)*x + (x2-x1)*y + x1 , ...
              (y3-y1)*x + (y2-y1)*y + y1];
% ... Et son jacobien
J_F = [ x3-x1 , x2-x1 ; y3-y1 , y2-y1];
Jt_inv = inv(J_F');
det_J = det(J_F);

% Poids et points
w1 = -9/32; w234 = 25/96;

X1 = F_l(1/3,1/3);
X2 = F_l(1/5,1/5);
X3 = F_l(1/5,3/5);
X4 = F_l(3/5,1/5);

% Integrale de A
Aq = w1*A(X1(1), X1(2), Atype) + ...
     w234*(A(X2(1), X2(2), Atype) + ...
           A(X3(1), X3(2), Atype) + ...
           A(X4(1), X4(2), Atype);


norm_ref = zeros(3, 2); %dérivées des vecteurs de la base sur ce triangle
norm_ref(1, :) = [-1, -1]; 
norm_ref(2, :) = [1, 0];
norm_ref(3, :) = [0, 1];

dW_1 = WS1_1*norm_ref(1, :) + ...  
       WS2_1*norm_ref(2, :) + ...
       WS3_1*norm_ref(3, :);
dW_2 = WS1_2*norm_ref(1, :) + ...
       WS2_2*norm_ref(2, :) + ...
       WS3_2*norm_ref(3, :);

Id = eye(2);
vsums = zeros(2);
vsums(1,:) = Id(1,:) + dW_1; % e_1 + dw1
vsums(2,:) = Id(2,:) + dW_2; % e_2 + dw2

% calcul de la matrice de raideur
% -------------------------------
Ahel = zeros(2);
for i=1:2
    for j=1:2
        Ahel(i,j) = ( (Aq*Jt_inv*(vsums(i,:))')' * ...
                     (Jt_inv*vsums(j,:)') ) * ...
            abs(det_J);
    end; % j
end; % i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
