function val = A(x,y, Atype, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
if (Atype == 1)
    val = [1 0; 0 1];
elseif (Atype == 2)
    val = [1 0; 0 2];
elseif (Atype == 3)
    val = [2+2*sin(2*pi*epsilon*x) 0 ; 0 4];
elseif (Atype == 4)
    val = [2+2*sin(2*pi*epsilon*x) 0 ; 0 4+sin(2*pi*epsilon*x)];
elseif (Atype == 5)
    val = (2*sin(2*pi*epsilon*x) + 2)*(sin(2*pi*epsilon*y) + 4)*eye(2);
elseif (Atype == 6)
    val =  sin(10*x)*eye(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
