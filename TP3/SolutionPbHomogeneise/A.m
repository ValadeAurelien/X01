function val = A(x, y, Atype, epsilon)
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
if epsilon>=0
    if (Atype == 1)
        val = [1 0; 0 1];
    elseif (Atype == 2)
        val = [1 0; 0 3];
    elseif (Atype == 3)
    val = [2+sin(2*pi*x*epsilon) 0 ; 0 4];
    elseif (Atype == 4)
        val = [2+sin(2*pi*x*epsilon) 0 ; 0 4+sin(2*pi*x*epsilon)];
    elseif (Atype == 5)
        val = (sin(2*pi*x*epsilon) + 2)*(sin(2*pi*y*epsilon) + 4)*eye(2);
    elseif (Atype == 6)
        val =  sin(10*x*epsilon)*eye(2);
    end
else
    if (Atype == 1)
        val = [1 0; 0 1];
    elseif (Atype == 2)
        val = [1 0; 0 3];
    elseif (Atype == 3)
    val = [sqrt(3) 0 ; 0 4];
    elseif (Atype == 4)
        val = [sqrt(3) 0 ; 0 4];
    elseif (Atype == 5)
        val = [ 4*sqrt(3) 0 ; 0 2*sqrt(15) ]
    elseif (Atype == 6)
        val =  sin(10*x*epsilon)*eye(2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
