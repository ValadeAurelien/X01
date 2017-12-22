function val = A(x, y, Atype, epsilon, ...
                 macro_x, macro_y, macro_Atype, macro_epsilon)
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

%%% Partie micro  %%%
if epsilon>=0
    if (Atype == 1)
        val = [1 0; 0 1];
    elseif (Atype == 2)
        val = [1 0; 0 2];
    elseif (Atype == 3)
    val = [2+sin(2*pi*x*epsilon) 0 ; 0 4];
    elseif (Atype == 4)
        val = [2+sin(2*pi*x*epsilon) 0 ; 0 4+sin(2*pi*x*epsilon)];
    elseif (Atype == 5)
        val = (sin(2*pi*x*epsilon) + 2)*(sin(2*pi*y*epsilon) + 4)*eye(2);
    elseif (Atype == 6)
        val =  sin(10*x*epsilon)*eye(2);
    end
else %%% On veut récupérer la matrice homogène %%%
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

%%% Partie macro  %%%
%%% Inutilisée jusqu'à présent  %%%
if (macro_Atype == 1)
    val = val * [1 0; 0 1];
elseif (macro_Atype == 2)
    val = val * [1 0; 0 3];
elseif (macro_Atype == 3)
    val = val * [2+sin(2*pi*macro_x*macro_epsilon) 0 ; 0 4];
elseif (macro_Atype == 4)
    val = val * [2+sin(2*pi*macro_x*macro_epsilon) 0 ; 0 4+sin(2*pi*macro_x*macro_epsilon)];
elseif (macro_Atype == 5)
    val = val * (sin(2*pi*macro_x*macro_epsilon) + 2)*(sin(2*pi*macro_y*macro_epsilon) + 4)*eye(2);
elseif (macro_Atype == 6)
    val = val * sin(10*macro_x*macro_epsilon)*eye(2);
end

