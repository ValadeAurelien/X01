function val = f(x, y, Atype, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%val = ones(size(x));
if (Atype == 1)
    val = -2*pi^2*sin(pi*x).*sin(pi*y);
elseif (Atype == 2)
    val = -3*pi^2*sin(pi*x).*sin(pi*y);
elseif (Atype == 3)
    val = 2*pi^2*(-epsilon*(sin(2*pi*x/epsilon) + 3).*sin(pi*x) + ...
                  2*cos(pi*x).*cos(2*pi*x/epsilon)).*sin(pi*y)/epsilon;
elseif (Atype == 4)
    val = pi^2*(-3*epsilon*(sin(2*pi*x/epsilon) + 2).*sin(pi*x) + ...
                4*cos(pi*x).*cos(2*pi*x/epsilon)).*sin(pi*y)/epsilon;
elseif (Atype == 5)
    val = 4*pi^2*(-epsilon*(sin(2*pi*x/epsilon) + 1)*(sin(2*pi*x/epsilon) + 4).*sin(pi*x) + ...
                  (sin(2*pi*x/epsilon) + 1).*cos(pi*x).*cos(2*pi*x/epsilon) + ...
                  (sin(2*pi*x/epsilon) + 4).*cos(pi*x).*cos(2*pi*x/epsilon)).*sin(pi*y)/epsilon;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%