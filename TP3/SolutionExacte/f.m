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
    val = 2*pi^2*sin(pi*x).*sin(pi*y);
elseif (Atype == 2)
    val = 3*pi^2*sin(pi*x).*sin(pi*y);
elseif (Atype == 3)
    val = pi^2*(-2*epsilon*cos(pi*x).*cos(2*pi*epsilon*x) + ...
                sin(pi*x).*sin(2*pi*epsilon*x) + 6*sin(pi*x)).*sin(pi*y);
elseif (Atype == 4)
    val = 2*pi^2*(-epsilon*cos(pi*x).*cos(2*pi*epsilon*x) + ...
                  sin(pi*x).*sin(2*pi*epsilon*x) + 3*sin(pi*x)).*sin(pi*y);
elseif (Atype == 5)
    val = 2*pi^2*...
          (-epsilon*(sin(2*pi*epsilon*x) + 2).*sin(pi*x).*cos(pi*y).*cos(2*pi*epsilon*y) - ...
           epsilon*(sin(2*pi*epsilon*y) + 4).*sin(pi*y).*cos(pi*x).*cos(2*pi*epsilon*x) + ...
           (sin(2*pi*epsilon*x) + 2)*(sin(2*pi*epsilon*y) + 4).*sin(pi*x).*sin(pi*y));
elseif (Atype == 6) 
    val = 2*pi*(pi*sin(10*x).*sin(pi*x) - 5*cos(10*x).*cos(pi*x)).*sin(pi*y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%