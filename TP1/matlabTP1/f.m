function val = f(x, y, Acst, alpha)
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
if (Acst)
    val = (1+5*pi^2)*cos(pi*x).*cos(2*pi*y); % A = 1
else
    val = pi^2*alpha*sin(pi*x).*sin(pi*alpha*y).*cos(2*pi*y).*cos(pi*alpha*x) - ...
          2*pi^2*alpha*sin(2*pi*y).*sin(pi*alpha*x).*cos(pi*x).*cos(pi*alpha*y) + ...
          5*pi^2*(sin(pi*alpha*x).*sin(pi*alpha*y) + 2).*cos(pi*x).*cos(2*pi*y) + ...
          cos(pi*x).*cos(2*pi*y);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%