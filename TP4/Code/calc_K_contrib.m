function K_kl = calc_K_contrib(Atype, epsilon, ...
                               macro_x, macro_y, macro_Atype, macro_epsilon, ...
                               macro_S1, macro_S2, macro_S3, ...
                               Nbpt, Nbtri, Coorneu, Numtri, Nbaretes, PP)


%%% constructions de la matrice de rigidité autour de x_{T, k} %%%
KK = sparse(Nbpt,Nbpt);
for l=1:Nbtri
    S1=Coorneu(Numtri(l, 1), :);
    S2=Coorneu(Numtri(l, 2), :);
    S3=Coorneu(Numtri(l, 3), :);
    Kel=matK_elem_pbcell(S1, S2, S3, Atype, epsilon, ...
                         macro_x, macro_y, macro_Atype, ...
                         macro_epsilon);
    for i=1:3
        I = Numtri(l, i);
        for j=1:3
            J = Numtri(l, j);
            KK(I, J) = KK(I, J) + Kel(i, j);
        end
    end     
end

%%% Projections des éléments de la macro-base locale sur la micro-base %%%
%%% Semble bien marcher -> Décommenter plus loin pour voir les résultats  %%%
x1 = macro_S1(1); y1 = macro_S1(2);
x2 = macro_S2(1); y2 = macro_S2(2);
x3 = macro_S3(1); y3 = macro_S3(2);
D = (x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
lambda = @(xa, xb, ya, yb, XY, D, mx, my) 1/D * ...
         ( ya - yb ) * ( XY(:, 1) + repmat(mx - xb, Nbpt, 1) ) - ...
         ( xa - xb ) * ( XY(:, 2) + repmat(my - yb, Nbpt, 1) );
macro_sol_1 = lambda(x2, x3, y2, y3, Coorneu, macro_x, macro_y, D); 
macro_sol_2 = lambda(x3, x1, y3, y1, Coorneu, macro_x, macro_y, D); 
macro_sol_3 = lambda(x1, x2, y1, y2, Coorneu, macro_x, macro_y, D); 


macro_sols = [macro_sol_1, macro_sol_2, macro_sol_3];

%%% Système contraint %%%
Lconstr = [ zeros(Nbpt, 3) ;
            PP'*macro_sols ];
Kconstr = [ KK PP                ;
            PP' sparse(size(PP, 2), size(PP, 2)) ];

Uconst = Kconstr \ Lconstr;

%%% On récupère la partie qui nous intéresse  %%%
Usol = Uconst(1:Nbpt, :);

K_kl = Usol'*KK*Usol;


%%% Décommenter pour voir les résultats  %%%

% $$$ figure();
% $$$ set(gca,'DataAspectRatio',[1,1,1])
% $$$ xlabel('x');
% $$$ ylabel('y');
% $$$ plot([x1 x2 x3], ...
% $$$      [y1 y2 y3], 'o', 'Color', 'r');
% $$$ hold on;
% $$$ text(x1, y1, '1');
% $$$ text(x2, y2, '2');
% $$$ text(x3, y3, '3');
% $$$ plot([macro_x], ...
% $$$      [macro_y], 'o', 'Color', 'b');
% $$$ plot(macro_x + Coorneu(:,1), macro_y + Coorneu(:,2), '.', 'Color', 'k');
% $$$ hold off;
% $$$ affiche(macro_sol_1, Numtri, Coorneu, 'w_{1, lin} ');
% $$$ affiche(macro_sol_2, Numtri, Coorneu, 'w_{2, lin} ');
% $$$ affiche(macro_sol_3, Numtri, Coorneu, 'w_{3, lin} ');
% $$$ affiche(Usol(:,1)- macro_sol_1, Numtri, Coorneu, 'sol 1');
% $$$ affiche(Usol(:,2)- macro_sol_2, Numtri, Coorneu, 'sol 2');
% $$$ affiche(Usol(:,3)- macro_sol_3, Numtri, Coorneu, 'sol 3');
% $$$ 
% $$$ keyboard;

