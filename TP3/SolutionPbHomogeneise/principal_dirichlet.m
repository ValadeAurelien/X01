% ====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%c
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================

% PARAMETRES A RENTRER 
validation = 1; 
validation_curves = 1;
visualisation = 2; % 0, 1 (juste les diff), 2 (tout), 3 (avec A et F)
Atype = 4; % 1(Id), 2(i), 3(ii), 4(iii), 5(iv) 6(test)
eta = 0.01;
epsilon = -1;
msh_pbcell=13;
msh=14;


hs = [0.79432823472428150207, 0.63095734448019324943, 0.501187233627272285, 0.39810717055349725077, 0.3162277660168379332, 0.25118864315095801111, 0.19952623149688796014, 0.15848931924611134852, 0.12589254117941672104, 0.1, 0.07943282347242815021, 0.06309573444801932494, 0.0501187233627272285, 0.03981071705534972508, 0.03162277660168379332, 0.02511886431509580111, 0.01995262314968879601, 0.01584893192461113485, 0.0125892541179416721, 0.01];
namemshs = {'geoms/geomCarre1.msh', 'geoms/geomCarre2.msh', 'geoms/geomCarre3.msh', 'geoms/geomCarre4.msh', 'geoms/geomCarre5.msh', 'geoms/geomCarre6.msh', 'geoms/geomCarre7.msh', 'geoms/geomCarre8.msh', 'geoms/geomCarre9.msh', 'geoms/geomCarre10.msh', 'geoms/geomCarre11.msh', 'geoms/geomCarre12.msh', 'geoms/geomCarre13.msh', 'geoms/geomCarre14.msh', 'geoms/geomCarre15.msh', 'geoms/geomCarre16.msh', 'geoms/geomCarre17.msh', 'geoms/geomCarre18.msh', 'geoms/geomCarre19.msh', 'geoms/geomCarre20.msh'};
namemshs_pbcell = {'geoms_per/geomCarre1.msh', 'geoms_per/geomCarre2.msh', 'geoms_per/geomCarre3.msh', 'geoms_per/geomCarre4.msh', 'geoms_per/geomCarre5.msh', 'geoms_per/geomCarre6.msh', 'geoms_per/geomCarre7.msh', 'geoms_per/geomCarre8.msh', 'geoms_per/geomCarre9.msh', 'geoms_per/geomCarre10.msh', 'geoms_per/geomCarre11.msh', 'geoms_per/geomCarre12.msh', 'geoms_per/geomCarre13.msh', 'geoms_per/geomCarre14.msh', 'geoms_per/geomCarre15.msh', 'geoms_per/geomCarre16.msh', 'geoms_per/geomCarre17.msh', 'geoms_per/geomCarre18.msh', 'geoms_per/geomCarre19.msh', 'geoms_per/geomCarre20.msh'};

close all;
disp(sprintf('Solution homogène', epsilon));
[UUh, MM, KKh] = ...
        principal_dirichlet_aux(hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
                                Atype, -1, eta); 

[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemshs{msh});
affiche(UUh, Numtri,Coorneu, ...
        sprintf('Solution homogène', epsilon));
oldKK = KKh;
for epsilon=logspace(-4,3,4)
    disp(sprintf('epsilon = %g', epsilon));
    [UU, MM, KK] = principal_dirichlet_aux( ...
        hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
        Atype, epsilon, eta); 
    disp(sum(sum(KK-oldKK)));
    oldKK=KK;
    if (visualisation>0)
        affiche(abs(UUh-UU), Numtri,Coorneu, ...
                sprintf('DIFF eps : %g', epsilon));
        if(visualisation>1)
             affiche(UU, Numtri,Coorneu, ...
                     sprintf('eps : %g', epsilon));
        end
    end
end


% $$$ rap_L2=zeros(20);
% $$$ rap_H1=zeros(20);
% $$$     for msh=14:14
% $$$         %affichemaillage(namemshs{msh}, hs(msh));
% $$$         [diff_L2, norm_ex_L2, diff_H1, norm_ex_H1, UU, UU_exact, MM, KK] = ...
% $$$             principal_dirichlet_aux(hs(msh), namemshs{msh}, visualisation, ...
% $$$                                     validation, Atype, epsilon, eta);
% $$$         disp([msh, hs(msh), diff_L2/norm_ex_L2, diff_H1/norm_ex_H1]);
% $$$         rap_L2 = [rap_L2 diff_L2/norm_ex_L2];
% $$$         rap_H1 = [rap_H1 diff_H1/norm_ex_H1];
% $$$     end
% $$$ 
% $$$     if (validation_curves)
% $$$         figure(msh+1);
% $$$         loglog(-hs, rap_L2);
% $$$         title('Normalized error in L2 as a function of h (!!! log-log !!!)', ...
% $$$               'FontSize', 25);
% $$$         xlabel('h','FontSize',20);
% $$$         ylabel('L2 Normalized error','FontSize',20);
% $$$         figure(msh+2);
% $$$         loglog(-hs, rap_H1);
% $$$         title('Normalized error in H1 (semi) as a function of h (!!! log-log !!!)', ...
% $$$               'FontSize', 25);
% $$$         xlabel('h','FontSize',20);
% $$$         ylabel('H1 Normalized error','FontSize',20);
% $$$     end