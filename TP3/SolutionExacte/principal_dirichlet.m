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
close all;
clear all;
clc;

% PARAMETRES A RENTRER 
validation = 1; 
validation_curves = 1;
visualisation = 0; % 0, 1 (juste les diff), 2 (tout), 3 (avec A et F)
Atype = 3; % 1(Id), 2(i), 3(ii), 4(iii), 5(iv) 6(test)
vareps=1;
varmsh=0;

hs = [0.79432823472428150207, 0.63095734448019324943, 0.501187233627272285, 0.39810717055349725077, 0.3162277660168379332, 0.25118864315095801111, 0.19952623149688796014, 0.15848931924611134852, 0.12589254117941672104, 0.1, 0.07943282347242815021, 0.06309573444801932494, 0.0501187233627272285, 0.03981071705534972508, 0.03162277660168379332, 0.02511886431509580111, 0.01995262314968879601, 0.01584893192461113485, 0.0125892541179416721, 0.01];
namemshs = {'geoms/geomCarre1.msh', 'geoms/geomCarre2.msh', 'geoms/geomCarre3.msh', 'geoms/geomCarre4.msh', 'geoms/geomCarre5.msh', 'geoms/geomCarre6.msh', 'geoms/geomCarre7.msh', 'geoms/geomCarre8.msh', 'geoms/geomCarre9.msh', 'geoms/geomCarre10.msh', 'geoms/geomCarre11.msh', 'geoms/geomCarre12.msh', 'geoms/geomCarre13.msh', 'geoms/geomCarre14.msh', 'geoms/geomCarre15.msh', 'geoms/geomCarre16.msh', 'geoms/geomCarre17.msh', 'geoms/geomCarre18.msh', 'geoms/geomCarre19.msh', 'geoms/geomCarre20.msh'};

if (vareps)
    nbeps=20;
    epsilons=logspace(0, 2.5, nbeps);
    rap_L2=zeros(nbeps);
    rap_H1=zeros(nbeps);

    msh=10;
    for eps=1:nbeps
        %affichemaillage(namemshs{msh}, hs(msh));
        [diff_L2, norm_ex_L2, diff_H1, norm_ex_H1] = ...
            principal_dirichlet_aux(hs(msh), namemshs{msh}, visualisation, ...
                                    validation, Atype, epsilons(eps));
        disp([msh, hs(msh), eps, diff_L2/norm_ex_L2, diff_H1/norm_ex_H1]);
        rap_L2(eps) = diff_L2/norm_ex_L2;
        rap_H1(eps) = diff_H1/norm_ex_H1;
    end

    msh=17;
    rap2_L2=zeros(nbeps);
    rap2_H1=zeros(nbeps);
    for eps=1:nbeps
        %affichemaillage(namemshs{msh}, hs(msh));
        [diff_l2, norm_ex_l2, diff_h1, norm_ex_h1] = ...
            principal_dirichlet_aux(hs(msh), namemshs{msh}, visualisation, ...
                                    validation, Atype, epsilons(eps));
        disp([msh, hs(msh), eps, diff_l2/norm_ex_l2, diff_h1/norm_ex_h1]);
        rap2_L2(eps) = diff_l2/norm_ex_l2;
        rap2_H1(eps) = diff_h1/norm_ex_h1;
    end

    if (validation_curves)
        figure();
        loglog(epsilons, rap_L2);
        hold on;
        loglog(epsilons, rap2_L2);
        hold off;
        title('normalized error in L2 as a function of h ', ...
              'fontsize', 20);
        xlabel('1/eps','fontsize',20);
        ylabel('l2 normalized error','fontsize',20);
        grid on;
        
        figure();
        loglog(epsilons, rap_H1);
        hold on;
        loglog(epsilons, rap2_H1);
        hold off;
        title('Normalized error in H1 as a function of h ', ...
              'FontSize', 20);
        xlabel('1/eps','FontSize',20);
        ylabel('H1 Normalized error','FontSize',20);
        grid on;
    end
elseif (varmsh)
    nbmsh=13;
    mshmin=5;
    mshmax=mshmin+nbmsh-1;

    rap_L2=zeros(nbmsh);
    rap_H1=zeros(nbmsh);
    
    eps=1;
    for msh=mshmin:mshmax
        %affichemaillage(namemshs{msh}, hs(msh));
        [diff_L2, norm_ex_L2, diff_H1, norm_ex_H1] = ...
            principal_dirichlet_aux(hs(msh), namemshs{msh}, visualisation, ...
                                    validation, Atype, eps);
        disp([msh, hs(msh), eps, diff_L2/norm_ex_L2, diff_H1/norm_ex_H1]);
        rap_L2(msh-mshmin+1) = diff_L2/norm_ex_L2;
        rap_H1(msh-mshmin+1) = diff_H1/norm_ex_H1;
    end

    rap2_L2=zeros(nbmsh);
    rap2_H1=zeros(nbmsh);
    
    eps=20;
    for msh=mshmin:mshmax
        %affichemaillage(namemshs{msh}, hs(msh));
        [diff_l2, norm_ex_l2, diff_h1, norm_ex_h1] = ...
            principal_dirichlet_aux(hs(msh), namemshs{msh}, visualisation, ...
                                    validation, Atype, eps);
        disp([msh, hs(msh), eps, diff_l2/norm_ex_l2, diff_h1/norm_ex_h1]);
        rap2_L2(msh-mshmin+1) = diff_l2/norm_ex_l2;
        rap2_H1(msh-mshmin+1) = diff_h1/norm_ex_h1;
    end

    if (validation_curves)
        figure();
        loglog(1./hs(mshmin:mshmax), rap_L2);
        hold on;
        loglog(1./hs(mshmin:mshmax), rap2_L2);
        title('normalized error in L2 as a function of h ', ...
              'fontsize', 20);
        xlabel('1/msh','fontsize',20);
        ylabel('l2 normalized error','fontsize',20);
        grid on;
        hold off;        
        figure();
        loglog(1./hs(mshmin:mshmax), rap_H1);
        hold on;                   
        loglog(1./hs(mshmin:mshmax), rap2_H1);
        title('Normalized error in H1 as a function of h ', ...
              'FontSize', 20);
        xlabel('1/msh','FontSize',20);
        ylabel('H1 Normalized error','FontSize',20);
        grid on;
        hold off;
    end
end