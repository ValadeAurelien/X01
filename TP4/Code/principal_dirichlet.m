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

%%% ENSEMBLES DES MAILLAGES ET DES PAS DE TEMPS %%%
hs = [0.79432823472428150207, 0.63095734448019324943, 0.501187233627272285, 0.39810717055349725077, 0.3162277660168379332, 0.25118864315095801111, 0.19952623149688796014, 0.15848931924611134852, 0.12589254117941672104, 0.1, 0.07943282347242815021, 0.06309573444801932494, 0.0501187233627272285, 0.03981071705534972508, 0.03162277660168379332, 0.02511886431509580111, 0.01995262314968879601, 0.01584893192461113485, 0.0125892541179416721, 0.01];
namemshs = {'geoms/geomCarre1.msh', 'geoms/geomCarre2.msh', 'geoms/geomCarre3.msh', 'geoms/geomCarre4.msh', 'geoms/geomCarre5.msh', 'geoms/geomCarre6.msh', 'geoms/geomCarre7.msh', 'geoms/geomCarre8.msh', 'geoms/geomCarre9.msh', 'geoms/geomCarre10.msh', 'geoms/geomCarre11.msh', 'geoms/geomCarre12.msh', 'geoms/geomCarre13.msh', 'geoms/geomCarre14.msh', 'geoms/geomCarre15.msh', 'geoms/geomCarre16.msh', 'geoms/geomCarre17.msh', 'geoms/geomCarre18.msh', 'geoms/geomCarre19.msh', 'geoms/geomCarre20.msh'};
namemshs_pbcell = {'geoms_per/geomCarre1.msh', 'geoms_per/geomCarre2.msh', ...
                   'geoms_per/geomCarre3.msh', 'geoms_per/geomCarre4.msh', 'geoms_per/geomCarre5.msh', 'geoms_per/geomCarre6.msh', 'geoms_per/geomCarre7.msh', 'geoms_per/geomCarre8.msh', 'geoms_per/geomCarre9.msh', 'geoms_per/geomCarre10.msh', 'geoms_per/geomCarre11.msh', 'geoms_per/geomCarre12.msh', 'geoms_per/geomCarre13.msh', 'geoms_per/geomCarre14.msh', 'geoms_per/geomCarre15.msh', 'geoms_per/geomCarre16.msh', 'geoms_per/geomCarre17.msh', 'geoms_per/geomCarre18.msh', 'geoms_per/geomCarre19.msh', 'geoms_per/geomCarre20.msh'};

%%% PARAMETRES A RENTRER %%%
msh=10;
msh_pbcell=11;
validation = 1; 
validation_curves = 1;
visualisation = 1; % 0, 1 (juste les diff), 2 (tout), 3 (avec A et
                   % F)
fourpointsKquad=0;
macro_Atype = 1;  % 1(Id), 2(i), 3(ii), 4(iii), 5(iv) 6(test)
macro_epsilon = 1; % eps<0 --> homo % pas important ici. 
                   % ATTENTION le epsilon du code est l'inverse du epsilon théorique
Atype = 4;
epsilon = floor(10/hs(msh));
bac_a_sable = 1;
vareps=0; % variation de epsilon et comparaison avec le cas homo
vareta=0; % variation et convergence de Ahom en fonction de eta
varmsh=0; % variation et convergence de Ahom en fonction du pas



%%% BAC A SABLE %%%
if (bac_a_sable)
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemshs{msh});
    [UU, ~, ~] =  principal_dirichlet_aux(hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
                                          macro_Atype, macro_epsilon, ...
                                          Atype, epsilon, fourpointsKquad);
    affiche(UU, Numtri, Coorneu, 'test');
    %%% PREMIÈRE EXPÉRIENCE %%%
elseif (vareps)
    disp(sprintf('Solution homogène'));
    [UUh, MMh, KKh] = ...
        principal_dirichlet_aux(hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
                                Atype, -1, eta); 

    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemshs{msh});
    affiche(UUh, Numtri,Coorneu, ...
            sprintf('Solution homogène'));

    msh=17;
    nbptseps=20;
    L2diffs = zeros(nbptseps, 1);
    H1diffs = zeros(nbptseps, 1);
    epsilons=1:nbptseps; %logspace(0,1,nbptseps); % epsilons entiers pris entre 0 et 20 ici
    for neps=1:nbptseps
        epsilon=epsilons(neps);
        disp(sprintf('epsilon = %g', epsilon));
        [UU, MM, KK] = principal_dirichlet_aux( ...
            hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
            Atype, epsilon, eta); 
        if (visualisation>0)
            affiche(abs(UUh-UU), Numtri,Coorneu, ...
                    sprintf('DIFF eps : %g', epsilon));
            %plot(Coorneu(:,1)', abs(UUh-UU)', '.');
            if(visualisation>1)
                affiche(UU, Numtri,Coorneu, ...
                        sprintf('eps : %g', epsilon));
            end
        end
        L2diffs(neps) = (UU-UUh)'*MMh*(UU-UUh) ;
        H1diffs(neps) = (UU-UUh)'*KKh*(UU-UUh);
    end
    
    %%% PAREIL AVEC UN MAILLAGE PLUS GROSSIER %%%
    msh=15;
    disp(sprintf('Solution homogène'));
    [UUh, MMh, KKh] = ...
        principal_dirichlet_aux(hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
                                Atype, -1, eta); 

    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(namemshs{msh});

    L2_2_diffs = zeros(nbptseps, 1);
    H1_2_diffs = zeros(nbptseps, 1);
    epsilons=1:nbptseps; %logspace(0,1,nbptseps);
    for neps=1:nbptseps
        epsilon=epsilons(neps);
        disp(sprintf('epsilon = %g', epsilon));
        [UU, MM, KK] = principal_dirichlet_aux( ...
            hs(msh), namemshs{msh}, namemshs_pbcell{msh_pbcell}, ...
            Atype, epsilon, eta); 
        if (visualisation>0)
            affiche(abs(UUh-UU), Numtri,Coorneu, ...
                    sprintf('DIFF eps : %g', epsilon));
            %plot(Coorneu(:,1)', abs(UUh-UU)', '.');
            if(visualisation>1)
                affiche(UU, Numtri,Coorneu, ...
                        sprintf('eps : %g', epsilon));
            end
        end
        L2_2_diffs(neps) = (UU-UUh)'*MMh*(UU-UUh) ;
        H1_2_diffs(neps) = (UU-UUh)'*KKh*(UU-UUh);
    end
    
    figure();
    loglog(epsilons, L2diffs);
    hold on;
    loglog(epsilons, L2_2_diffs);
    hold off;
    grid on;
    title('L2');
    figure();
        loglog(epsilons, H1diffs);
    hold on;
    loglog(epsilons, H1_2_diffs);
    hold off;
    grid on;
    title('H1');
    
    
    %%% SECONDE EXPÉRIENCE %%%
elseif (vareta)
    nbeta=20;
    etas=logspace(-6,0,nbeta);
    norms_Ah = zeros(nbeta,1);
    for eta=1:nbeta 
        eta
        Ah = calc_A_homo(namemshs_pbcell{msh_pbcell}, Atype, eta);
        norms_Ah(eta) = norm(Ah);
    end
    figure();
    loglog(etas, norms_Ah);
    xlabel('eta');
    ylabel('norme Ah');
    
elseif (varmshbcell)
    nbmsh=10;
    mshmin=10;
    mshmax=mshmin+nbmsh;
    for msh=mshmin;mshmax
        Ah = calc_A_homo(namemsh{msh_pbcell}, Atype, eta);
        norms_Ah(eta) = norm(Ah);
    end
    figure();
    loglog(1./etas, norms_Ah);
    xlabel('1/eta');
    ylabel('norme Ah');
end
