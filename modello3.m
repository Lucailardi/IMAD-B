%% --------------------------------------------------------------------- %%
%        Modello polinomiale di Grasso corporeo vs ICM vs Età
% -----------------------------------------------------------------------
% IDENTIFICAZIONE DEL MODELLO:
% Avendo due variabili indipendenti la matrice phi deve tenere conto di
% ogni monomio combinazione delle due variabili che abbia come grado quello
% del polinomio preso a modello: per esempio per il grado 1 la phi avrà tre
% colonne: una di 1, una con eta e una con ICM. Al grado due
% avrà sei colonne, cioè le tre precedenti più ogni combinazione di eta e
% ICM di grado due: eta^2, ICM^2, eta*ICM

% caricamento dei dati
% load dati_iden.mat
% load dati_val.mat

eta_iden = dati_iden(:,2) ;
IMC_iden = dati_iden(:,15);
grasso_corp_iden = dati_iden(:,1);
eta_val = dati_val(:,2);
IMC_val = dati_val(:,15);
grasso_corp_val = dati_val(:,1);

% -----------------------------------------------
% plot 3D
figure('Name','Grasso corporeo vs ICM ed Età')
plot3(eta_iden, IMC_iden, grasso_corp_iden, 'o');
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or');
grid on;
title("Grasso corporeo vs IMC e Età");
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");


%% ----------------------
% Variabili utili
nIden = length(grasso_corp_iden);
nVal = length(dati_val);

% Per poter plottare i modelli ottenuti
eta_grid = linspace(20, 80, 100);
IMC_grid = linspace(15, 40, 100);
[E, I] = meshgrid(eta_grid, IMC_grid);

% Per incolonnare in un vettore unico le colonne della matrice si usa (:)
E_vett = E(:);    % età
I_vett = I(:);    % ICM


%% ----------------------
%   Polinomio grado 1
% ----------------------

phi1 = [ones(nIden, 1) eta_iden IMC_iden];
[thetaLS1, std_theta1] = lscov(phi1, grasso_corp_iden);

q1 = length(thetaLS1);

ystima1 = phi1*thetaLS1;
epsilon1 = grasso_corp_iden-ystima1;
SSR1 = epsilon1'*epsilon1;

phi1_grid = [ones(length(E_vett), 1) E_vett I_vett];
superficie1 = phi1_grid*thetaLS1;
superficie1_matrix = reshape(superficie1, size(E));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 1')
plot3(eta_iden,IMC_iden,grasso_corp_iden, 'o')
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or')
hold on
mesh(E,I, superficie1_matrix)
hold on
grid on
title("Grasso corporeo in funzione di età e IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione", "polinomio di grado 1")


%% ----------------------
%   Polinomio grado 2
% ----------------------
phi2 = [ones(nIden, 1), eta_iden, IMC_iden, eta_iden.^2, IMC_iden.^2, eta_iden.*IMC_iden];
[thetaLS2, std_theta2] = lscov(phi2, grasso_corp_iden);

q2 = length(thetaLS2);

ystima2 = phi2*thetaLS2;
epsilon2 = grasso_corp_iden-ystima2;
SSR2 = epsilon2'*epsilon2;

phi2_grid = [ones(length(E_vett), 1), E_vett, I_vett, E_vett.^2, I_vett.^2, E_vett.*I_vett];
superficie2 = phi2_grid*thetaLS2;
superficie2_matrix = reshape(superficie2, size(E));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 2')
plot3(eta_iden,IMC_iden,grasso_corp_iden, 'o')
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or')
hold on
mesh(E,I,superficie2_matrix)
grid on
title("Grasso corporeo in funzione di età e IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione", "polinomio di grado 2")


%% ----------------------
%   Polinomio grado 3
% ----------------------
phi3 = [ones(nIden, 1), eta_iden, IMC_iden, eta_iden.^2, IMC_iden.^2, eta_iden.*IMC_iden, eta_iden.^3, IMC_iden.^3, (eta_iden.^2).*IMC_iden, eta_iden.*(IMC_iden.^2)];
[thetaLS3, std_theta3] = lscov(phi3, grasso_corp_iden);

q3 = length(thetaLS3);

ystima3 = phi3*thetaLS3;
epsilon3 = grasso_corp_iden-ystima3;
SSR3 = epsilon3'*epsilon3;

phi3_grid = [ones(length(E_vett), 1), E_vett, I_vett, E_vett.^2, I_vett.^2, E_vett.*I_vett, E_vett.^3, I_vett.^3, (E_vett.^2).*I_vett, E_vett.*(I_vett.^2)];
superficie3 = phi3_grid*thetaLS3;
superficie3_matrix = reshape(superficie3, size(E));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 3')
plot3(eta_iden, IMC_iden, grasso_corp_iden, 'o')
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or')
hold on
mesh(E, I, superficie3_matrix)
grid on
title("Grasso corporeo in funzione di età e IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione", "polinomio di grado 3")


%% ----------------------
%   Polinomio grado 4
% ----------------------
phi4 = [ones(nIden, 1),  eta_iden, IMC_iden, eta_iden.^2, IMC_iden.^2, eta_iden.*IMC_iden, eta_iden.^3, IMC_iden.^3, (eta_iden.^2).*IMC_iden, eta_iden.*(IMC_iden.^2),eta_iden.^4, IMC_iden.^4, (eta_iden.^3).*IMC_iden, eta_iden.*(IMC_iden.^3), (eta_iden.^2).*(IMC_iden.^2)];
[thetaLS4, std_theta4] = lscov(phi4, grasso_corp_iden);

q4 = length(thetaLS4);

ystima4 = phi4*thetaLS4;
epsilon4 = grasso_corp_iden-ystima4;
SSR4 = epsilon4'*epsilon4;

phi4_grid = [ones(length(E_vett), 1), E_vett, I_vett, E_vett.^2, I_vett.^2, E_vett.*I_vett, E_vett.^3, I_vett.^3, (E_vett.^2).*I_vett, E_vett.*(I_vett.^2), ...
    E_vett.^4, I_vett.^4, (E_vett.^3).*I_vett, E_vett.*(I_vett.^3), (E_vett.^2).*(I_vett.^2)];
superficie4 = phi4_grid*thetaLS4;
superficie4_matrix = reshape(superficie4, size(E));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 4')
plot3(eta_iden, IMC_iden, grasso_corp_iden, 'o')
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or')
hold on
mesh(E, I, superficie4_matrix)
grid on
title("Grasso corporeo in funzione di età e IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione",  "Dati di validazione", "polinomio di grado 4")


%% ----------------------
%   Polinomio grado 5
% ----------------------

phi5 = [ones(nIden, 1), eta_iden, IMC_iden, eta_iden.^2, IMC_iden.^2, eta_iden.*IMC_iden, eta_iden.^3, IMC_iden.^3, (eta_iden.^2).*IMC_iden, eta_iden.*(IMC_iden.^2), ...
    eta_iden.^4, IMC_iden.^4, (eta_iden.^3).*IMC_iden, eta_iden.*(IMC_iden.^3), (eta_iden.^2).*(IMC_iden.^2), eta_iden.^5, IMC_iden.^5, (eta_iden.^4).*IMC_iden, eta_iden.*(IMC_iden.^4), ...
    (eta_iden.^3).*(IMC_iden.^2), (eta_iden.^2).*(IMC_iden.^3)];
[thetaLS5, std_theta5] = lscov(phi5, grasso_corp_iden);

q5 = length(thetaLS5);

ystima5 = phi5*thetaLS5;
epsilon5 = grasso_corp_iden-ystima5;
SSR5 = epsilon5'*epsilon5;

phi5_grid = [ones(length(E_vett), 1), E_vett, I_vett, E_vett.^2, I_vett.^2, E_vett.*I_vett, E_vett.^3, I_vett.^3, (E_vett.^2).*I_vett, E_vett.*(I_vett.^2), ...
    E_vett.^4, I_vett.^4, (E_vett.^3).*I_vett, E_vett.*(I_vett.^3), (E_vett.^2).*(I_vett.^2), E_vett.^5, I_vett.^5, (E_vett.^4).*I_vett, E_vett.*(I_vett.^4), ...
    (E_vett.^3).*(I_vett.^2), (E_vett.^2).*(I_vett.^3)];
superficie5 = phi5_grid*thetaLS5;
superficie5_matrix = reshape(superficie5, size(E));

figure('Name','Polinomio di grado 5')
plot3(eta_iden, IMC_iden, grasso_corp_iden, 'o')
hold on
plot3(eta_val, IMC_val, grasso_corp_val, 'or')
hold on
mesh(E, I, superficie5_matrix)
grid on
title("Grasso corporeo in funzione di età e IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione", "polinomio di grado 5")


%% ======================================================= %%

% Test F
alpha = 0.05;

% Confronto fra polinomio di grado 1 e di grado 2
f_alpha2 = finv(1-alpha, 1, nIden-q2);
f2 = (nIden-q2)*(SSR1 - SSR2)/SSR2;

if f2 < f_alpha2
    % f2 < f_alpha2 quindi procediamo

    % Confronto fra polinomio di grado 2 e di grado 3
    f_alpha3 = finv(1-alpha, 1, nIden-q3);
    f3 = (nIden-q3)*(SSR2 - SSR3)/SSR3;

    if f3 < f_alpha3
        % f3 < f_alpha3 quindi procediamo
        
        % Confronto fra polinomio di grado 3 e di grado 4
        f_alpha4 = finv(1-alpha, 1, nIden-q4);
        f4 = (nIden-q4)*(SSR3 - SSR4)/SSR4;
        
        if f4 < f_alpha4
            % f4 < f_alpha4 quindi procediamo
            
            % Confronto fra polinomio di grado 4 e di grado 5
            f_alpha5 = finv(1-alpha, 1, nIden-q5);
            f5 = (nIden-q5)*(SSR4 - SSR5)/SSR5;

        end 
    end 
end 


%% ======================================================= %%

% Test oggettivi
% FPE, AIC, MLE

% Polinomio di primo grado
FPE1 = (nIden+q1)/(nIden-q1)*SSR1;
AIC1 = 2*q1/nIden + log(SSR1);
MDL1 = log(nIden)*q1/nIden + log(SSR1);

% Polinomio di secondo grado
FPE2 = (nIden+q2)/(nIden-q2)*SSR2;
AIC2 = 2*q2/nIden + log(SSR2);
MDL2 = log(nIden)*q2/nIden + log(SSR2);

% Polinomio di terzo grado
FPE3 = (nIden+q3)/(nIden-q3)*SSR3;
AIC3 = 2*q3/nIden + log(SSR3);
MDL3 = log(nIden)*q3/nIden + log(SSR3);

% Polinomio di quarto grado
FPE4 = (nIden+q4)/(nIden-q4)*SSR4;
AIC4 = 2*q4/nIden + log(SSR4);
MDL4 = log(nIden)*q4/nIden + log(SSR4);

% Polinomio di quinto grado
FPE5 = (nIden+q5)/(nIden-q5)*SSR5;
AIC5 = 2*q5/nIden + log(SSR5);
MDL5 = log(nIden)*q5/nIden + log(SSR5);

%% Crossvalidazione
% Polinomio grado 1
phiVal1 = [ones(nVal, 1) eta_val IMC_val];
ystimaVal1 = phiVal1*thetaLS1;
epsilonVal1 = grasso_corp_val - ystimaVal1;
SSRVal1 = epsilonVal1'*epsilonVal1;

% Polinomio di grado 2
phiVal2 = [ones(nVal, 1), eta_val, IMC_val, eta_val.^2, IMC_val.^2, eta_val.*IMC_val];
ystimaVal2 = phiVal2*thetaLS2;
epsilonVal2 = grasso_corp_val - ystimaVal2;
SSRVal2 = epsilonVal2'*epsilonVal2;

% Polinomio di grado 3
phiVal3 = [ones(nVal, 1), eta_val, IMC_val, eta_val.^2, IMC_val.^2, eta_val.*IMC_val, ...
    eta_val.^3, IMC_val.^3, (eta_val.^2).*IMC_val, eta_val.*(IMC_val.^2)];
ystimaVal3 = phiVal3*thetaLS3;
epsilonVal3 = grasso_corp_val - ystimaVal3;
SSRVal3 = epsilonVal3'*epsilonVal3;

% Polinomio di grado 4
phiVal4 = [ones(nVal, 1), eta_val, IMC_val, eta_val.^2, IMC_val.^2, eta_val.*IMC_val, ...
    eta_val.^3, IMC_val.^3, (eta_val.^2).*IMC_val, eta_val.*(IMC_val.^2), ...
    eta_val.^4, IMC_val.^4, (eta_val.^3).*IMC_val, eta_val.*(IMC_val.^3), (eta_val.^2).*IMC_val.^2];
ystimaVal4 = phiVal4*thetaLS4;
epsilonVal4 = grasso_corp_val - ystimaVal4;
SSRVal4 = epsilonVal4'*epsilonVal4;

% Polinomio di grado 5
phiVal5 = [ones(nVal, 1), eta_val, IMC_val, eta_val.^2, IMC_val.^2, eta_val.*IMC_val, ...
    eta_val.^3, IMC_val.^3, (eta_val.^2).*IMC_val, eta_val.*(IMC_val.^2), ...
    eta_val.^4, IMC_val.^4, (eta_val.^3).*IMC_val, eta_val.*(IMC_val.^3), (eta_val.^2).*IMC_val.^2, ...
    eta_val.^5, IMC_val.^5, (eta_val.^4).*IMC_val, eta_val.*(IMC_val.^4), (eta_val.^3).*(IMC_val.^2), ...
    (eta_val.^2).*(IMC_val.^3)];
ystimaVal5 = phiVal5*thetaLS5;
epsilonVal5 = grasso_corp_val - ystimaVal5;
SSRVal5 = epsilonVal5'*epsilonVal5;

%% Conclusione
% AIC minimo è quello del modello 3
% FPE minimo è quello del modello 3
% MDL minimo è quello del modello 1

% Fisher ci dava come migliore il modello 1

% Crossvalidazione da SSR minimo per il modello 1 (con molto margine
% rispetto agli altri soprattutto 3, 4, 5)
