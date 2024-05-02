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
load dati_iden.mat
load dati_val.mat

eta = dati_iden(:,2) ;
IMC = dati_iden(:,15);
grassocorporeo = dati_iden(:,1);

% -----------------------------------------------
% plot 3D
figure('Name','Grasso corporeo vs ICM ed Età')
plot3(IMC, eta, grassocorporeo, 'o');
hold on
plot3(dati_val(:, 15), dati_val(:,2), dati_val(:, 1), 'or');
grid on;
title("Grasso corporeo vs IMC e Età");
xlabel("IMC [kg*m\^(-2)]");
ylabel("Età [anni]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");


%% ----------------------
% Variabili utili
n = length(grassocorporeo);
nVal = length(dati_iden);

% Per poter plottare i modelli ottenuti
eta_grid = linspace(0, 100, 100);
IMC_grid = linspace(0, 50, 100);
[G, P] = meshgrid(eta_grid, IMC_grid);

% Per incolonnare in un vettore unico le colonne della matrice si usa (:)
G_vett = G(:);    % età
P_vett = P(:);    % ICM


%% ----------------------
%   Polinomio grado 1
% ----------------------

phi1 = [ones(n, 1) eta IMC];
[thetaLS1, std_theta1] = lscov(phi1, grassocorporeo);

q1 = length(thetaLS1);

ystima1 = phi1*thetaLS1;
epsilon1 = grassocorporeo-ystima1;
SSR1 = epsilon1'*epsilon1;

phi1_grid = [ones(length(G_vett), 1) G_vett P_vett];
superficie1 = phi1_grid*thetaLS1;
superficie1_matrix = reshape(superficie1, size(G));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 1')
plot3(eta,IMC,grassocorporeo, 'o')
hold on
mesh(G,P, superficie1_matrix)
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("dati", "polinomio di grado 1")


%% ----------------------
%   Polinomio grado 2
% ----------------------
phi2 = [ones(n, 1), eta, IMC, eta.^2, IMC.^2, eta.*IMC];
[thetaLS2, std_theta2] = lscov(phi2, grassocorporeo);

q2 = length(thetaLS2);

ystima2 = phi2*thetaLS2;
epsilon2 = grassocorporeo-ystima2;
SSR2 = epsilon2'*epsilon2;

phi2_grid = [ones(length(G_vett), 1), G_vett, P_vett, G_vett.^2, P_vett.^2, G_vett.*P_vett];
superficie2 = phi2_grid*thetaLS2;
superficie2_matrix = reshape(superficie2, size(G));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 2')
plot3(eta,IMC,grassocorporeo, 'o')
hold on
mesh(G,P,superficie2_matrix)
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("dati", "polinomio di grado 2")


%% ----------------------
%   Polinomio grado 3
% ----------------------
phi3 = [ones(n, 1), eta, IMC, eta.^2, IMC.^2, eta.*IMC, eta.^3, IMC.^3, (eta.^2).*IMC, eta.*(IMC.^2)];
[thetaLS3, std_theta3] = lscov(phi3, grassocorporeo);

q3 = length(thetaLS3);

ystima3 = phi3*thetaLS3;
epsilon3 = grassocorporeo-ystima3;
SSR3 = epsilon3'*epsilon3;

phi3_grid = [ones(length(G_vett), 1), G_vett, P_vett, G_vett.^2, P_vett.^2, G_vett.*P_vett, G_vett.^3, P_vett.^3, (G_vett.^2).*P_vett, G_vett.*(P_vett.^2)];
superficie3 = phi3_grid*thetaLS3;
superficie3_matrix = reshape(superficie3, size(G));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 3')
plot3(eta, IMC, grassocorporeo, 'o')
hold on
mesh(G, P, superficie3_matrix)
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("dati", "polinomio di grado 3")


%% ----------------------
%   Polinomio grado 4
% ----------------------
phi4 = [ones(n, 1),  eta, IMC, eta.^2, IMC.^2, eta.*IMC, eta.^3, IMC.^3, (eta.^2).*IMC, eta.*(IMC.^2),eta.^4, IMC.^4, (eta.^3).*IMC, eta.*(IMC.^3), (eta.^2).*(IMC.^2)];
[thetaLS4, std_theta4] = lscov(phi4, grassocorporeo);

q4 = length(thetaLS4);

ystima4 = phi4*thetaLS4;
epsilon4 = grassocorporeo-ystima4;
SSR4 = epsilon4'*epsilon4;

phi4_grid = [ones(length(G_vett), 1), G_vett, P_vett, G_vett.^2, P_vett.^2, G_vett.*P_vett, G_vett.^3, P_vett.^3, (G_vett.^2).*P_vett, G_vett.*(P_vett.^2), ...
    G_vett.^4, P_vett.^4, (G_vett.^3).*P_vett, G_vett.*(P_vett.^3), (G_vett.^2).*(P_vett.^2)];
superficie4 = phi4_grid*thetaLS4;
superficie4_matrix = reshape(superficie4, size(G));

% ---------------------------------------
% plot
figure('Name','Polinomio di grado 4')
plot3(eta, IMC, grassocorporeo, 'o')
hold on
mesh(G, P, superficie4_matrix)
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("Età [anni]");
ylabel("IMC [kg*m\^(-2)]");
zlabel("Grasso corporeo [kg]");
legend("dati", "polinomio di grado 4")


%% ----------------------
%   Polinomio grado 5
% ----------------------

phi5 = [ones(n, 1), eta, IMC, eta.^2, IMC.^2, eta.*IMC, eta.^3, IMC.^3, (eta.^2).*IMC, eta.*(IMC.^2), ...
    eta.^4, IMC.^4, (eta.^3).*IMC, eta.*(IMC.^3), (eta.^2).*(IMC.^2), eta.^5, IMC.^5, (eta.^4).*IMC, eta.*(IMC.^4), ...
    (eta.^3).*(IMC.^2), (eta.^2).*(IMC.^3)];
[thetaLS5, std_theta5] = lscov(phi5, grassocorporeo);

q5 = length(thetaLS5);

ystima5 = phi5*thetaLS5;
epsilon5 = grassocorporeo-ystima5;
SSR5 = epsilon5'*epsilon5;

phi5_grid = [ones(length(G_vett), 1), G_vett, P_vett, G_vett.^2, P_vett.^2, G_vett.*P_vett, G_vett.^3, P_vett.^3, (G_vett.^2).*P_vett, G_vett.*(P_vett.^2), ...
    G_vett.^4, P_vett.^4, (G_vett.^3).*P_vett, G_vett.*(P_vett.^3), (G_vett.^2).*(P_vett.^2), G_vett.^5, P_vett.^5, (G_vett.^4).*P_vett, G_vett.*(P_vett.^4), ...
    (G_vett.^3).*(P_vett.^2), (G_vett.^2).*(P_vett.^3)];
superficie5 = phi5_grid*thetaLS5;
superficie5_matrix = reshape(superficie5, size(G));

figure('Name','Polinomio di grado 5')
plot3(eta, IMC, grassocorporeo, 'o')
hold on
mesh(G, P, superficie5_matrix)
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("peso corporeo")
ylabel("ETA")
zlabel("IMC")
legend("dati", "polinomio di grado 5")


%% ======================================================= %%
% (VB) !!!!!! ricopiare questa parte nel main per i 
% modelli con variabile singola

% Test F
alpha = 0.05;

% Confronto fra polinomio di grado 1 e di grado 2
f_alpha2 = finv(1-alpha, 1, n-q2);
f2 = (n-q2)*(SSR1 - SSR2)/SSR2;

if f2 < f_alpha2
    % f2 < f_alpha2 quindi procediamo

    % Confronto fra polinomio di grado 2 e di grado 3
    f_alpha3 = finv(1-alpha, 1, n-q3);
    f3 = (n-q3)*(SSR2 - SSR3)/SSR3;

    if f3 < f_alpha3
        % f3 < f_alpha3 quindi procediamo
        
        % Confronto fra polinomio di grado 3 e di grado 4
        f_alpha4 = finv(1-alpha, 1, n-q4);
        f4 = (n-q4)*(SSR3 - SSR4)/SSR4;
        
        if f4 < f_alpha4
            % f4 < f_alpha4 quindi procediamo
            
            % Confronto fra polinomio di grado 4 e di grado 5
            f_alpha5 = finv(1-alpha, 1, n-q5);
            f5 = (n-q5)*(SSR4 - SSR5)/SSR5;

        end 
    end 
end 


%% ======================================================= %%
% (VB) !!!!!! ricopiare questa parte nel main per i 
% modelli con variabile singola

% FPE, AIC, MLE

% Polinomio di primo grado
FPE1 = (n+q1)/(n-q1)*SSR1;
AIC1 = 2*q1/n + log(SSR1);
MDL1 = log(n)*q1/n + log(SSR1)

% Polinomio di secondo grado
FPE2 = (n+q2)/(n-q2)*SSR2;
AIC2 = 2*q2/n + log(SSR2);
MDL2 = log(n)*q2/n + log(SSR2);

% Polinomio di terzo grado
FPE3 = (n+q3)/(n-q3)*SSR3;
AIC3 = 2*q3/n + log(SSR3);
MDL3 = log(n)*q3/n + log(SSR3);

% Polinomio di quarto grado
FPE4 = (n+q4)/(n-q4)*SSR4;
AIC4 = 2*q4/n + log(SSR4);
MDL4 = log(n)*q4/n + log(SSR4);

% Polinomio di quinto grado
FPE5 = (n+q5)/(n-q5)*SSR5;
AIC5 = 2*q5/n + log(SSR5);
MDL5 = log(n)*q5/n + log(SSR5);

% AIC minimo è quello del modello 3
% FPE minimo è quello del modello 4
% MDL minimo è quello del modello 3

% Fisher ci dava come migliore il modello 4
