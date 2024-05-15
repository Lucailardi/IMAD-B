%% ------------------------------------
% Variabili utili
% --------------------
% grasso
grasso_corp_iden = dati_iden(:, 1);
grasso_corp_val = dati_val(:, 1);

% IMC
IMC_iden = dati_iden(:, 15);
IMC_val = dati_val(:, 15);

% eta
eta_iden = dati_iden(:, 2);
eta_val = dati_val(:, 2);

% griglia
IMC_grid = linspace(0, 45, 100)';
eta_grid = linspace(15, 90, 100)';

% Variabili utili
n = length(dati_iden(:,1));
nVal = length(dati_iden);
n_iden = length(dati_puliti(:, 1)) - n_val;


%% --------------------------------------------------------------------- %%
%              Modelli polinomiali di Grasso corporeo vs IMC
% -----------------------------------------------------------------------

% plot di  Grasso corporeo vs IMC
figure('Name','Grasso corporeo vs IMC') 
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
grid on;
title("Grasso corporeo vs IMC");
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");

% --------------------------------
%       Modello costante
% --------------------------------
phi0 = [ones(n_iden, 1)];
[thetaLS0, ~] = lscov(phi0, grasso_corp_iden);

q0 = length(thetaLS0);

ystima0 = phi0*thetaLS0;
epsilon0 = grasso_corp_iden-ystima0;
SSR0 = epsilon0'*epsilon0;

phi_grid = ones(length(IMC_grid), 1);
curva0 = phi_grid*thetaLS0;
hold on
plot(IMC_grid, curva0)

% --------------------------------
%        Modello affine
% --------------------------------
phi1 = [ones(n_iden, 1), IMC_iden];
[thetaLS1, ~] = lscov(phi1, grasso_corp_iden);

q1 = length(thetaLS1);

ystima1 = phi1*thetaLS1;
epsilon1 = grasso_corp_iden-ystima1;
SSR1 = epsilon1'*epsilon1;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid] ;
curva1 = phi_grid*thetaLS1;
hold on
plot(IMC_grid, curva1)

% --------------------------------
%      Modello quadratico
% --------------------------------
phi2 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2];
[thetaLS2, ~] = lscov(phi2, grasso_corp_iden);

q2 = length(thetaLS2);

ystima2 = phi2*thetaLS2;
epsilon2 = grasso_corp_iden-ystima2;
SSR2 = epsilon2'*epsilon2;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2] ;
curva2 = phi_grid*thetaLS2;
hold on
plot(IMC_grid, curva2)

% --------------------------------
%          Modello cubico
% --------------------------------
phi3 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2, IMC_iden.^3];
[thetaLS3, ~] = lscov(phi3, grasso_corp_iden);

q3 = length(thetaLS3);

ystima3 = phi3*thetaLS3;
epsilon3 = grasso_corp_iden-ystima3;
SSR3 = epsilon3'*epsilon3;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2, IMC_grid.^3] ;
curva3 = phi_grid*thetaLS3;
hold on
plot(IMC_grid, curva3)

% --------------------------------
%         Modello grado 4
% --------------------------------
phi4 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2, IMC_iden.^3, IMC_iden.^4];
[thetaLS4, ~] = lscov(phi4, grasso_corp_iden);

q4 = length(thetaLS4);

ystima4 = phi4*thetaLS4;
epsilon4 = grasso_corp_iden-ystima4;
SSR4 = epsilon4'*epsilon4;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2, IMC_grid.^3, IMC_grid.^4] ;
curva4 = phi_grid*thetaLS4;

hold on
xlim([15, 40]);
ylim([0, 50]);
plot(IMC_grid, curva4)

legend("Dati di identificazione", "Dati di validazione", "grado 0", "grado 1", "grado 2", "grado 3", "grado 4", 'Location', 'northwest')

%% --------------------------------
%        Plot diviso
% --------------------------------
figure('Name','Modelli Polinomiali IMC')

subplot 251
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(IMC_grid, curva0)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
legend('dati di identificazione')
title('grado 0')

subplot 256
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
hold on
plot(IMC_grid, curva0)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
legend('dati di validazione')
title('grado 0')

subplot 252
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(IMC_grid, curva1)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 1')

subplot 257
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
hold on
plot(IMC_grid, curva1)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 1')

subplot 253
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(IMC_grid, curva2)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 2')

subplot 258
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
hold on
plot(IMC_grid, curva2)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 2')

subplot 254
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(IMC_grid, curva3)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 3')

subplot 259
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
hold on
plot(IMC_grid, curva3)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 3')

subplot 255
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(IMC_grid, curva4)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 4')

subplot(2,5,10)
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
hold on
plot(IMC_grid, curva4)
grid on
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
xlim([15, 40]);
ylim([0, 50]);
title('grado 4')


sgtitle('Modelli polinomiali - grasso vs IMC')

%% ======================================================= %%
% FPE, AIC, MLE

% Polinomio di primo grado
FPE0 = (n+q0)/(n-q0)*SSR0;
AIC0 = 2*q0/n + log(SSR0);
MDL0 = log(n)*q0/n + log(SSR0);

% Polinomio di primo grado
FPE1 = (n+q1)/(n-q1)*SSR1;
AIC1 = 2*q1/n + log(SSR1);
MDL1 = log(n)*q1/n + log(SSR1);

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


% ----------------------------------------------------------------
% AIC minimo è quello del modello 3
clc
disp('Modelli polinomiali di IMC vs Grasso corporeo')
disp('---------------------------------------------------')
AIC = [AIC0 AIC1 AIC2 AIC3 AIC4];
k = min(AIC);
AIC_min = find(AIC == k);
disp(['AIC minimo è quello associato al modello ', num2str(AIC_min-1)]);

% FPE minimo è quello del modello 3
FPE = [FPE0 FPE1 FPE2 FPE3 FPE4];
clear k
k = min(FPE);
FPE_min = find(FPE == k);
disp(['FPE minimo è quello associato al modello ', num2str(FPE_min-1)]);

% MDL minimo è quello del modello 3
MDL = [MDL0 MDL1 MDL2 MDL3 MDL4];
clear k
k = min(MDL);
MDL_min = find(MDL == k);
disp(['MDL minimo è quello associato al modello ', num2str(MDL_min-1)]);


%% ======================================================= %%
% Fisher ci dava come migliore il modello 2
% Test F
alpha = 0.05;
disp('---------------------------------------------------')

% Confronto fra polinomio di grado 0 e di grado 1
f_alpha1 = finv(1-alpha, 1, n-q1);
f1 = (n-q1)*(SSR0 - SSR1)/SSR1;
disp('Test di fisher migliore è associato al modello 1');


if f1 < f_alpha1
    % Confronto fra polinomio di grado 1 e di grado 2
    f_alpha2 = finv(1-alpha, 1, n-q2);
    f2 = (n-q2)*(SSR1 - SSR2)/SSR2;
    disp('Test di fisher migliore è associato al modello 2');

    if f2 < f_alpha2
    % f2 < f_alpha2 quindi procediamo

    % Confronto fra polinomio di grado 2 e di grado 3
    f_alpha3 = finv(1-alpha, 1, n-q3);
    f3 = (n-q3)*(SSR2 - SSR3)/SSR3;
    disp('Test di fisher migliore è associato al modello 3');

        if f3 < f_alpha3
            % f3 < f_alpha3 quindi procediamo
            
            % Confronto fra polinomio di grado 3 e di grado 4
            f_alpha4 = finv(1-alpha, 1, n-q4);
            f4 = (n-q4)*(SSR3 - SSR4)/SSR4;
            disp('Test di fisher migliore è associato al modello 4');
            
            
        end 
    end 

end 

%% =================================================================== %%
%      Cross Validazione
% =================================================================== %%

% Polinomio grado 0
phiVal0 = [ones(n_val, 1)];
ystimaVal0 = phiVal0*thetaLS0;
epsilonVal0 = grasso_corp_val - ystimaVal0;
SSRVal0 = epsilonVal0'*epsilonVal0;

% Polinomio grado 1
phiVal1 = [ones(n_val, 1) IMC_val];
ystimaVal1 = phiVal1*thetaLS1;
epsilonVal1 = grasso_corp_val - ystimaVal1;
SSRVal1 = epsilonVal1'*epsilonVal1;

% Polinomio di grado 2
phiVal2 = [ones(n_val, 1), IMC_val, IMC_val.^2];
ystimaVal2 = phiVal2*thetaLS2;
epsilonVal2 = grasso_corp_val - ystimaVal2;
SSRVal2 = epsilonVal2'*epsilonVal2;

% Polinomio di grado 3
phiVal3 = [ones(n_val, 1), IMC_val, IMC_val.^2, IMC_val.^3];
ystimaVal3 = phiVal3*thetaLS3;
epsilonVal3 = grasso_corp_val - ystimaVal3;
SSRVal3 = epsilonVal3'*epsilonVal3;

% Polinomio di grado 4
phiVal4 = [ones(n_val, 1), IMC_val, IMC_val.^2, IMC_val.^3, IMC_val.^4];
ystimaVal4 = phiVal4*thetaLS4;
epsilonVal4 = grasso_corp_val - ystimaVal4;
SSRVal4 = epsilonVal4'*epsilonVal4;

% SSRv minimizzato dal modello di grado 2, il modello di grado 1 segue
SSRv = [SSRVal0 SSRVal1 SSRVal2 SSRVal3 SSRVal4];
clear k
k = min(SSRv);
SSRv_min = find(SSRv == k);
disp(' ')
disp('---------------------------------------------------')
disp('    Crossvalidazione')
disp('---------------------------------------------------')
disp(['SSRv minimizzato dal modello di grado ', num2str(SSRv_min-1)]);


%% ------------------------------------------------
% pulizia delle variabili
clear AIC AIC1 AIC2 AIC3 AIC4 
clear FPE FPE1 FPE2 FPE3 FPE4 
clear MDL MDL1 MDL2 MDL3 MDL4 
clear k q0 q1 q2 q3 q4

clear epsilon0 epsilon1 epsilon2 epsilon3 epsilon4 
clear i j phi0 phi1 phi2 phi3 phi4
clear SSR0 SSR1 SSR2 SSR3 SSR4
clear stddev0 stddev1 stddev2 stddev3 stddev4
clear thetaLS0 thetaLS1 thetaLS2 thetaLS3 thetaLS4
clear ystima0 ystima1 ystima2 ystima3 ystima4
clear curva0 curva1 curva2 curva3 curva4


