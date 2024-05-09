% Questa parte è separata dal resto delle elaborazioni per i modelli
% monovariati per evitare di sovrapporsi sul file principale.
% Per questo manca la parte di caricamento e pulizia dati, e va eseguito
% con il workspace già pieno di tutto il necessario.
% Andrà poi messo in coda ai modelli polinomiali monovariati

%% Grasso corporeo vs IMC
%% Test F
alpha = 0.05;

% Confronto fra polinomio di grado 0 e di grado 1
f_alpha1 = finv(1-alpha, 1, n_iden-q1);
f1 = (n_iden-q1)*(SSR0 - SSR1)/SSR1;

if f1 < f_alpha1
    % f1 < f_alpha1 quindi procediamo

    % Confronto fra polinomio di grado 1 e di grado 2
    f_alpha2 = finv(1-alpha, 1, n_iden-q2);
    f2 = (n_iden-q2)*(SSR1 - SSR2)/SSR2;

    if f2 < f_alpha2
        % f2 < f_alpha2 quindi procediamo
        
        % Confronto fra polinomio di grado 2 e di grado 3
        f_alpha3 = finv(1-alpha, 1, n_iden-q3);
        f3 = (n_iden-q3)*(SSR2 - SSR3)/SSR3;
        
        if f3 < f_alpha3
            % f3 < f_alpha3 quindi procediamo
            
            % Confronto fra polinomio di grado 4 e di grado 5
            f_alpha4 = finv(1-alpha, 1, n_iden-q4);
            f4 = (n_iden-q4)*(SSR3 - SSR4)/SSR4;

        end 
    end 
end 

% Modello consigliato è di grado 0? Strano

%% Test oggettivi
% Polinomio di grado 0
FPE0 = (n_iden+q0)/(n_iden-q0)*SSR0;
AIC0 = 2*q0/n_iden + log(SSR0);
MDL0 = log(n_iden)*q0/n_iden + log(SSR0);

% Polinomio di primo grado
FPE1 = (n_iden+q1)/(n_iden-q1)*SSR1;
AIC1 = 2*q1/n_iden + log(SSR1);
MDL1 = log(n_iden)*q1/n_iden + log(SSR1);

% Polinomio di secondo grado
FPE2 = (n_iden+q2)/(n_iden-q2)*SSR2;
AIC2 = 2*q2/n_iden + log(SSR2);
MDL2 = log(n_iden)*q2/n_iden + log(SSR2);

% Polinomio di terzo grado
FPE3 = (n_iden+q3)/(n_iden-q3)*SSR3;
AIC3 = 2*q3/n_iden + log(SSR3);
MDL3 = log(n_iden)*q3/n_iden + log(SSR3);

% Polinomio di quarto grado
FPE4 = (n_iden+q4)/(n_iden-q4)*SSR4;
AIC4 = 2*q4/n_iden + log(SSR4);
MDL4 = log(n_iden)*q4/n_iden + log(SSR4);

% FPE minimizzato da modello grado 1
% AIC minimizzato da modello grado 1
% MDL minimizzato da modello grado 1

%% Crossvalidazione
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

%% Grasso corporeo vs Età
% Test F
alpha = 0.05;

% Confronto fra polinomio di grado 0 e di grado 1
f_alpha1 = finv(1-alpha, 1, n_iden-q1);
f1 = (n_iden-q1)*(SSR0 - SSR1)/SSR1;

if f1 < f_alpha1
    % f1 < f_alpha1 quindi procediamo

    % Confronto fra polinomio di grado 1 e di grado 2
    f_alpha2 = finv(1-alpha, 1, n_iden-q2);
    f2 = (n_iden-q2)*(SSR1 - SSR2)/SSR2;

    if f2 < f_alpha2
        % f2 < f_alpha2 quindi procediamo
        
        % Confronto fra polinomio di grado 2 e di grado 3
        f_alpha3 = finv(1-alpha, 1, n_iden-q3);
        f3 = (n_iden-q3)*(SSR2 - SSR3)/SSR3;
        
        if f3 < f_alpha3
            % f3 < f_alpha3 quindi procediamo
            
            % Confronto fra polinomio di grado 4 e di grado 5
            f_alpha4 = finv(1-alpha, 1, n_iden-q4);
            f4 = (n_iden-q4)*(SSR3 - SSR4)/SSR4;

        end 
    end 
end 

% Fisher dice modello grado 0 ????

%% Test oggettivi
% Polinomio di grado 0
FPE0 = (n_iden+q0)/(n_iden-q0)*SSR0;
AIC0 = 2*q0/n_iden + log(SSR0);
MDL0 = log(n_iden)*q0/n_iden + log(SSR0);

% Polinomio di primo grado
FPE1 = (n_iden+q1)/(n_iden-q1)*SSR1;
AIC1 = 2*q1/n_iden + log(SSR1);
MDL1 = log(n_iden)*q1/n_iden + log(SSR1);

% Polinomio di secondo grado
FPE2 = (n_iden+q2)/(n_iden-q2)*SSR2;
AIC2 = 2*q2/n_iden + log(SSR2);
MDL2 = log(n_iden)*q2/n_iden + log(SSR2);

% Polinomio di terzo grado
FPE3 = (n_iden+q3)/(n_iden-q3)*SSR3;
AIC3 = 2*q3/n_iden + log(SSR3);
MDL3 = log(n_iden)*q3/n_iden + log(SSR3);

% Polinomio di quarto grado
FPE4 = (n_iden+q4)/(n_iden-q4)*SSR4;
AIC4 = 2*q4/n_iden + log(SSR4);
MDL4 = log(n_iden)*q4/n_iden + log(SSR4);

% FPE minimizzato da modello grado 3
% AIC minimizzato da modello grado 3
% MDL minimizzato da modello grado 1

%% Crossvalidazione
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

% SSRv minimizzato dal modello di grado 0, il modello di grado 1 segue
