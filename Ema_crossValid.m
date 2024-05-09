%% -------------------------------------------------------------------- %%
%                      PROGETTO IMAD.B 2023-2024 
% --------------------------------------------------------------------- %%
% Antonini Emanuele Elia
% Bianchi Vittoria
% Caiulo Alessio
% Ilardi Luca

clear all
close all
clc


%% Crossvalidazione devo usare dati_val
eta1=dati_val(:,2) ;
IMC1 =dati_val(:,15);
grassocorporeo1=dati_val(:,1);
figure(14)
plot3(eta1, IMC1, grassocorporeo1, 'o')
grid on
title("Grasso in funzione del numero di ETA e della IMC")
xlabel("grassocorporeo")
ylabel("ETA")
zlabel("ICM")
hold on 
plot3(eta1, IMC1,grassocorporeo1, 'rx')


%% Polinomio di grado 1
phiVal1 = [ones(nVal, 1), eta1, IMC1];
ystimaVal1 = phiVal1*thetaLS1;
epsilonVal1 = grassocorporeo1- ystimaVal1;
SSRVal1 = epsilonVal1'*epsilonVal1;



% %% Polinomio di grado 2
% phiVal2 = [ones(nVal, 1), giriValidazione, pressioneValidazione, giriValidazione.^2, pressioneValidazione.^2, giriValidazione.*pressioneValidazione];
% ystimaVal2 = phiVal2*thetaLS2;
% epsilonVal2 = rendimentoValidazione - ystimaVal2;
% SSRVal2 = epsilonVal2'*epsilonVal2;
% 
% %% Polinomio di grado 3
% phiVal3 = [ones(nVal, 1), giriValidazione, pressioneValidazione, giriValidazione.^2, pressioneValidazione.^2, giriValidazione.*pressioneValidazione, ...
%     giriValidazione.^3, pressioneValidazione.^3, (giriValidazione.^2).pressioneValidazione, giriValidazione.(pressioneValidazione.^2)];
% ystimaVal3 = phiVal3*thetaLS3;
% epsilonVal3 = rendimentoValidazione - ystimaVal3;
% SSRVal3 = epsilonVal3'*epsilonVal3;
% 
% %% Polinomio di grado 4
% phiVal4 = [ones(nVal, 1), giriValidazione, pressioneValidazione, giriValidazione.^2, pressioneValidazione.^2, giriValidazione.*pressioneValidazione, ...
%     giriValidazione.^3, pressioneValidazione.^3, (giriValidazione.^2).pressioneValidazione, giriValidazione.(pressioneValidazione.^2), ...
%     giriValidazione.^4, pressioneValidazione.^4, (giriValidazione.^3).pressioneValidazione, giriValidazione.(pressioneValidazione.^3), (giriValidazione.^2).*pressioneValidazione.^2];
% ystimaVal4 = phiVal4*thetaLS4;
% epsilonVal4 = rendimentoValidazione - ystimaVal4;
% SSRVal4 = epsilonVal4'*epsilonVal4;
% 
% %% Polinomio di grado 4
% phiVal5 = [ones(nVal, 1), giriValidazione, pressioneValidazione, giriValidazione.^2, pressioneValidazione.^2, giriValidazione.*pressioneValidazione, ...
%     giriValidazione.^3, pressioneValidazione.^3, (giriValidazione.^2).pressioneValidazione, giriValidazione.(pressioneValidazione.^2), ...
%     giriValidazione.^4, pressioneValidazione.^4, (giriValidazione.^3).pressioneValidazione, giriValidazione.(pressioneValidazione.^3), (giriValidazione.^2).*pressioneValidazione.^2, ...
%     giriValidazione.^5, pressioneValidazione.^5, (giriValidazione.^4).pressioneValidazione, giriValidazione.(pressioneValidazione.^4), (giriValidazione.^3).*(pressioneValidazione.^2), ...
%     (giriValidazione.^2).*(pressioneValidazione.^3)];
% ystimaVal5 = phiVal5*thetaLS5;
% epsilonVal5 = rendimentoValidazione - ystimaVal5;
% SSRVal5 = epsilonVal5'*epsilonVal5;
% 
% %% Fine
% % Poichè SSRVal minore è quello del modello di terzo grado, plottiamo i
% % dati di validaizone su di esso
% figure(15)
% hold on
% mesh(G, P, superficie3_matrix)
% legend("dati di training", "dati di validazione", "polinomio di grado 3")
% 



% decsione modello di regressione stepwise regression
%% Variabili utili




counter1 = 1;
counter2 = 1;

%normalize data
[m n] = size(dati_iden);
numvar = n-1;
y = dati_iden(:,n);
x = dati_iden(:,1:numvar);
% %Check for interaction
% for i = 1:numvar
% for j = counter1:numvar
% x2(:,counter2) = x(:,i).*x(:,j);
% counter2 = counter2+1;
% end
% counter1 = counter1+1;
% end
%
