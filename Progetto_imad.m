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

% load dei dati
% dati = xlsread('BodyFatData.xlsx');
% dati = readcell('BodyFatData.xlsx');
tabella = readtable('BodyFatData.xlsx');

matrice = table2array(tabella);   % converto la tabella in matrice
dati = matrice(:,[2:end]);        % matrice dei dati SENZA la variabile target
target = matrice(:,1);            % SOLO variabile target

% set del nome delle variabili
nomi = {'Eta','Peso','Altezza','Collo','Petto','Addome','Anca','Coscia',...
    'Ginocchio','Caviglia','Bicipite','Avambraccio','Polso','BMI','RVF'};
nom_target = {'Grasso'};

%% --------------------------------------------------------------------- %%
%                      visualizzazione grafica
% ---------------------------------------------------------------------- %%
% Utlizziamo diverse tipologie di grafico per visualizzare i dati e
% identificare visivamente eventuali outliers

% --------------------------------
% Grafico dei dati raccolti

figure('Name', 'Dati')
for i = 1:size(dati,2)

    subplot(4,4,i)
    plot(dati(:,i), target, '*')
    hold on
    plot(mean(dati(:,i)),mean(target(:,1)),'*r')
    xlabel([nomi{i}])
    ylabel([nom_target{1}])
    title([nomi{i}, 'vs' ,nom_target{1}])

end

% --------------------------------
% Istogrami di tutte le variabili

figure('Name','Istogramma')
for i = 1:size(dati,2)

    subplot(4,4,i)
    histogram(dati(:,i))
    title([nomi{i}])
    sgtitle('Istogramma di tutte le variabili')

end

% --------------------------------
% BoxPlot di tutte le variabili

figure('Name','BoxPlot')
for i = 1:size(dati,2)

    subplot(4,4,i)
    boxplot(dati(:,i))
    title([nomi{i}])
    sgtitle('BoxPlot di tutte le variabili')

end 

% ----------------------------------------------------------------
% !!! ATTENZIONE 
% il corrplot è commentato perche ci mette 6 ore a runnare
% ----------------------------------------------------------------
% CorrPlot con tutte le variabili (compresa quella di outcome)
% figure('Name','CorrPlot')
% corrplot(matrice)
% title('corrPlot')


% --------------------------------
% HeatMap

% calcolo della matrice di correlazione
R1 = corr(matrice,'rows','complete','Type','Kendall');  

% plot della heatmap
figure('Name','Heatmap')
h1 = heatmap(R1);
h1.XDisplayLabels = {'Grasso',nomi{1,:}};
h1.YDisplayLabels = {'Grasso',nomi{1,:}};
title('HeatMap')

% --------------------------------
% Coefficienti di correlazione
% Il range del coefficiente sta tra -1 e 1
% quindi è giusto che quello dell'altezza sia negativo

x = [2:1:size(matrice,2)];        % calcolo dei coefficienti di correlazione
cor1 = R1(x,1);
corr1 = flip(sort(cor1));         % ordino i coefficienti in senso decrescente

% estrazione dei nomi in ordine
labels = cell(1,size(dati,2));    % preistanzio

for i = 1:size(cor1)
    for j = 1:size(corr1)
        if corr1(i)==cor1(j)
            labels{1,i} = nomi{1,j};
        end 
    end 
end 

% bar plot con inddici di correlazione
figure('Name','Coefficienti di correlazione')
bar(x,corr1,'BarWidth',0.2);
set(gca,'xticklabel',labels)
title('Coefficienti di correlazione')

% pulizia delle variabili
clear cor1 corr1 R1 h1 x

%% --------------------------------------------------------------------- %%
%    Gestione degli Outliers
% --------------------------------
% Siccome si tratta di pochi dati ~ 5% allora elimino le righe
outliers = isoutlier(matrice);
Any = any(outliers,2);

n_out = size(find(Any==1),1);
perc = (n_out*100)/size(dati,1);
disp('percentuale di outliers')
disp(perc)

% pulizia da outliers
dove = find(Any==1);
dati_puliti = matrice;        % matrice dei dati ma senza outliers
dati_puliti(dove,:) = [];

% salviamo la matrice 
% save('dati_puliti.mat',"dati_puliti");

% pulizia delle variabili
clear outliers Any n_out perc dove 


%% --------------------------------------------------------------------- %%
%   divisione del dataset 70-30
% --------------------------------

% prima devo randomizzare
rng(1)                                 % fisso il seed
ind = randperm(size(dati_puliti,1));   % permuto gli indici in maniera casuale

n_val = round(30*size(dati_puliti,1)/100);            % il 30% va nel set di validazione
n_iden = length(dati_puliti(:, 1)) - n_val;
dati_val = dati_puliti(ind([1:n_val]),:);
dati_iden = dati_puliti(ind([1+n_val:end]),:);

% salviamo la matrice
% save('dati_val.mat',"dati_val")
% save('dati_iden.mat',"dati_iden")

% pulizia delle variabili
clear ind n_val 


%% --------------------------------------------------------------------- %%
%     Modelli polinomiali
% --------------------------------

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


%% --------------------------------------------------------------------- %%
%              Modelli polinomiali di Grasso corporeo vs IMC
% -----------------------------------------------------------------------

% plot di  Grasso corporeo vs IMC
figure('Name','Grasso corporeo vs IMC') 
plot(IMC_iden, grasso_corp_iden, 'x')
hold on
plot(IMC_val, grasso_corp_val, 'xr')
grid on;
title("Grasso corporeo vs IMC")
xlabel("IMC [kg*m\^(-2)]")
ylabel("Grasso corporeo [kg]")

% --------------------------------
%       Modello costante
% --------------------------------
phi0 = [ones(n_iden, 1)];
[thetaLS0, stddev0] = lscov(phi0, grasso_corp_iden);

q0 = length(thetaLS0);

ystima0 = phi0*thetaLS0;
epsilon0 = grasso_corp_iden-ystima0;
SSR0 = epsilon0'*epsilon0;

phi_grid = ones(length(IMC_grid), 1);
curva = phi_grid*thetaLS0;
hold on
plot(IMC_grid, curva)

% --------------------------------
%        Modello affine
% --------------------------------
phi1 = [ones(n_iden, 1), IMC_iden];
[thetaLS1, stddev1] = lscov(phi1, grasso_corp_iden);

q1 = length(thetaLS1);

ystima1 = phi1*thetaLS1;
epsilon1 = grasso_corp_iden-ystima1;
SSR1 = epsilon1'*epsilon1;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid] ;
curva = phi_grid*thetaLS1;
hold on
plot(IMC_grid, curva)

% --------------------------------
%      Modello quadratico
% --------------------------------
phi2 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2];
[thetaLS2, stddev2] = lscov(phi2, grasso_corp_iden);

q2 = length(thetaLS2);

ystima2 = phi2*thetaLS2;
epsilon2 = grasso_corp_iden-ystima2;
SSR2 = epsilon2'*epsilon2;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2] ;
curva = phi_grid*thetaLS2;
hold on
plot(IMC_grid, curva)

% --------------------------------
%          Modello cubico
% --------------------------------
phi3 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2, IMC_iden.^3];
[thetaLS3, stddev3] = lscov(phi3, grasso_corp_iden);

q3 = length(thetaLS3);

ystima3 = phi3*thetaLS3;
epsilon3 = grasso_corp_iden-ystima3;
SSR3 = epsilon3'*epsilon3;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2, IMC_grid.^3] ;
curva = phi_grid*thetaLS3;
hold on
plot(IMC_grid, curva)

% --------------------------------
%         Modello grado 4
% --------------------------------
phi4 = [ones(n_iden, 1), IMC_iden, IMC_iden.^2, IMC_iden.^3, IMC_iden.^4];
[thetaLS4, stddev4] = lscov(phi4, grasso_corp_iden);

q4 = length(thetaLS4);

ystima4 = phi4*thetaLS4;
epsilon4 = grasso_corp_iden-ystima4;
SSR4 = epsilon4'*epsilon4;

phi_grid = [ones(length(IMC_grid), 1), IMC_grid, IMC_grid.^2, IMC_grid.^3, IMC_grid.^4] ;
curva = phi_grid*thetaLS4;
hold on
xlim([15, 40]);
ylim([0, 50]);
plot(IMC_grid, curva)

legend("Dati di identificazione", "Dati di validazione", "grado 0", "grado 1", "grado 2", "grado 3", "grado 4", 'Location', 'northwest')

% --------------------------------
% pulizia delle variabili
clear epsilon0 epsilon1 epsilon2 epsilon3 epsilon4 
clear i j phi0 phi1 phi2 phi3 phi4
clear q0 q1 q2 q3 q4 
clear SSR0 SSR1 SSR2 SSR3 SSR4
clear stddev0 stddev1 stddev2 stddev3 stddev4
clear thetaLS0 thetaLS1 thetaLS2 thetaLS3 thetaLS4
clear ystima0 ystima1 ystima2 ystima3 ystima4
clear curva


%% --------------------------------------------------------------------- %%
%                   plot di  Grasso corporeo vs Età
% -----------------------------------------------------------------------

figure('Name','Grasso corporeo vs Età')
% equivalentemente si può usare scatter(arg1, arg2, [], 'colore')
plot(eta_iden, grasso_corp_iden, 'x'); 
hold on
plot(eta_val, grasso_corp_val, 'xr');
grid on;
title("Grasso corporeo vs Età");
xlabel("Età [anni]");
ylabel("Grasso corporeo [kg]");

% --------------------------------
%       Modello costante
% --------------------------------
phi0 = [ones(n_iden, 1)];
[thetaLS0, stddev0] = lscov(phi0, grasso_corp_iden);

q0 = length(thetaLS0);

ystima0 = phi0*thetaLS0;
epsilon0 = grasso_corp_iden-ystima0;
SSR0 = epsilon0'*epsilon0;

phi_grid = ones(length(eta_grid), 1);
curva = phi_grid*thetaLS0;
hold on
plot(eta_grid, curva)

% --------------------------------
%        Modello affine
% --------------------------------
phi1 = [ones(n_iden, 1), eta_iden];
[thetaLS1, stddev1] = lscov(phi1, grasso_corp_iden);

q1 = length(thetaLS1);

ystima1 = phi1*thetaLS1;
epsilon1 = grasso_corp_iden-ystima1;
SSR1 = epsilon1'*epsilon1;

phi_grid = [ones(length(eta_grid), 1), eta_grid] ;
curva = phi_grid*thetaLS1;
hold on
plot(eta_grid, curva)

% --------------------------------
%      Modello quadratico
% --------------------------------
phi2 = [ones(n_iden, 1), eta_iden, eta_iden.^2];
[thetaLS2, stddev2] = lscov(phi2, grasso_corp_iden);

q2 = length(thetaLS2);

ystima2 = phi2*thetaLS2;
epsilon2 = grasso_corp_iden-ystima2;
SSR2 = epsilon2'*epsilon2;

phi_grid = [ones(length(eta_grid), 1), eta_grid, eta_grid.^2] ;
curva = phi_grid*thetaLS2;
hold on
plot(eta_grid, curva)

% --------------------------------
%        Modello cubico
% --------------------------------
phi3 = [ones(n_iden, 1), eta_iden, eta_iden.^2, eta_iden.^3];
[thetaLS3, stddev3] = lscov(phi3, grasso_corp_iden);

q3 = length(thetaLS3);

ystima3 = phi3*thetaLS3;
epsilon3 = grasso_corp_iden-ystima3;
SSR3 = epsilon3'*epsilon3;

phi_grid = [ones(length(eta_grid), 1), eta_grid, eta_grid.^2, eta_grid.^3] ;
curva = phi_grid*thetaLS3;
hold on
plot(eta_grid, curva)

% --------------------------------
%       Modello grado 4
% --------------------------------
phi4 = [ones(n_iden, 1), eta_iden, eta_iden.^2, eta_iden.^3, eta_iden.^4];
[thetaLS4, stddev4] = lscov(phi4, grasso_corp_iden);

q4 = length(thetaLS4);

ystima4 = phi4*thetaLS4;
epsilon4 = grasso_corp_iden-ystima4;
SSR4 = epsilon4'*epsilon4;

phi_grid = [ones(length(eta_grid), 1), eta_grid, eta_grid.^2, eta_grid.^3, eta_grid.^4] ;
curva = phi_grid*thetaLS4;
hold on
xlim([15, 75]);
ylim([0, 45])
plot(eta_grid, curva)

legend("Dati di identificazione", "Dati di validazione", "grado 0", "grado 1", "grado 2", "grado 3", "grado 4", 'Location', 'northwest')% 

% --------------------------------
% pulizia delle variabili
clear epsilon0 epsilon1 epsilon2 epsilon3 epsilon4 
clear i j phi0 phi1 phi2 phi3 phi4
clear q0 q1 q2 q3 q4 
clear SSR0 SSR1 SSR2 SSR3 SSR4
clear stddev0 stddev1 stddev2 stddev3 stddev4
clear thetaLS0 thetaLS1 thetaLS2 thetaLS3 thetaLS4
clear ystima0 ystima1 ystima2 ystima3 ystima4
clear curva

%% --------------------------------------------------------------------- %%
% % plot di  grasso corporeo vs IMC e Età
% figure(7)
% plot3(dati_iden(:, 15), dati_iden(:, 2), dati_iden(:, 1), 'o');
% hold on
% plot3(dati_val(:, 15), dati_val(:,2), dati_val(:, 1), 'or');
% grid on;
% title("Grasso corporeo vs IMC e Età");
% xlabel("IMC [kg*m\^(-2)]");
% ylabel("Età [anni]");
% zlabel("Grasso corporeo [kg]");
% legend("Dati di identificazione", "Dati di validazione");