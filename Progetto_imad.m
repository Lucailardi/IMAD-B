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


%% --------------------------------------------------------------------- %%
%                      visualizzazione grafica
% ---------------------------------------------------------------------- %%
% Utlizziamo diverse tipologie di grafico per visualizzare i dati e
% identificare visivamente eventuali outliers

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

x = [2:1:16];                     % calcolo dei coefficienti di correlazione
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


%% --------------------------------------------------------------------- %%
% Gestione degli Outliers
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


%% --------------------------------------------------------------------- %%
% randomizzazione e divisione del dataset 70-30
rng(1)                          % fisso il seed
ind = randperm(size(dati_puliti,1));   % permuto gli indici in maniera casuale

n_val = round(30*size(dati_puliti,1)/100);            % il 30% va nel set di validazione
n_iden = length(dati_puliti(:, 1)) - n_val;
dati_val = dati_puliti(ind([1:n_val]),:);
dati_iden = dati_puliti(ind([1+n_val:end]),:);

%% --------------------------------------------------------------------- %%
% plot di  Grasso corporeo vs IMC
figure(5)
plot(dati_iden(:, 15), dati_iden(:, 1), 'x');
hold on
plot(dati_val(:, 15), dati_val(:, 1), 'xr');
grid on;
title("Grasso corporeo vs IMC");
xlabel("IMC [kg*m\^(-2)]");
ylabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");

%% --------------------------------------------------------------------- %%
% plot di  Grasso corporeo vs Età
figure(6)
plot(dati_iden(:, 2), dati_iden(:, 1), 'x'); % equivalentemente si può usare scatter(arg1, arg2, [], 'colore')
hold on
plot(dati_val(:, 2), dati_val(:, 1), 'xr');
grid on;
title("Grasso corporeo vs Età");
xlabel("Età [anni]");
ylabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");

%% --------------------------------------------------------------------- %%
% plot di  grasso corporeo vs IMC e Età
figure(7)
plot3(dati_iden(:, 15), dati_iden(:, 2), dati_iden(:, 1), 'o');
hold on
plot3(dati_val(:, 15), dati_val(:,2), dati_val(:, 1), 'or');
grid on;
title("Grasso corporeo vs IMC e Età");
xlabel("IMC [kg*m\^(-2)]");
ylabel("Età [anni]");
zlabel("Grasso corporeo [kg]");
legend("Dati di identificazione", "Dati di validazione");