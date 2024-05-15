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
outliers = isoutlier(dati);
close all

% #D95319 aracione
% #7E2F8E viola
% #77AC30 verde

figure('Name', 'Dati')
for i = 1:size(dati,2)

    subplot(4,4,i)
    plot(dati(:,i), target, '*','MarkerSize',7)
    hold on
    plot(mean(dati(:,i)),mean(target(:,1)),'r*','MarkerSize',7)

    xlabel([nomi{i}])
    ylabel([nom_target{1}])
    title([nomi{i}, 'vs' ,nom_target{1}])

end
legend('dati','media')
sgtitle('Scatter plot di tutte le variabili')

k=1;
figure('Name', 'scatter - outliers')
for i = 1:size(dati,2)

    subplot(4,4,i)
    plot(dati(:,i), target, '*','MarkerSize',7)
    hold on
    plot(mean(dati(:,i)),mean(target(:,1)),'r*','MarkerSize',7)
    hold on
    for j=1:size(dati,1)
        if(outliers(j,i)==1)
            % dove(k,1) = j;
            % dove(k,2) = i;
            % k = k+1;
            plot(dati(j,i),target(j),'*','Color',"#77AC30",'MarkerSize',7)
        end 
    end 

    xlabel([nomi{i}])
    ylabel([nom_target{1}])
    title([nomi{i}, 'vs' ,nom_target{1}])

end
legend('dati','media','outliers')
sgtitle('Scatter plot di tutte le variabili')

% pulizia
clear outliers

%% --------------------------------
% Istogrami di tutte le variabili

figure('Name','Istogramma')
for i = 1:size(dati,2)

    subplot(4,4,i)
    histogram(dati(:,i))
    title([nomi{i}])
    sgtitle('Istogramma di tutte le variabili')

end

%% --------------------------------
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
clear ind  


%% --------------------------------------------------------------------- %%
%     Modelli polinomiali
% --------------------------------
run modelliIMC.m
run modelliEta.m


%% --------------------------------------------------------------------- %%
%        Modello polinomiale di Grasso corporeo vs ICM vs Età
% -----------------------------------------------------------------------
run modello3.m


%% --------------------------------------------------------------------- %%
%       Stepwise Regression
% -----------------------------------------------------------------------
run stepwise_copy.m




