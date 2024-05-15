%% ===================================================
%   RETI NEURALI
% ----------------
clear all
close all
clc

load dati_iden.mat
load dati_val.mat
load('dati_puliti.mat');

% -----------------------------------------------------------------------
% RBE posiziona e allena la rete:
% viene posizionato un neurone per ogni esempio del training set all'interno 
% dello strato nascosoto, vengono automaticamente quindi collocati i centri 
% tutti i nostri neuroni. --> non ho un ulteriore passagio di training 
% ma posso passare direttamente alla valutazione
            
% Vengono calcolati i pesi sinaptici anche tra lo strato della funzione di 
% base radiale e lo strato d'uscita della rete

% RELU È una funzione semplice e computazionalmente efficiente. 
% Risolve il problema del gradiente sparito per i valori positivi, ma 
% potrebbe causare "morte del neurone" per i valori negativi (i neuroni 
% che non attivano mai). Può essere più efficace in reti neurali profonde.
% -----------------------------------------------------------------------


% target = dati_iden(:, 1)';     % trainig set totale
% tSet = dati_iden(:,2:end)';
% testSet = dati_val(:,2:end)';  % test set
% y_target = dati_val(:,1)';
% 
% % training e validazione
% n_val = round(15*size(tSet,2)/100);   
% n_train = length(tSet) - n_val;
% 
% valSet = dati_iden([1:n_val],:);
% trSet = dati_iden([1+n_val:end],:);
% 
% varVal = valSet(:,2:end)';
% targVal = valSet(:,1)';
% 
% varTrain = trSet(:,2:end)';
% targTrain = trSet(:,1)';


%% ===================================================
%   prova disperata
% ---------------------
% In input do i risultati della regressione
% AddomeRVF + PolsoRVF + EtaCoscia + EtaCollo + PettoRVF + EtaBicipite + EtaGinocchio

nom_reg = {'AddomeRVF','PolsoRVF','EtaCoscia','EtaCollo','PettoRVF','EtaBicipite','EtaGinocchio'};
nom_par = {'Grasso'};

load('X.mat');  % set di identificazione
load('Xv.mat'); % set di validazione

X = cell2mat(X_);
Xvv = cell2mat(Xv);

data = [X(:,2:5); Xvv(:,2:5)];
grasso = dati_puliti(:,1);
data =[grasso, data];

clear grasso X Xvv X_ Xv

%% --------------------------------
%   CHATGPT
% ------------
% Carica il dataset 
% data = dati_puliti;
% clear dati_puliti

% Dividi il dataset in set di addestramento e set di test 
% (ad esempio, 80% per l'addestramento e 20% per il test)
train_ratio = 0.8;
num_samples = size(data, 1);
num_train = floor(train_ratio * num_samples);

train_data = data(1:num_train, :);
test_data = data(num_train+1:end, :);

% test set che non viene normalizzato
X_test = test_data(:, 2:end); 
y_test = test_data(:, 1);  

clear num_train num_samples

% -------------------
%   Normalizzazione
% -------------------
% proviamo con una normalizzazione del dataset - solo il training set
% Calcola il minimo e il massimo per ogni variabile
min_vals = min(train_data);
max_vals = max(train_data);

% Normalizza i dati utilizzando il Min-Max normalization
normalized_data = (train_data - min_vals) ./ (max_vals - min_vals);

% Estrai le variabili di input e l'output di grasso corporeo per l'addestramento
X_train = normalized_data(:, 2:end); 
y_train = normalized_data(:, 1);   

% --------------------------
%   Senza Normalizzazione
% --------------------------
% X_train = train_data(:, [end-1]); 
% y_train = train_data(:, 1); 


% -------------------
%       Rete
% -------------------

% set dei neuroni nei diversi strati della rete
neuroni = [20 20];
% spread è la sigma della gaussiana
spread = 0.05;

clear net
net = feedforwardnet(neuroni);
net.layers{1}.transferFcn = 'poslin';  % funzione relu like
net.trainParam.lr = 0.01;              % tasso di apprendimento
net.trainParam.goal = 0.05;             % mse desiderato
net.trainParam.BN = 'false';           % Normalizza i dati prima dell'addestramento
net.trainParam.BI = 'true';            % Inizializza i pesi in modo casuale

% net = fitnet(num_hidden_units);
% net = newrbe(X_train', y_train',spread);    % no addestramento
% net = newgrnn(X_train', y_train',spread);     % no addestramento

% Addestra la rete neurale
net = train(net, X_train', y_train');


% -----------------------------------
%   comportamento sul TRAINING set
% -----------------------------------
y_pred = net(X_train');
erryt = y_train' - round(y_pred);
errRatio = sum(abs(erryt))/numel(y_train);
clc
fprintf('La rete commette il %2.2f percento di errori sul TRAINING set \n', 100*errRatio)

mse = mean((y_train' - y_pred).^2); 
% MSE = sum((y_train' - y_pred).^2)/length(y_pred);
fprintf('Errore quadratico medio sul TRAINING set: %f\n', mse);
% fprintf('Errore quadratico medio sul TRAINING set: %f\n', MSE);


% -------------------------------
%   comportamento sul TEST set
% -------------------------------
y_predTEST = net(X_test');

% visualizzo gli errori sul test set: devo usare un arrotondamento
% delle uscite in modo che possano esere confrontate con i nostri
% target (se non metto round ho il 100% di errore)

errytTEST = y_test' - round(y_predTEST);

% prestazioni della rete
errRatioTEST = sum(abs(errytTEST))/numel(y_test);
fprintf('La rete commette il %2.2f percento di errori sul TEST set\n', 100*errRatioTEST)

% Calcola l'errore quadratico medio
mseTEST = mean((y_test' - y_predTEST).^2); 
% MSETEST = sum((y_test' - y_predTEST).^2)/length(y_predTEST);

% Visualizza i risultati
fprintf('Errore quadratico medio sul TEST set: %f\n', mseTEST);
% fprintf('Errore quadratico medio sul TEST set: %f\n', MSETEST);

% ------------------------------------
%            grafici 
% ------------------------------------

figure('Name','Visualizzazione')
subplot 121
plot(y_pred')
hold on
plot(y_train)
legend('output rete','y vera')
title('TRAINING set')

subplot 122
plot(y_predTEST')
hold on
plot(y_test)
legend('output rete','y vera')
title('TEST set')


% goodness of fit 
figure('Name','Goodness of fit')
subplot 121
x = 0:0.01:length(y_pred);
scatter(y_pred',y_train)    
hold on
plot(x,x*0.5,'-r')
axis equal
xlabel('output rete')
ylabel('valore vero')
title('TRAINING set')

subplot 122
xT = 0:0.01:length(y_predTEST);
scatter(y_predTEST',y_test)    
hold on
plot(xT,xT*0.5,'-r')
axis equal
xlabel('output rete')
ylabel('valore vero')
title('TEST set')

