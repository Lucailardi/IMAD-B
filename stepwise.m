%% -------------------------------------------------------------------- %%
%                               REGRESSIONE
% --------------------------------------------------------------------- %%
clear 
close all
clc

% caricamento dei dati già puliti da outliers e già divisi in set di
% identificazione e di validazione
load("dati_puliti.mat")
load("dati_iden.mat")
load("dati_val.mat")

alpha = 0.05;

% considerando solo i dati di identificazione, divido i regressori (covariate)
% dal parametro target che è il grasso corporeo
regressori = dati_iden(:,2:end); % covariate (regressori=x)
parametri = dati_iden(:,1);        % target

% set di validazione
regres = dati_val(:,[2:end]); % covariate (regressori=x)
Y_val = dati_val(:,1);        % target


nom_reg = {'Eta','Peso','Altezza','Collo','Petto','Addome','Anca','Coscia',...
    'Ginocchio','Caviglia','Bicipite','Avambraccio','Polso','BMI','RVF'};
nom_par = {'Grasso'};


%% --------------------------------------------------------------------- %%
%     M --> matrice dei regressori candidati - identificazione
%    VAL -> matrice dei regressori candidati - validazione 
% ---------------------------------------------------------------------- %
% non abbiamo variabili dummy 

% candidati singoli
for i = 1:size(regressori, 2)
    M{i} = regressori(:,i);
    reg{i} = nom_reg{i};
end

% candidati del secondo ordine
clear i j
for i = 1:size(regressori, 2)
    for j = 1:size(regressori, 2)
        M{length(M) + 1} = M{i}.*M{j};
        reg{length(M)} = strcat(nom_reg{i},nom_reg{j});
    end 
end

% mi salvo i nomi di tutti i regressori
tutti = reg;

% ------------------------------------------------
% faccio la stessa cosa col set di validazione
% ------------------------------------------------
% candidati singoli
clear i j
for i = 1:size(regres, 2)
    VAL{i} = regres(:,i);
end

% candidati del secondo ordine
clear i j
for i = 1:size(regres, 2)
    for j = 1:size(regres, 2)
        VAL{length(VAL) + 1} = VAL{i}.*VAL{j};
    end 
end


%% -----------------------------------------------------------------------%
% con la funzione matlab 'stepwise'
% ----------------------------------------------

% richiede le variabili in forma matriciale enon di cell array
MM = [M{:}]; 

% [b,se,pval,finalmodel,stats] = stepwisefit(MM,parametri);
[b,~,~,finalmodel,~] = stepwisefit(MM,parametri);
index = find(finalmodel==1);
beta = b(index);

disp('---------------------')
disp('modello stepwisefit')
for i=1:length(index)
    disp([tutti(index(i))])
end 

% estraggo solo i regressori che fanno parte del modello
Xsw = cell(1,length(index));    % identificazione
Xv_sw = cell(1,length(index));  % validazione

for i =1:length(index)
    Xsw{1,i} = M{1,index(i)};       % identificazione
    Xv_sw{1,i} = VAL{1,index(i)};   % validazione
end 

% applico la regressione 
Ysw = [Xsw{:}]*beta;
Yv_sw = [Xv_sw{:}]*beta;



%% --------------------------------------------------------------------- %%
%                Stepwise 'a mano'
% -------------------------------------------------
% Trovare migliore modello di regressione
p = length(reg);

% inizializzazione
Y = parametri;
n = size(Y,1);
c = 1;
X_ = {};
X_{1,1} = ones(n,1);
C = M;      % matrice dei candidati
n_reg = {};
R2a = [];
Ftest = [];


while c <= p
    m = size(C,2);
    for i = 1:m
        Xp = [[X_{:}] C{i}];
        B_hat = pinv(Xp'*Xp)*Xp'*Y;
        Y_hat = Xp*B_hat;

        res = Y-Y_hat;
        SSR = sum((Y_hat-mean(Y)).^2);
        SSE = sum(res.^2);
        SST = sum((Y-mean(Y)).^2);

        k = (size(Xp,2)-1);     % numero di regressori
        n = (size(Xp,1));       % numero di osservazioni
        R2 = SSR/SST;
        R2av(i) = R2 - (k*(1-R2)/(n-k-1));

        Ftest(i)=(SSR/k)/(SSE/(n-k-1));
    end

    [Fmax, indF] = max(Ftest);
    R2a = R2av(indF);
    Chi{1,c} = reg{indF};
     if c == 1
        c = 2;
        best_reg = C{indF};
        n_reg = {n_reg{:} reg{indF}};
        R2ac = R2a;
        C = C(:, 1:length(C)~=indF);
        reg = reg(:, 1:length(reg)~=indF);
        X_{1,c} = best_reg(:,:);
    else
        if R2a > R2ac
            best_reg = C{indF};
            n_reg = {n_reg{:} reg{indF}};
            R2ac = R2a;
            C = C(:, 1:length(C)~=indF);
            reg = reg(:, 1:length(reg)~=indF);
            % crea un vettore di valori booleani che rappresenta quali colonne 
            % devono essere mantenute e quali eliminate
            % ~=indF restituisce 'true' per tutti gli indici che NON corrispondono a indF 
            % e 'false' per l'indice indF.
            c = c+1;
            X_{1,c} = best_reg(:,:);
        else
            break
        end
     end
     R2av = [];
     Ftest = [];
end

% coefname = n_reg;
coefname = Chi;
Beta = inv([X_{:}]'*[X_{:}])*[X_{:}]'*Y;

B_hat = [Beta];      % coefficienti beta
Y_hat = [X_{:}]*B_hat;

% ================================================================== %
res = Y-Y_hat;
SSR = sum((Y_hat-mean(Y)).^2);
SSE = sum(res.^2);
SST = sum((Y-mean(Y)).^2);

k = size(X_,2)-1;     % numero di regressori
n = size(Y,1);        % numero di osservazioni

se = SSE/(n-k-1);     % stima della varianza degli errori

mat = inv([X_{:}]'*[X_{:}]);

for l = 1:length(B_hat)
    s=sqrt(mat(l,l)*se);    % varianza dei beta            
    CV(l)=abs(s/B_hat(l));  % CV=stddev/beta
    t=B_hat(l)/s;           % statistica test t
    
    t_crit=tinv(alpha/2, n-k-1);
    
    % intervalli di confidenza    
    I_inf(l)=(B_hat(l)'+t_crit*s)'; 
    I_sup(l)=(B_hat(l)'-t_crit*s)';

end
% ================================================================== %

disp(['Migliori regressori a mano: '])
for i = 1:size(coefname,2)
    if i < size(coefname,2)
        disp([coefname{1,i}, ' + '])
    else 
        disp([coefname{1,i}])
    end 
end 
disp('--------------------------')
disp('Intervalli di confidenza: ')
for g = 1:length(B_hat)
    disp(['Beta',num2str(g-1), ': ', num2str(B_hat(g)), ' ∈ [', num2str(I_inf(g)), '; ', num2str(I_sup(g)),']'])
end
disp('--------------------------')


%% ---------------------------------------------------------------------- %
% Verificare la significativit`a statistica dei 
% regressori inclusi nel modello.

% Calcolo statistica F
F = (SSR/k)/(SSE/(n-k-1));
% trovo il p_value
p_value = 1 - fcdf(F, k, n-k-1);
alfa=0.05;

if p_value<=alfa
    disp(['Rifiuto H0 ->Esiste almeno un regressore utile...' ...
        'per predirre la variabile dipendente ']);
else
    disp(['Non Rifiuto H0 ->Nessuno dei regressori è decisamente...' ...
        ' utle per predirre la variabile dipendente']);
end
disp('==========================')


%% --------------------------------------------------------------------- %%
%      'a mano' sul set di validazione
% ----------------------------------------
clear indice

% trovo gli indici dei migliori regressori
for i=1:length(tutti)
    for j=1:length(Chi)

        if (strcmp(Chi{1,j},tutti{1,i})==true)
            indice(j) = i;
        end 
    end 
    
end 

% estraggo le colonne giuste dal set di validazione
clear i j
Xv{1,1} = ones(length(VAL{1,1}),1);
for k = 1:length(indice)-1
    Xv{1,k+1} = VAL{1,indice(k)};
end 

% applico la regressione a mano al set di validazione
Y_hat_val = [Xv{:}]*B_hat;
xx = linspace(0,100,length(Y_hat_val));  % serve per il plot


%% --------------------------------------------------------------------- %%
% come si comporta la predizione sul set di validazione?
% spoiler: molto male

x = linspace(0,100,length(Y_hat));

figure('Name','Regressione')
subplot 211
plot(x, Y, 'xr')
hold on
plot(x, Y_hat, 'bo')
hold on 
plot(x, Ysw, 'm^')
title('Set di identificazione')
legend('identificazione','regressione a mano','regressione con funzione')

subplot 212
plot(xx, Y_val, 'xr')
hold on
plot(xx, Y_hat_val, 'bo')
hold on
plot(xx, Yv_sw, 'm^')
title('Set di validazione')
legend('identificazione','regressione a mano','regressione con funzione')


%% ===================================================================== %%
% Predizione miei parametri

% nom_reg = {'Eta','Peso','Altezza','Collo','Petto','Addome','Anca','Coscia',...
%    'Ginocchio','Caviglia','Bicipite','Avambraccio','Polso','BMI','RVF'};

% %riguarda
% mydata = [23, 50, 162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% 
% grasso = 0;
% for i =1:length(mydata)
%     grasso = grasso + Beta(i)*mydata(i);
% end 
