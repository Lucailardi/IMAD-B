% Quando si hanno molte variabili e si fa pulizia dei dati, bisogna stare
% attenti non solo agli outlier per la singola variabile, ma anche quelli
% che sono outlier rispetto alla correlazione fra due variabili (per
% esempio: il peso è 70 kg, l'altezza 1.90m, sono entrambi valori che
% singolarmente sono ragionevoli ma divergono dalla correlazione tipica.
% Una soluzione può essere la distanza di Mahalabis che descrive la
% distanza di un punto su uno scatterplot dal centro della nuvola.

% Imputazione: assegnare valore posticcio a una misura di un soggetto che è
% outlier solo in quella misura. Idea è che, soprattutto se ho pochi dati, 
% non voglio buttare via tutto il soggetto, ma sostituire il valore outlier
% con uno plausibile. 
% Operazione delicata perché sto falsificando dataset. Si può usare il 
% valore atteso di quella variabile condizionato dalle altre. La media va
% bene come imputazione solo se le variabili sono fra loro incorrelate, ma
% in generale il valore atteso condizionato il resto delle variabili
% performa meglio: se vi è un soggetto con peso di 200kg e alto 2m,
% sicuramente 200 non va bene, ma non possiamo assegnarli il valor medio di
% 80 kg.