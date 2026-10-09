# Dati e confrontabilità

I dati non sono incorporati nel pacchetto. Il downloader è esplicito e scrive
un manifest locale con URL, data del controllo, dimensione e hash SHA-256.

## Card

Fonte di distribuzione del dataset Wooldridge, Boston College.

https://fmwww.bc.edu/ec-p/data/wooldridge/card.dta

Alternativa in Stata, `bcuse card, clear`. Controlli usati nel benchmark
`exper expersq black smsa south smsa66 reg662-reg669`, oltre all'intercetta.
Campione atteso 3010, coefficiente 2SLS con `nearc4` circa 0.131504.
Il benchmark numerico da solo non valida l'identificazione economica.

## Angrist–Krueger, estratto minimo

Archivio originale di Joshua Angrist al MIT.

https://economics.mit.edu/sites/default/files/publications/asciiqob.zip

Le cinque colonne sono log salario settimanale, istruzione, anno di nascita,
trimestre di nascita, Stato di nascita. Lo schema è verificabile nel programma
originale https://economics.mit.edu/sites/default/files/inline-files/ak91.sas

La costruzione dei 180 strumenti è documentata da
https://economics.mit.edu/sites/default/files/publications/Jive180.sas

Il software trattiene i nati 1930–1939 e i codici FIPS di 50 Stati più DC.
Il campione effettivo deve essere registrato dal run, non imposto per far
coincidere un numero di osservazioni. Le matrici devono avere rango pieno.
Se alcuni gruppi sono assenti in un sottocampione, il software non deve
continuare a chiamare una matrice di rango inferiore «180 strumenti».

Le interazioni anno–trimestre sono costruite con Q1–Q3 come base, mentre
il programma Jive180 usa Q2–Q4. In presenza degli effetti fissi per anno i
due insiemi descrivono lo stesso spazio strumentale; la differenza di base
non cambia il modello. Per lo Stato è omesso Michigan, come nel file originale.
Nel disegno controllato anche i confronti 3 e 30 includono gli effetti di
Stato. Nel file originale questi entrano soltanto nella parte a 180 strumenti.

## Estratto completo

`python3 python/download_data.py ak-full` scarica
https://economics.mit.edu/sites/default/files/inline-files/NEW7080_1.rar

Estrarre il contenuto con un programma RAR e aprire il dataset Stata originale.
Il file include censimenti e coorti diversi. Non usarlo tutto come campione
1930–1939 del censimento 1980. Per la Table V usare le ricodifiche ufficiali
https://economics.mit.edu/sites/default/files/inline-files/QOB%20Table%20V.do

Il programma originale mappa v4 in EDUC, v9 in LWKLYWGE, v16 in CENSUS,
v18 in QOB e v27 in YOB, e contiene le ricodifiche dei controlli. Per un
confronto con covariate aggiuntive preparare esplicitamente un dataset comune,
aggiornare `GEL_W` e conservare il campione identico tra metodi. Il pacchetto
non indovina la colonna dello Stato di nascita nell'estratto completo.

## Memoria e tempi

Con 330 mila osservazioni e 240 momenti una matrice double piena pesa circa
634 MB decimali. I momenti non vengono materializzati per l'intero campione.
Il blocco 4096 per 240 pesa circa 7.9 MB. Restano in memoria le variabili
Stata, eventuali viste o copie Mata di X e Z e vettori di pesi. I prodotti
quadratici interni richiedono lavoro dell'ordine n*m² per iterazione Newton.
La gestione a blocchi limita la memoria, non elimina il costo aritmetico.
Il pacchetto non offre una promessa di tempo su Angrist–Krueger.
