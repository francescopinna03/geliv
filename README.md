# GELIV, laboratorio Mata per strumenti deboli

Versione 0.1.3, 9 ottobre 2026. Codice di ricerca per Stata 16 o successivo.

Implementa ET, EL e CUE da zero in Mata per il modello IV lineare. Il duale usa
Newton smorzato e derivate analitiche. Il problema esterno usa `optimize()` con
BFGS, gradiente dell'inviluppo e più punti iniziali. Controlli e intercetta
restano nel modello e nei momenti. I momenti vengono costruiti solo a blocchi.

Le formule e gli algoritmi sono stati controllati numericamente in Python, con
il rapporto in `output/math_validation.json`. Dall'8 ottobre 2026 i moduli Mata
e il comando `geliv` sono compilati ed eseguiti in Stata/BE 19.0, dove superano
il self-test sintetico, il confronto del CUE sovraidentificato con `ivreg2` e i
gate di accettazione su Card, descritti in dettaglio in `docs/VALIDATION.md`.
Le stime complete su Angrist e Krueger e la simulazione calibrata restano da
eseguire, e il pacchetto non contiene tabelle empiriche che le anticipino.

## Avvio

Estrarre il pacchetto e impostarne la cartella come directory corrente in Stata.

```stata
cd "/percorso/geliv_toolkit"
do do/00_setup.do
do do/01_selftest.do
```

I benchmark richiedono pacchetti esterni, da installare una volta.

```stata
ssc install ivreg2
ssc install ranktest
ssc install avar
ssc install weakiv
ssc install bcuse
do do/10_card.do
do do/11_weakiv_benchmark.do
```

Per Angrist–Krueger, scaricare l'estratto ufficiale dal terminale e poi avviare
la preparazione in Stata. Il download calcola l'hash del file sorgente.

```bash
python3 python/download_data.py ak-minimal
```

```stata
do do/20_ak_prepare.do
do do/21_ak_estimate.do
do do/22_ak_weakiv.do
do do/40_geometry.do
do do/30_mc.do
```

`30_mc.do` contiene una configurazione pilota da 20 repliche. Modificarla prima
di una simulazione finale. Con 330 mila osservazioni e 240 momenti la stima
annidata può essere costosa. Non eseguire la simulazione grande prima dei test
sintetici e di Card. Non sono richieste operazioni di rete al motore Mata.
Su Stata/BE 19.0 una singola stima CUE con 3 strumenti esclusi ha richiesto circa
sei minuti, e il costo di ogni valutazione cresce con il quadrato del numero di
momenti, per cui le specificazioni con 180 strumenti vanno pianificate di
conseguenza.

## Interfaccia

```stata
geliv lwage, endog(educ) controls(exper expersq black smsa south smsa66 ///
    reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669) ///
    instruments(nearc4 nearc2) method(et) weights(p_et)
matrix list e(diagnostics)
matrix list e(starts)
```

`method()` accetta `et`, `el`, `cue`. `block(4096)` regola la memoria temporanea.
`noconstant` elimina l'intercetta. `from(M)` fornisce una matrice di punti iniziali,
uno per riga, nell'ordine endogene, controlli, intercetta. I coefficienti sono
espressi nelle unità originali. Sono accettate più endogene nel motore GEL.
Le variabili devono essere già espanse, senza sintassi factor-variable.
Non si eliminano silenziosamente strumenti collineari. Controllare il rango.

L'interfaccia Mata usa `gel_data_init(y,X,Z,gamma,block)` e `gel_fit(d,starts)`.
`X` include endogene, controlli e costante, mentre `Z` include strumenti
esclusi, controlli e costante. I tre valori di `gamma` sono 0, -1 e 1. Non è un motore
generico per qualsiasi funzione non lineare dei momenti.

Il risultato Mata contiene coefficienti, moltiplicatori nelle unità originali,
pesi impliciti, covarianza convenzionale, residui numerici e log dei punti
iniziali. `e(diagnostics)` contiene nell'ordine

| Campo | Significato |
|---|---|
| converged | convergenza esterna e certificazione numerica del duale |
| q | criterio medio, per ET uguale alla divergenza KL |
| Jgel | criterio GEL normalizzato, moltiplicato per 2n |
| Jdistance | 2n volte il criterio usato dall'ottimizzatore |
| inner_error | norma massima del gradiente interno, normalizzata anche per EL |
| outer_error | norma massima del gradiente esterno nelle coordinate scalate |
| moment_error | norma massima dei momenti pesati, strumenti scalati |
| ess | 1/somma(pesi²), solo per pesi non negativi |
| min_nw, max_nw | estremi del rapporto tra peso e 1/n |
| negative | numero di pesi negativi, possibile per CUE |
| N, m, p | osservazioni, momenti totali, coefficienti totali |

`e(starts)` registra codice di errore, flag di convergenza, criterio e norma del
gradiente per ogni partenza. Il minimo selezionato è il migliore tra le soluzioni
trovate, non un certificato di ottimalità globale. La matrice `e(V)` è la
covarianza HC0 di primo ordine sotto identificazione forte e momenti validi.
Non rende validi gli intervalli Wald sotto identificazione debole.

## Cosa viene testato

1. `01_selftest.do` verifica Newton CUE contro forma chiusa, gradienti esterni,
   invarianza rispetto al blocco, momenti pesati, identificazione esatta e
   inversione AR. Se disponibile, confronta CUE sovraidentificato con `ivreg2`.
2. `10_card.do` controlla N=3010 e il coefficiente 2SLS 0.131504 con i controlli
   dichiarati. ET, EL e CUE devono coincidere nel caso esattamente identificato.
   Il secondo strumento `nearc2` crea un test numerico sovraidentificato.
   Non implica una nuova difesa economica della sua validità.
3. CUE va confrontato con `ivreg2, cue robust`, non con il CUE omoschedastico.
   La tolleranza sul coefficiente di istruzione è 1e-6.
4. `11_weakiv_benchmark.do` confronta le statistiche canoniche omoschedastiche
   con `weakiv`, dopo averle riportate alla forma del moltiplicatore di Lagrange
   usata da quel comando, AR/(1 + AR/n) e analogamente per K. Il p-value CLR usa simulazione
   condizionale, quindi ha una tolleranza Monte Carlo, non sei cifre garantite.
   Se un assert fallisce, il gate resta aperto: investigare prima di procedere.

## Inferenza inclusa e confini

| Funzione | Oggetto | Ipotesi e limiti |
|---|---|---|
| `gel_full_tests(d,theta0)` | AR/S e K di tipo GMM, intero vettore strutturale | eteroschedasticità, osservazioni indipendenti, m fisso, rango regolare del Jacobiano corretto |
| `gel_null_grid(d,grid)` | criterio GEL del null, senza sottrarre il minimo | p=1, nessun nuisance, chi-quadro con m gradi di libertà |
| `gel_rf_test(r,beta,draws)` | AR, K e CLR canonici per una endogena | ramo omoschedastico; covarianza reduced-form stimata, CLR asintotico |
| `gel_ar_set(r,alpha)` | inversione AR F mediante disequazione quadratica | code infinite determinate analiticamente, errori gaussiani omoschedastici |
| `gel_rf_grid(r,grid,draws)` | inversione numerica K/CLR e AR | bordi della griglia non certificano code infinite |
| `weakiv` negli script | benchmark e inferenza scalare robusta con controlli | implementazione esterna, opzioni da mantenere coerenti |

Non viene inventato un CLR universale sostituendo ET o EL al criterio CUE.
Non è implementato un CLR eteroschedastico autonomo o cluster-robust in Mata.
Non è implementata l'inversione GEL profilata sui numerosi coefficienti dei
controlli. Per Card/AK l'inferenza robusta scalare con controlli è affidata al
benchmark `weakiv, robust` dopo una stima IV robusta. Il ramo autonomo canonical
serve per il confronto omoschedastico e per la simulazione calibrata gaussiana.

Un fallimento del duale è un esito numerico non risolto, non un rifiuto statistico
automatico. Il codice richiede radici accurate ma non contiene un certificatore
LP dell'appartenenza all'interno dell'inviluppo convesso.

## Scelte empiriche

Il confronto principale usa lo stesso campione e gli stessi controlli, inclusi
anno e Stato di nascita, con 3, 30 e 180 strumenti esclusi. Sono rispettivamente
3 trimestri, 30 interazioni trimestre–anno, e queste ultime più 150 interazioni
trimestre–Stato. Quarto trimestre e Michigan fungono da categorie base.
Con 59 controlli e intercetta, i momenti totali sono 63, 90 e 240.

Questa è un'estensione controllata del disegno classico. Non è una replica
letterale della Table V, né una riproduzione numerica dei risultati di Bound,
Jaeger e Baker. Aggiungere gli strumenti può cambiare anche la ponderazione
dell'effetto dell'istruzione se gli effetti sono eterogenei. La distanza tra
stime reali non identifica da sola il bias.

La simulazione fissa il disegno degli strumenti, calibra una direzione del primo
stadio e la covarianza degli errori sui dati, poi impone un coefficiente vero
esplicito e una concentrazione totale 10, 30 o 100. La stessa replica è usata
per tutte le specificazioni. Riporta bias, RMSE, errore mediano, copertura Wald,
copertura AR/K/CLR e fallimenti. Bias e RMSE sono condizionati al successo,
mentre la copertura include anche limiti inferiori e superiori che trattano i fallimenti.
La copertura dei test robusti è una proprietà della procedura di test del modello,
non sei diverse proprietà degli stimatori puntuali.

Leggere `docs/MATHEMATICS.md` per formule, normalizzazioni e interpretazione.
Leggere `docs/DATA.md` per fonti e per l'estratto completo.
`docs/VALIDATION.md` separa le verifiche eseguite da quelle ancora da eseguire
e riporta l'audit effettivamente eseguito sull'estratto MIT.
