# Stato della validazione

Verifiche in Python del 7 ottobre 2026, verifiche in Stata dell'8 e 9 ottobre 2026.

## Eseguito in Python

`python/verify_math.py` supera 35 controlli matematici e numerici, più 3009
confronti puntuali fra inversione analitica AR e test F. I controlli includono
gradienti e hessiane interni, gradienti dell'inviluppo, primali indipendenti
risolti con SLSQP, momenti pesati, CUE Newton contro forma chiusa, invarianza
rispetto a trasformazioni invertibili degli strumenti, uguaglianza in esatta
identificazione, identità pitagorica ET e accumulazione a blocchi. La simulazione
controlla inoltre il livello di K a concentrazione nulla e del CLR condizionale
a quattro valori di TT. Soglie ed errori sono nel JSON incluso.

`python/audit_ak_design.py` è stato eseguito sull'estratto minimo originale MIT.
Il checksum dell'archivio scaricato è

    0b6c505e52fd6c9f6194c356c91688acf2f8937b5023a84e754e942c9a2f2dcc

Le osservazioni sono 329509 prima e dopo i filtri dichiarati. Il campione è
comune alle tre specificazioni. I controlli sono 59 più intercetta.

| Strumenti esclusi | Momenti totali | Rango strumentale | Coefficiente 2SLS istruzione | F omoschedastico primo stadio |
|---:|---:|---:|---:|---:|
| 3 | 63 | 63 | 0.107693705 | 36.036354 |
| 30 | 90 | 90 | 0.086930200 | 5.363572 |
| 180 | 240 | 240 | 0.092818056 | 2.582341 |

Il coefficiente OLS con gli stessi controlli è 0.067338971. Queste sono stime
Python da prodotti incrociati sparsi, non output di Stata, Mata o GEL.
La traiettoria 2SLS non è monotona. Il dato serve anche a respingere una
narrativa empirica prestabilita. Non misura il bias senza conoscere il vero
coefficiente strutturale.

## Eseguito in Stata

L'8 ottobre 2026, su Stata/BE 19.0 per macOS, i due moduli Mata e il comando
`geliv` sono stati compilati senza errori. Il self-test sintetico di
`01_selftest.do`, su 800 osservazioni con due strumenti esclusi, supera tutti i
controlli, tra cui il confronto fra Newton e forma chiusa del CUE, la verifica
alle differenze finite dei gradienti esterni di ET, EL e CUE, l'invarianza
rispetto alla dimensione del blocco, la coincidenza con il 2SLS nel caso
esattamente identificato e la concordanza fra inversione analitica dell'AR e
griglia su 400 punti. Lo stesso file confronta il CUE sovraidentificato con
`ivreg2, cue robust`, e i due coefficienti differiscono di 1.1e-10.

Su Card (1995), con 3010 osservazioni, `10_card.do` riproduce il coefficiente
2SLS 0.131504 e verifica che ET, EL e CUE coincidano con esso nel caso
esattamente identificato, con pesi impliciti uniformi. Con `nearc4` e `nearc2`
come strumenti esclusi il CUE di `geliv` e quello di `ivreg2` coincidono su tutti
i 17 coefficienti entro 1e-6, con il coefficiente dell'istruzione pari a
0.1623756, mentre ET ed EL convergono a 0.1622988 e 0.1622307 con errore sui
momenti inferiore a 1e-15 e nessun peso negativo.

`11_weakiv_benchmark.do` confronta le statistiche canoniche omoschedastiche al
valore nullo 0.1 con quelle di `weakiv`. Il modulo calcola AR e K con la varianza
del residuo della forma ridotta completa, mentre `weakiv` adotta la forma del
moltiplicatore di Lagrange, costruita sul residuo dei soli controlli, e le due
versioni sono legate dall'identità AR_LM = AR/(1 + AR/n), che vale allo stesso
modo per K. Dopo la conversione AR, K e la statistica CLR coincidono con `weakiv`
alla settima cifra significativa, con valori 2.8329641, 1.4888246 e 1.6032545,
e il p-value del CLR differisce di 0.002, entro la tolleranza Monte Carlo della
simulazione condizionale.

Su Angrist e Krueger (1991), `20_ak_prepare.do` ricostruisce in Stata lo stesso
campione dell'audit Python, con 329509 osservazioni e insiemi di 3, 30 e 180
strumenti esclusi. Con 3 strumenti esclusi il CUE di `geliv` stima il
coefficiente dell'istruzione a 0.108965 e l'ET a 0.1089056, entrambi convergenti,
con statistica J pari a 3.10 su 2 gradi di libertà e pesi impliciti compresi fra
0.96 e 1.07 volte 1/n. La stima CUE ha richiesto 362 secondi su Stata/BE.

## Non ancora eseguito

Restano da eseguire le stime ET, EL e CUE con 30 e 180 strumenti esclusi, il
confronto sul campione completo previsto da `22_ak_weakiv.do`, l'analisi dei pesi
impliciti di `40_geometry.do` e la simulazione calibrata di `30_mc.do`. Al costo
osservato, le specificazioni con 180 strumenti richiedono molte ore per stimatore,
e per ridurre quel costo è prevista una variante del motore che aggrega i
prodotti esterni per cella di anno, trimestre e Stato di nascita, lasciando
invariati gli stimatori. I gate superati certificano il motore sui casi descritti
e non si estendono da soli ai risultati su Angrist e Krueger, che vanno letti solo
dopo l'esecuzione dei controlli previsti da ciascuno script.
