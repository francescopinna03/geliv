# Definizioni e decisioni matematiche

## 1. Criterio e segni

Il vettore dei coefficienti strutturali, di dimensione p, è indicato con theta.
Il vettore degli strumenti inclusi ed esclusi, di dimensione m, è z_i.
Il vettore dei regressori, di dimensione p, è x_i. I momenti sono

    g_i(theta) = z_i (y_i - x_i' theta).

I moltiplicatori, di dimensione m, sono lambda. Con le funzioni concave
rho_ET(v)=-exp(v), rho_EL(v)=log(1-v), rho_CUE(v)=-(1+v)^2/2,
il criterio medio normalizzato è

    q_gamma(theta) = sup_lambda mean[rho_gamma(lambda'g_i)-rho_gamma(0)].

Si minimizza q rispetto a theta. EL richiede 1-lambda'g_i>0.
Il limite gamma=-1 della parametrizzazione Cressie–Read va inteso al netto
della costante divergente indipendente da parametri e moltiplicatori.

Il software minimizza il negativo dell'obiettivo interno concavo. Per ET usa
la trasformazione logaritmica monotona dell'esponenziale. A theta fissato,

    L(theta,lambda) = log(mean(exp(lambda'g_i(theta)))),
    q_ET_log(theta) = - min_lambda L(theta,lambda).

Il gradiente interno è la media dei momenti sotto i pesi esponenziali; l'hessiana
è la loro covarianza pesata. Questo permette Newton stabile con log-sum-exp.
Il criterio GEL originale è 1-exp(-q_ET_log). È una trasformazione strettamente
crescente, quindi i minimizzatori coincidono. Si salvano entrambi i criteri.

Per EL l'obiettivo interno convesso è -mean(log(1-lambda'g_i)); per CUE è
mean(lambda'g_i + (lambda'g_i)^2/2). Le rispettive hessiane sono positive
semidefinite. La convessità non garantisce esistenza di un ottimo finito.

La media campionaria dei momenti è gbar. La matrice dei momenti secondi non
centrati è Omega2=mean(g_i g_i'). Se invertibile,

    lambda_CUE = -solve(Omega2,gbar),
    q_CUE = (gbar' solve(Omega2,gbar))/2.

La covarianza centrata differisce da Omega2 fuori dal null. Non sono
intercambiabili per una pretesa uguaglianza numerica tra implementazioni.

## 2. Gradiente esterno

La derivata del residuo rispetto a theta è -x_i. Si definisce lo scalare
t_i=z_i'lambda. Le tre derivate dei criteri effettivamente minimizzati sono

    ET log     sum_i pi_i x_i t_i,
    EL         mean_i x_i t_i / (1-u_i t_i),
    CUE        mean_i (1+u_i t_i) x_i t_i.

Il codice ottimizza coefficienti scalati con la norma RMS dei regressori;
scala anche gli strumenti per le loro norme RMS. Le trasformazioni sono
invertibili, fisse rispetto a theta e non cambiano gli stimatori. I risultati
sono restituiti nelle unità originali. Nessuna FWL viene applicata a GEL.

La regolarizzazione di Newton, se necessaria per un passo, riguarda solo
l'hessiana usata per proporre una direzione. Non si aggiunge una penalità al
criterio e non si altera il sistema dei momenti. La convergenza richiede
comunque il gradiente del problema non penalizzato. I sistemi singolari
nel CUE chiuso o nella covarianza non vengono risolti con una pseudoinversa
silenziosa. Newton EL controlla anche il momento normalizzato, per evitare
una falsa convergenza quando tutti i pesi non normalizzati tendono a zero.

## 3. Probabilità e divergenze

La distribuzione empirica uniforme assegna 1/n a ogni osservazione. I pesi sono

    ET      pi_i = exp(v_i) / sum_j exp(v_j),
    EL      pi_i = [1/(1-v_i)] / sum_j [1/(1-v_j)],
    CUE     pi_i = (1+v_i) / sum_j (1+v_j).

Per ET ed EL, sotto fattibilità interna, sono probabilità. CUE non impone
1+v_i>=0. Può produrre pesi firmati. Non ha allora un'interpretazione come
distribuzione, né un'effettiva numerosità probabilistica dei pesi.

ET minimizza KL(pi||uniforme), EL minimizza KL(uniforme||pi), soggetti a
somma dei pesi uguale a uno e momenti pesati nulli. Il duale ET senza un
moltiplicatore esplicito per la normalizzazione produce la stessa soluzione
tramite la trasformazione logaritmica sopra. Al suo ottimo interno

    KL(pi_ET||uniforme) = -log(mean(exp(v_i))).

Per EL, all'ottimo interno, la somma dei pesi non normalizzati è n, e

    KL(uniforme||pi_EL) = mean(log(1-v_i)).

Il caso quadratico del GEL è equivalente al CUE con momenti secondi non
centrati. Un primale quadratico con vincolo di positività esplicito è diverso
quando il CUE non vincolato produce pesi negativi.

## 4. Identità pitagorica, con il suo dominio corretto

Fissare theta. Il politopo delle probabilità che soddisfano quei momenti è
una famiglia affine intersecata con il simplesso. Sia pi* la proiezione ET
dell'uniforme, interna al simplesso. Per ogni altra distribuzione fattibile r,

    KL(r||uniforme) = KL(r||pi*) + KL(pi*||uniforme).

Infatti log(pi*_i/(1/n))=lambda'g_i-log(mean exp(lambda'g_i)). La media del
primo termine sotto r è zero, come sotto pi*. L'identità segue sottraendo
le due divergenze. Non confronta automaticamente proiezioni calcolate a theta
diversi. La successiva minimizzazione in theta attraversa famiglie affini
diverse, la cui unione può non essere convessa.

In presenza di momenti annidati e a theta fissato si può costruire una
decomposizione informativa aggiungendo gruppi di vincoli. Con stime theta
diverse, chiamare tale differenza una decomposizione pitagorica esatta sarebbe
sbagliato. Analogamente, pesi estremi segnalano osservazioni necessarie alla
riconciliazione dei momenti, non invalidità individuale dello strumento.

## 5. Covarianza convenzionale e sovraidentificazione

Il Jacobiano medio, una matrice m per p, è D=-mean(z_i x_i'). Sotto
identificazione forte e corretta specificazione,

    V(theta_hat) = inverse(D' inverse(Omega2) D) / n.

È la covarianza efficiente di primo ordine condivisa dai tre stimatori.
Non incorpora correzioni higher-order né robustezza alla misspecificazione.
Sotto strumenti deboli i normali intervalli Wald possono avere copertura
errata anche quando questa matrice è invertibile.

Al parametro stimato e sotto regolarità con m fisso, il criterio GEL
moltiplicato per 2n ha riferimento chi-quadro(m-p). Per ET si restituisce
sia 2n*(1-exp(-KL)) sia 2n*KL. Sono asintoticamente equivalenti sotto il null,
non numericamente identici in un campione finito. Un test di distanza piccolo
non prova la validità economica degli strumenti.

## 6. Null strutturale e strumenti deboli

Per il null sull'intero vettore theta0 si valuta un sistema di m restrizioni
sulla media dei momenti. Il riferimento della statistica di tipo AR/S è
chi-quadro(m), senza richiedere identificazione forte di theta. Si mantiene
fisso m nell'approssimazione, con varianza dei momenti non degenere e CLT.
Si può usare il criterio GEL al null, senza sottrarre il suo minimo su theta.
Il rapporto profilato che sottrae il minimo non diventa automaticamente
chi-quadro(p) sotto strumenti deboli.

Per la statistica K robusta del modulo si indica con Omega la covarianza
centrata dei momenti al null. Per ogni coefficiente j, la covarianza tra
la colonna j del Jacobiano individuale e i momenti è una matrice m per m,
indicata con C_j. Il Jacobiano corretto ha colonne

    Dtilde_j = D_j - C_j inverse(Omega) gbar.

Si calcolano

    AR = n*gbar' inverse(Omega) gbar,
    K = n*gbar' inverse(Omega) Dtilde
          inverse(Dtilde' inverse(Omega) Dtilde)
          Dtilde' inverse(Omega) gbar.

Il riferimento K è chi-quadro(p), sotto le condizioni di rango e limite
congiunto appropriate. Non è il quadrato del semplice gradiente GEL con
l'hessiana osservata invertita. La correzione di covarianza è essenziale.
La funzione rifiuta la singolarità della matrice proiettata. Questi sono test
dell'intero vettore: inserire stime dei nuisance non conserva automaticamente
la giustificazione appena data.

## 7. Ramo canonico AR, K e CLR con una endogena

Per questo ramo soltanto si proiettano outcome, endogena e strumenti esclusi
fuori dai controlli esogeni W. Sia Y la matrice dei due outcome ridotti,
nell'ordine y e x. Sia Z la matrice degli strumenti esclusi residualizzati.
Le due matrici due per due di prodotti quadratici sono

    R = Y' Z inverse(Z'Z) Z' Y,
    E = Y'Y - R.

Il numero di strumenti esclusi è k, quello dei controlli inclusa l'intercetta
è l. Il modulo usa Sigma=E/(n-l-k). Al null beta, i vettori due per uno sono
b=(1,-beta)' e a=(beta,1)'. Le statistiche sufficienti scalari sono

    SS = b' R b / (b' Sigma b),
    TT = a' inverse(Sigma) R inverse(Sigma) a / (a' inverse(Sigma) a),
    ST = b' R inverse(Sigma) a / sqrt((b' Sigma b)*(a' inverse(Sigma) a)).

K=ST²/TT. Il criterio LR condizionale è

    CLR = (SS-TT + sqrt((SS-TT)²+4*ST²))/2.

La distribuzione condizionale è simulata tenendo TT fisso. Si estraggono due
variabili indipendenti U1~chi-quadro(1) e U2~chi-quadro(k-1), si sostituiscono
SS=U1+U2 e ST²=TT*U1. Il p-value registra anche l'errore Monte Carlo. Una stessa
matrice di estrazioni viene riutilizzata per tutta la griglia. Con covarianza
nota e normalità la costruzione canonica è esatta. Sostituendo la covarianza
con quella stimata dal reduced form, il CLR implementato è asintotico.

Il test classico AR usa SS/k e il riferimento F(k,n-l-k), esatto sotto errori
gaussiani omoschedastici. La sua inversione equivale a

    (1,-beta) [R - k*Fcrit*E/(n-l-k)] (1,-beta)' <= 0.

Il codice risolve la disequazione quadratica. Distingue insieme vuoto,
intervallo limitato, semiretta, due semirette e intera retta. Questo consente
di certificare code infinite senza scambiare il bordo di una griglia per
infinito. K e CLR rimangono griglie finite; la loro topologia tra nodi o fuori
dal dominio non è certificata. Un insieme AR vuoto può indicare incompatibilità
dei momenti, non una identificazione estremamente precisa.

## 8. Cosa afferma davvero Newey–Smith

Newey e Smith (2004) eliminano per la famiglia GEL una componente del bias
asintotico dovuta alla correlazione fra momenti e Jacobiano. EL elimina anche
una componente collegata alla stima della matrice dei pesi. Questo non dice
che ogni bias di ET, EL e CUE sia indipendente dal numero di momenti, né fornisce
una garanzia uniforme nell'asintotica weak-IV o many-IV. Il progetto deve
misurare prestazioni e fallimenti, non assumere la vittoria di GEL.

La simulazione è un esperimento sulle prestazioni finite, ispirato da queste
espansioni. Non è una verifica del teorema in un regime diverso dalle sue
ipotesi. Vanno distinti m fisso, m crescente, forza complessiva costante e
forza per strumento costante. Il programma implementa la concentrazione
complessiva fissata nel modello più ampio; la concentrazione proiettata sui
sottoinsiemi degli strumenti può essere inferiore.

## Fonti primarie

- Newey e Smith (2004), Econometrica 72, 219–255.
  https://doi.org/10.1111/j.1468-0262.2004.00482.x
- Kleibergen (2005), Econometrica 73, 1103–1123.
  https://doi.org/10.1111/j.1468-0262.2005.00610.x
- Moreira (2003), A Conditional Likelihood Ratio Test for Structural Models.
  https://doi.org/10.1111/1468-0262.00438
- Kitamura e Stutzer (1997), Econometrica 65, 861–874.
  An Information-Theoretic Alternative to Generalized Method of Moments Estimation.
- Finlay, Magnusson e Schaffer, codice e documentazione di weakiv.
  https://github.com/kfinlay/weakiv
- StataCorp, manuale Mata optimize().
  https://www.stata.com/manuals/m-5optimize.pdf
- Archivio originale Angrist.
  https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive
