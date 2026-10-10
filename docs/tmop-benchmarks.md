# Benchmark numerici pubblicati per TMOP

Nota per scegliere test riproducibili del Target-Matrix Optimization Paradigm (TMOP) da reimplementare nel confronto di GCEE con altri stimatori a priori usati dentro la r-adattività. I nomi dei metodi, i titoli degli articoli e gli identificatori di codice restano in inglese.

Le fonti sono i testi aperti (arXiv, rapporti OSTI, PDF IMR, header del miniapp MFEM su `master` consultato nell’ottobre 2026). Dove un dettaglio non è nel testo letto, è marcato come incerto. Non ho potuto aprire `paperA.tex`.

## Come sono costruiti i target

In TMOP, in ogni punto di campionamento, \(A\) è lo Jacobiano riferimento→fisico e \(W\) è lo Jacobiano bersaglio. Si minimizza un funzionale di \(\mu(T)\) con \(T = A W^{-1}\), muovendo i nodi e lasciando invariata la topologia. La fattorizzazione usata nei lavori successivi (Dobrev et al., *hr-Adaptivity…*, §2.2) scrive

\[
W = \zeta\, R\, Q\, D
\]

con \(\zeta\) size, \(R\) rotation, \(Q\) skew, \(D\) aspect ratio.

Due famiglie di test, da non mescolare nel confronto con GCEE:

- **Geometrici.** \(W\) è ideale (quadrato, cubo, simplesso equilatero), a size uniforme, oppure copia la size iniziale. Stress: forma, simmetria radiale, boundary layer, nodi poco spostati, untangling, ordine alto, bordi curvi.
- **r-adattività.** \(W\) dipende da un campo: funzione analitica, volume fraction, gradiente di una soluzione nota, oppure una misura di errore. È qui che uno stimatore a priori entra, di solito come size (e, se anisotropo, come aspect ratio o allineamento).

Nessun articolo TMOP che ho letto confronta due stimatori a priori tra loro. Il punto più vicino è *PDE-Constrained High-Order Mesh Optimization* (2026), che mette tre misure di errore dentro lo stesso obiettivo insieme alla qualità TMOP.

## Letteratura canonica

| Lavoro | Ruolo | URL |
| --- | --- | --- |
| Knupp, *Introducing the target-matrix paradigm for mesh optimization via node-movement*, Eng. Comput. 28:419–429, 2012. DOI [10.1007/s00366-011-0230-1](https://doi.org/10.1007/s00366-011-0230-1) | Formulazione del paradigma. Implementazione in Mesquite. | paywall Springer; scheda Sandia |
| Dobrev, Knupp, Kolev, Mittal, Tomov, *The Target-Matrix Optimization Paradigm for High-Order Meshes*, SIAM J. Sci. Comput. 41(1):B50–B68, 2019. DOI [10.1137/18M1167206](https://doi.org/10.1137/18M1167206). arXiv [1807.09807](https://arxiv.org/abs/1807.09807), rapporto [OSTI](https://www.osti.gov/servlets/purl/1729749) | Estensione high-order, metriche, Newton, test geometrici. Codice in MFEM. | aperto (arXiv / OSTI) |
| Dobrev, Knupp, Kolev, Tomov, *Towards Simulation-Driven Optimization of High-Order Meshes by the Target-Matrix Optimization Paradigm*, 27th International Meshing Roundtable, Springer, 2019, pp. 285–302. PDF [IMR27](https://project.inria.fr/imr27/files/2018/09/1034.pdf) | Primi target da un campo discreto \(\eta\). | PDF IMR aperto; il DOI Springer non è stato verificato |
| Dobrev, Knupp, Kolev, Mittal, Rieben, Tomov, *Simulation-Driven Optimization of High-Order Meshes in ALE Hydrodynamics*, Computers & Fluids 208:104602, 2020. DOI [10.1016/j.compfluid.2020.104602](https://doi.org/10.1016/j.compfluid.2020.104602). arXiv [2001.11536](https://arxiv.org/abs/2001.11536) | r-adattività ALE in BLAST, trigger, limiting. | arXiv aperto; l’applicazione BLAST è chiusa |
| Dobrev, Knupp, Kolev, Mittal, Tomov, *hr-Adaptivity for nonconforming high-order meshes with the target matrix optimization paradigm*, Eng. Comput. 38(4):3721–3737, 2022. DOI [10.1007/s00366-021-01407-6](https://doi.org/10.1007/s00366-021-01407-6). arXiv [2010.02166](https://arxiv.org/abs/2010.02166) | hr-adattività; unico test con soluzione esatta e norma dell’energia. | arXiv aperto; algoritmo in MFEM |
| Kolev, Lazarov, Mittal, Schmidt, Tomov, *PDE-Constrained High-Order Mesh Optimization*, Eng. Comput. 42(1), 2026. DOI [10.1007/s00366-025-02250-9](https://doi.org/10.1007/s00366-025-02250-9). arXiv [2507.01917](https://arxiv.org/abs/2507.01917) | r-adattività vincolata alla PDE; confronto di tre misure di errore. | arXiv aperto; implementazione dichiarata in MFEM, non come miniapp documentato |

Articoli vicini, non usati qui come benchmark di stimatori: Knupp, Kolev, Mittal, Tomov, *Adaptive Surface Fitting and Tangential Relaxation for High-Order Mesh Optimization*, arXiv [2105.12165](https://arxiv.org/abs/2105.12165) (fitting di interfacce, Taylor–Green, Rayleigh–Taylor); Mittal, Dobrev, Kolev, Tomov, *High-Order Mesh r-Adaptivity with Tangential Relaxation and Guaranteed Mesh Validity*, arXiv [2601.17708](https://arxiv.org/abs/2601.17708) (rivisita pala e untangling).

### Knupp 2012 e i test precedenti

Il testo del 2012 è a pagamento e non l’ho letto per intero. L’abstract elenca classi di problemi, non mesh con nome: shape improvement, boundary-layer preservation con rimozione di sliver, high-order mesh improvement, edge-length equalization. L’implementazione citata è la libreria Mesquite.

I casi con nome che ho potuto leggere stanno in due presentazioni, non nel PDF della rivista:

- Formulazione del 20 febbraio 2007 ([OSTI 1143308](https://www.osti.gov/servlets/purl/1143308)), esperimenti MATLAB: mesh replication (\(T=I\), \(sI\), \(R\), \(sR\)), replication su domini scalati e molto diversi, barrier metrics, horseshoe, hybrid mesh, effetto dell’esponente \(p\) della power mean, pesi \(c_n\).
- Talk IMR 2010 ([OSTI 1121788](https://www.osti.gov/servlets/purl/1121788)), già in Mesquite: shape optimization su una sfera; elementi finiti quadratici (untangle con metrica d’area, poi shape-barrier, \(W\) equilatero costante, sei punti campione); custom smoother size+shape con \(W = LVS\), \(S=I\), size = lunghezza media del lato; sliver removal su una mesh tetraedrica viscosa CFD, con coefficiente di trade-off a soglia sul diedro.

Non ho verificato quali di questi casi siano ristampati, con le stesse figure, nell’articolo del 2012. Per un confronto con GCEE sono poco utili: non c’è una soluzione PDE, e le mesh del talk non risultano distribuite insieme al testo.

## Test geometrici (SISC 2019, §5)

Tutti in MFEM, funzionale variazionale, quadratura di Gauss–Lobatto, Newton a tolleranza relativa \(10^{-12}\), MINRES. I nodi sui bordi curvi sono fissi; gli altri possono scorrere se non escono dal dominio. Sono test di forma, non di stimatore.

**Turbine blade.** 2D, mesh \(Q_4\) `blade.mesh`. Target isotropo quadrato \(W_1\), metrica di forma \(\mu_2\), 6 punti di quadratura per direzione. Stress: curvatura interna (i lati interni partono dritti), vertici di valenza 5, rilassamento tangenziale sui lati allineati agli assi. Obiettivo da \(170.792\) a \(66.4977\) (circa \(-61\%\)); massimo di \(\mu_2\) da \(21.5\) a \(16.7\). Ripetizione con spostamento limitato, \(\delta=5000\) nel funzionale (4.2): obiettivo a \(126.28\) (circa \(-26\%\)), massimo di \(\mu_2\) a \(6.12\). Sezione 5.1. Citazione: Dobrev et al. 2019, [arXiv:1807.09807](https://arxiv.org/abs/1807.09807).

**ICF, shape + equal size.** 2D, mesh \(Q_3\) perturbata di tipo inertial-confinement-fusion (il numero di elementi non è stampato nel testo). \(W_2\) a size uguale, metrica shape+size \(\mu_9\), 6 punti per direzione. Stress: simmetria radiale e size uniforme. Obiettivo da \(1.119\) a \(0.329\) (circa \(-70\%\)). Sezione 5.2.1.

**ICF, shape + initial size.** Stessa mesh. \(W_3 = \sqrt{\det A(x_0)}\, I\), ancora \(\mu_9\). Stress: migliorare la forma tenendo la size locale iniziale. Obiettivo da \(0.095\) a \(0.046\) (circa \(-51\%\)). Sezione 5.2.2.

**ICF, forma con nodi poco spostati.** \(\mu_1\) e \(W_1\). Con \(\delta=0\) l’obiettivo passa da \(0.189\) a \(0.141\) (circa \(-25\%\)) ma gli elementi si spostano e si ridimensionano. Con \(\delta=10\) la riduzione è circa \(2.76\%\) (\(0.189\) → \(0.184\)) e lo spostamento resta locale. Sezione 5.2.3. Utile come controllo “non voglio una r-adattività aggressiva”, inutile come banco per uno stimatore.

**ICF, anelli di size.** Combinazione \(\mu_1\) (forma, \(W_1\)) più peso spaziale \(\eta(r)\) su \(\mu_{77}\) (size), con target di size più piccola della media. \(\eta\) è una somma di tangenti iperboliche su due corone, \(r\in[0.16,0.17]\) e \([0.23,0.24]\). Obiettivo da \(7.291\) a \(1.421\) (circa \(-80\%\)). Sezione 5.2.4. È un target di size analitico, ancora geometrico: \(\eta\) non viene da un errore.

**Pinched sphere.** 3D. Metrica di forma \(\mu_{303}\), target \(W_1\). Mesh \(Q_4\) con 448 elementi, 6 punti per direzione: obiettivo da \(169.605\) a \(29.263\) (circa \(-82\%\)). Mesh \(Q_2\) con 229376 elementi, 4 punti per direzione. Una \(Q_2\) da 1.835 milioni di elementi è stata eseguita su 512 task MPI; il risultato non è mostrato. Sezione 5.3. Stress: forma interna 3D, ordine alto, bordo curvo fissato.

## Test di r-adattività guidata da un campo

### Banda sinusoidale (IMR 2019, §6.1)

2D, mesh cartesiana \(16\times 16\) di ordine 3. Il campo discreto \(\eta\in[0,1]\) è l’interpolante di

\[
\tanh\bigl(10(y_0-0.5)+\sin(4\pi x_0)+1\bigr)
- \tanh\bigl(10(y_0-0.5)+\sin(4\pi x_0)-1\bigr).
\]

Il testo letto non scrive l’intervallo del dominio; la funzione è centrata su \(y=0.5\). Gauss–Legendre, Newton a \(10^{-12}\), derivate di \(W\) poste a zero, \(\eta\) trasportato per advection remap. Tre costruzioni, stesso scopo (area circa tre volte più piccola dove \(\eta\approx 1\)):

| Variante | Target e metrica | Obiettivo |
| --- | --- | --- |
| Size | \(W_1\) isotropo con rapporto \(\alpha=7\), shape+size \(\mu_7\) | \(2.165\to 1.136\) (circa \(-47\%\)) |
| Shape | \(W_2\) diagonale con aspect ratio \(\beta=3\), shape \(\mu_{58}\) | \(1545.92\to 1207.33\) (circa \(-22\%\)) |
| Shape+size | \(W_3\) con \(\gamma=3\) e size \(s\), shape+size \(\mu_9\) | \(0.756\to 0.258\) (circa \(-65\%\)) |

È il test classico di r-adattività TMOP su un campo liscio. Non misura l’errore di una PDE. Citazione: Dobrev et al., IMR 2019, [PDF](https://project.inria.fr/imr27/files/2018/09/1034.pdf), §6.1.

### Gas impact, versione IMR (§6.2)

2D, dominio \([0,2]\times[0,2]\), mesh \(60\times 60\) di ordine 2, \(t=10\). Tre materiali (impactor \(0.15\le x\le 0.65\), \(0.9\le y\le 1.1\); wall \(0.1\le x\le 1.0\); resto background). Target isotropo da \(\eta=\max(\eta_{\mathrm{wall}},\eta_{\mathrm{imp}})\), \(\alpha=10\), metrica \(\mu_9\). Eseguito in BLAST. Stress: size che insegue interfacce mobili. Nessun errore di stimatore. Sorgente: Barlow, Hill, Shashkov, J. Comput. Phys. 276:92–135, 2014, per il problema idrodinamico.

### Triple Point, gas impact e ball impact (Computers & Fluids 2020, §5)

Stesso quadro, applicazione BLAST (chiusa). L’abstract cita “shock positions, material regions, known error estimates” come ingressi possibili. Nei test numerici l’ingresso è un indicatore di materiale o di interfaccia, non uno stimatore.

**Triple Point.** 2D, tre materiali, problema di Dobrev, Kolev, Rieben, SIAM J. Sci. Comput. 34:606–641, 2012. Il testo ALE non ristampa il dominio. Simulazione lagrangiana fino a \(t=3.5\) presa come riferimento (niente diffusione di interfaccia). Errore = differenza media della posizione di un materiale su 10000 punti equispaziati. Discretizzazione \(Q_3Q_2\). Adattività: shape+size \(\mu_7\), size piccola sulle interfacce, \(\alpha=10\). Trigger automatico con \(S_{22}=4\). Tabella 2: errore 0.033 (Eulerian), 0.022 (limiting \(a=1/2\)), 0.020 (\(a=1/3\)), 0.017 (adapted, auto trigger), 0.022 (adapted ogni 50 passi). Sezione 5.1.

**Gas impact, versione 2020.** Stesso dominio e materiali della versione IMR, ma mesh \(80\times 80\) di ordine 2 (non \(60\times 60\)), metrica \(\mu_7\) (non \(\mu_9\)), \(\eta\) che include anche un indicatore di interfaccia ricostruito, trigger \(S_{22}=4\). 49400 passi lagrangiani, 38 passi ALE, 760 passi di advection. Sezione 5.2. Le due versioni non sono lo stesso numero.

**Steel ball impact.** 2D assialsimmetrico e 3D, mesh non strutturate high-order. Sfera d’acciaio, raggio 5.5, velocità iniziale \(0.31\,\mathrm{cm}/\mu\mathrm{s}\); piastra di alluminio, raggio 24, spessore 2.5. Tempo finale \(t=80\). Size verso impactor, piastra e interfacce, \(\mu_7\). In 3D \(\alpha=12\) e \(S_{33}=8\). Sezione 5.3. Problema originale di Howell e Ball, J. Comput. Phys. 175:128–167, 2002. Chiuso (BLAST) e lontano da uno stimatore a priori ellittico.

### Poisson a fronte circolare (hr-adattività, §4.1)

2D, \(\Omega=[0,1]^2\). L’equazione scritta è \(\nabla^2 u = f\), con \(f\) scelto in modo che \(u\) sia nota. Il paper la chiama “the 2D benchmark from” Cerveny, Dobrev, Kolev, *Nonconforming mesh refinement for high-order finite elements*, SIAM J. Sci. Comput. 41(4):C367–C392, 2019. Non ho riletto Cerveny, quindi non confermo che la formula sia già lì.

Formula nel PDF arXiv, spezzata su più righe e incompleta in estrazione automatica. Il testo dice fronte circolare di raggio \(r=0.7\) e \(\alpha=200\), mesh di figura \(16\times 16\) di ordine 2. Le slide ICOSAHOM 2021 degli stessi autori ([PDF LLNL](https://computing.llnl.gov/sites/default/files/inline-files/icosahom2021_0.pdf)) danno

\[
u=\arctan\Bigl[\alpha\Bigl(\sqrt{(x-x_c)^2+(y-y_c)^2}-r\Bigr)\Bigr],
\]

con \(\alpha=200\), \((x_c,y_c)=(-0.05,-0.05)\), \(r=0.7\). La radice quadrata è coerente con la frase “radius \(r\)” e con la formula, questa sì esplicita, del paper 2026 (sotto), che usa \(\beta=20\) invece di \(\alpha=200\). Prima di implementare, conviene ricontrollare l’equazione (13) sul PDF, non su un estratto.

Target (15): shape+size, skew ideale \(\phi=\pi/2\), size \(\zeta\propto\|\nabla u\|\), aspect ratio \(\rho\propto \partial_x u/\partial_y u\), metrica \(\mu_9\). Il gradiente è quello della **soluzione esatta**, non di uno stimatore. Studio quantitativo: si parte da \(4\times 4\) e si raffina; un solo giro di \(h\) e di \(r\). Errore in norma dell’energia rispetto a \(u\). Il testo dice circa il \(66\%\) di gradi di libertà in meno con hr rispetto alla sola r, a parità di accuratezza. Sezione 4.1, figure 7–9.

Il §3.6 dice che lo stimatore TMOP usato per decidere il raffinamento \(h\) “can also be replaced or used in conjunction with other (error) estimators”. Tutti gli esperimenti del paper usano lo stimatore TMOP. Non c’è un confronto fra stimatori.

### Anello analitico, simplessi, derefinement (hr, §4.2–4.4)

Dominio \([0,1]^2\). Per la sola size:

\[
\eta=\tanh(\beta(r-0.2))-\tanh(\beta(r-0.3)),\quad
\zeta=\eta\delta+(1-\eta)\psi\delta,
\]

con \(\beta=30\), \(\delta=0.001\), \(\psi=10\), aspect ratio 1, skew \(\pi/2\). Metriche: \(\mu^h_{55}\) (size) e \(\mu^r_7\) (shape+size). Su \(8\times 8\) (\(N_E=64\)): \(h\) riduce \(F\) del \(40.36\%\) (\(N_E=484\)), \(r\) non riduce \(F\), hr del \(69.2\%\) (\(N_E=616\)). Su \(16\times 16\): \(h\) \(21.9\%\) (\(N_E=544\)), \(r\) \(51.8\%\), hr \(67.3\%\) (\(N_E=616\)).

Variante shape+size con \(\rho=f(r,\theta)\), metrica \(\mu^h_7\). Su \(8\times 8\): \(h\) \(92.7\%\) (\(N_E=572\)), \(r\) \(6.1\%\), hr \(93.5\%\) (\(N_E=704\)). Su \(16\times 16\): \(h\) \(65.6\%\) (\(N_E=600\)), \(r\) \(60.6\%\), hr \(85.3\%\) (\(N_E=716\)).

Simplessi, \(N_E=128\), \(\mu^r_9\) e \(\mu^h_{55}\), stessi target di size: \(h\) \(62.4\%\) (\(N_E=928\)), \(r\) \(43.9\%\), hr \(85.2\%\) (\(N_E=1100\)).

Derefinement: \(4\times 4\) raffinata quattro volte, \(N_E=4096\). La \(r\) riduce \(F\) del \(55.4\%\) senza togliere elementi; hr scende a \(N_E=664\) e riduce \(F\) del \(98.6\%\).

Questi test misurano quanto \(W\) è soddisfatto, non l’errore di una PDE. Sono il modo più pulito per verificare che un campo di size prodotto da GCEE entra in \(W\) senza trascinarsi dietro idrodinamica o stime a posteriori.

### Gas impact e triple point 3D dentro hr (§4.5)

Snapshot, non la simulazione intera. In 2D si prende la densità discreta del gas impact e se ne usa il gradiente per size e aspect ratio; metriche \(\mu^r_9\) e \(\mu^h_1\). In 3D, generalizzazione del triple point nel miniapp aperto [Laghos](https://github.com/ceed/Laghos): \(\Omega=[0,7]\times[0,3]\times[0,3]\), cinque regioni, tre materiali, snapshot a \(t=5\). L’indicatore è trasferito su una mesh esaedrica uniforme con \(N_E=896\). Target di forma ideale, size dal gradiente dell’indicatore, \(\mu^r_{321}\) e \(\mu^h_{315}\). hr porta \(N_E\) a 3724. Errore volumetrico dell’indicatore rispetto alla mesh Laghos: 1.71 (uniforme), 0.62 (solo \(r\)), 0.078 (hr). Il codice ALE multi-materiale che produce il 2D non è aperto; Laghos sì.

### Poisson, shock allineato, trave e shear wall (2026, §8)

Kolev, Lazarov, Mittal, Schmidt, Tomov, Eng. Comput. 2026, [arXiv:2507.01917](https://arxiv.org/abs/2507.01917). L’obiettivo è una combinazione convessa di una misura di errore della PDE (\(F_p\)) e di una metrica TMOP (\(F_\mu\)). Gradienti aggiunti, regolarizzazione per convoluzione, solutore MMA. Implementazione dichiarata in MFEM; non ho trovato un miniapp con questi casi nel header di `mesh-optimizer`.

Tre misure, §4, nessuna delle quali è uno stimatore a priori classico:

- variazione locale della soluzione (la didascalia della Fig. 5 la chiama “element-local solution variation”; il §8.2 la chiama anche “solution-average”), §4.1;
- funzionale di carico, minimizzazione indiretta della norma dell’energia, §4.2;
- continuità del gradiente, esplicitamente ispirata allo stimatore di Zienkiewicz e Zhu (1987), §4.3.

Il paper dice che \(F_p\) può incorporare stimatori a posteriori esistenti, e rimanda a lavori futuri l’integrazione di stimatori più avanzati. Il confronto pubblicato è fra queste tre misure, più la mesh uniforme.

**Poisson 2D/3D, fronte circolare.** \(\Omega\) unitario, Dirichlet non omogeneo,

\[
u=\arctan\bigl(\beta(r(x)-r_c)\bigr),\quad
x_c=y_c=z_c=-0.05,\; r_c=0.7,\; \beta=20.
\]

È lo stesso fronte del test hr, con pendenza \(\beta=20\) invece di \(\alpha=200\). Mesh quad da \(4\times 4\) raffinata fino a tre volte, ordini 1, 2 e 3. Qualità: shape \(\mu_2\), \(W=I\). La Fig. 6 mostra il tasso \(O(h^{p+1})\) e una riduzione dell’errore \(L^2\) di circa un fattore 10 rispetto alla mesh uniforme. Pesi, raggio del filtro e iterazioni sono tarati a mano. Poi triangoli (da \(4\times 4\) quad, ogni quad in quattro triangoli e ogni triangolo in quattro), esaedri \(8\times 8\times 8\), tetraedri (hex \(4\times 4\times 4\) spezzati in 24 tet). Infine una mesh non strutturata con ordine della geometria diverso dall’ordine della soluzione (\(p_u=2\)). Sezione 8.1.

**Shock obliquo con allineamento.** \(\Omega=[0,1]^2\),

\[
u=\arctan\bigl(\alpha(x-0.5-0.2(y-0.5))\bigr).
\]

Target di orientazione \(\theta=\pi y(1-y)\cos(2\pi x)\), aspect ratio 1, skew \(\pi/2\), metrica shape+orientation \(\nu_{107}\). Mesh quad \(16\times 16\). Il peso \(\alpha\) della misura di errore vale \(0\), \(4\cdot 10^3\), \(4\cdot 10^4\), \(4\cdot 10^6\). Errore \(L^2\): 0.0034 (iniziale), 0.0015, 0.0007, 0.00025, 0.00032. Il peso più alto peggiora la qualità e rialza l’errore. Sezione 8.2, Tabella 1. È il test in cui qualità geometrica ed errore della PDE sono in trade-off esplicito.

**Trave elastica.** Dominio \(L=1\), \(H=0.1\), mesh \(10\times 4\), \(E=1\), \(\nu=0.3\), pressione \(t=-1.0\cdot 10^{2}\) sulla faccia superiore, incastro a sinistra. Si minimizza il funzionale di carico. Riferimento: spazio arricchito, ordine 4, mesh raffinata due volte. Spostamento in punta, ordine 1: riferimento \((-0.0907,-1.3766)\), iniziale \((-0.0662,-0.9999)\), ottimizzato \((-0.0736,-1.1776)\). Ordine 2 già vicino al riferimento. Sezione 8.3, Tabella 2.

**Shear wall con foro.** \(\Omega=[0,1]^2\setminus[1/3,2/3]^2\), 512 elementi, soluzione di ordine 2. Riferimento: ordine 4, 131072 elementi. Il modulo dello spostamento massimo passa da 11.6611 a 11.7329, contro 11.7585 del riferimento. Sezione 8.3, Tabella 3.

## Miniapp MFEM

TMOP non è un example numerato (`ex1`, `ex2`, …). Sta in

- `miniapps/meshing/mesh-optimizer.cpp` (seriale),
- `miniapps/meshing/pmesh-optimizer.cpp` (parallelo, stessi casi con `mpirun`).

Repository [github.com/mfem/mfem](https://github.com/mfem/mfem), licenza BSD-3. Compilazione: `make mesh-optimizer` dentro `miniapps/meshing`. La mesh di default, se si omette `-m`, è `icf.mesh`. I comandi qui sotto sono quelli nell’header di `master`; possono spostarsi fra una release e l’altra. Miniapp vicini, non catalogati voce per voce: `pmesh-fitting.cpp`, `fit-node-position.cpp` (fitting di superfici con \(\mu_2\) / \(\mu_{302}\)).

`-tid` è il target, `-mid` la metrica. L’help a riga di comando documenta i target 1–5; i casi 6–11 ci sono nello `switch` ma non in quella stringa di help.

| tid | Costruzione | Cosa rappresenta |
| --- | --- | --- |
| 1 | `IDEAL_SHAPE_UNIT_SIZE` | Forma ideale, size unitaria |
| 2 | `IDEAL_SHAPE_EQUAL_SIZE` | Forma ideale, size uguale |
| 3 | `IDEAL_SHAPE_GIVEN_SIZE` | Forma ideale, size iniziale |
| 4 | `GIVEN_FULL` analitico (`HessianCoefficient`) | Jacobiano analitico in spazio fisico; la formula dipende da `-mid` |
| 5 | size discreta (`DiscreteAdaptTC`) | Forma ideale, size da `size_indicator` |
| 6 | size + aspect ratio 2D | Da `material_indicator_2d` (interfaccia sinusoidale ruotata di \(22.5^\circ\)) |
| 7 | aspect ratio discreto 3D | `discrete_aspr_3d`: lati \(1\), \(1+5y\), \(1+10z\) |
| 8 | size + orientazione 2D | Size costante \(0.01\), \(\theta=\pi y(1-y)\cos(2\pi x)\) |
| 9 | hr, anello, sola size | Size \(0.001\) / \(0.01\) fra i raggi \(0.15\) e \(0.35\) |
| 10 | hr, anello, size + aspect ratio | Come §4.2 del paper hr, raggi \(0.2\) e \(0.3\) |
| 11 | hr, onda sinusoidale ruotata | Size + aspect ratio, rotazione \(15.5^\circ\) |

`size_indicator` (tid 5) è una corona semicircolare: centro \((0,0.5)\) in 2D, raggi \(0.45\) e \(0.55\), `tanh` con sharpness 30. Il rapporto small/big è 9 in 2D e 27 in 3D. Non è la banda sinusoidale dell’IMR.

`HessianCoefficient` (tid 4): se la metrica non è 14, 36 o 85, corona circolare fra i raggi \(0.15\) e \(0.35\) centrata in \((0.5,0.5)\), con \(K_{00}\) variabile e \(K_{11}=1\) (shape anisotropa). Per mid 14 e 36, allineamento con lo stesso \(\theta\) del tid 8, scalato di \(0.1\). Per mid 85, shape+allineamento su un’interfaccia pesata.

Comandi dell’header, raggruppati. `square01.mesh` e `cube.mesh` sono il quadrato e il cubo unitari distribuiti col miniapp; `blade.mesh`, `icf.mesh`, `jagged.mesh`, `stretched2D.mesh`, `amr-quad-q2.mesh` stanno nella stessa directory. Le mesh marcate `mfem/data` stanno nel repository dati separato.

**r-adattività / size (i candidati per GCEE)**

- Forma analitica adattata: `mesh-optimizer -m square01.mesh -o 2 -rs 2 -mid 2 -tid 4 -ni 200 -bnd -qt 1 -qo 8`
- Size+orientazione analitica: `... -mid 14 -tid 4 -ni 100 -bnd -qt 1 -qo 8`
- Shape+orientazione (automatic differentiation): `... -o 3 -mid 85 -tid 4 -ni 100 -bnd -qt 1 -qo 8 -rtol 1e-6`
- hr: `-m square01.mesh -o 2 -tid 9 -ni 50 -li 20 -hmid 55 -mid 7 -hr` (size); `-tid 10` (size+aspect); `-tid 11 -hmid 58` (onda)
- Size discreta: `-m square01.mesh -o 2 -rs 2 -mid 94 -tid 5 -ni 50 -qo 4 -nor`
- Size discreta su mesh non conforme: `-m amr-quad-q2.mesh -o 2 -rs 2 -mid 94 -tid 5 -ni 50 -qo 4 -nor`
- Size discreta 3D: `-m cube.mesh -o 2 -rs 2 -mid 321 -tid 5 -ls 3 -nor -pa -rtol 1e-8`
- Size+aspect ratio da interfaccia: `-m square01.mesh -o 2 -rs 2 -mid 7 -tid 6 -ni 100`
- Size+orientazione discreta: `-m square01.mesh -o 2 -rs 2 -mid 36 -tid 8 -qo 4 -nor -rtol 1e-6`
- Aspect ratio 3D: `-m cube.mesh -o 2 -rs 2 -mid 302 -tid 7 -ni 20 -bnd -qt 1 -qo 8`
- Periodico + size: `periodic-square.mesh`, `periodic-tri.mesh`, `periodic-cube.mesh` con `-tid 5`

**Geometrici, untangling, limiting**

- Pala, forma: `-m blade.mesh -o 4 -mid 2 -tid 1 -ni 30 -ls 3 -art 1 -bnd -qt 1 -qo 8`. Limitata: stesso comando con `-lc 5000` (il \(\delta=5000\) del §5.1).
- ICF (mesh di default): equal size `-o 3 -mid 80 -bec -tid 2`; initial size `-o 3 -mid 9 -tid 3`; shape `-o 3 -mid 1 -tid 1`; shape limitata `-lc 10`; combo anelli `-cmb 1`. Attenzione: nel paper il caso equal-size e initial-size usa \(\mu_9\), mentre il miniapp usa `mid 80` e `mid 9`. La quadratura del paper è 6 punti per direzione; diversi comandi del miniapp usano `-qo 8`. I valori pubblicati dell’obiettivo non vanno attesi bit a bit.
- Sfera pizzicata: mesh `ball-pert.mesh` nel repo dati, `-o 4 -mid 303 -tid 1`.
- Untangling 2D: `-m jagged.mesh -o 2 -mid 22 -tid 1` (metrica di untangling) oppure `-mid 4 -btype 1` (shifted barrier). Untangling 3D: `cube-holes-inv.mesh`, `-mid 313`.
- Limiting adattivo: `-m stretched2D.mesh -o 2 -mid 2 -tid 1 -alc 0.5`.
- Misto tet/hex con limiting: `fichera-mixed-p2.mesh`, `-mid 301 -tid 1 -lc 0.25`.
- Non conforme, shape + equal size: `amr-quad-q2.mesh`, `-mid 9 -tid 2`.

Il peso degli anelli ICF del §5.2.4 è reimplementato come `weight_fun` in `mesh-optimizer.hpp`.

## Riproducibilità

| Banco | Codice | Mesh e dati |
| --- | --- | --- |
| Blade, ICF, sfera, jagged, square, cube, tid 4–11 | MFEM, aperto | Mesh nel miniapp o in `mfem/data`. I numeri del SISC 2019 non sono un gold file: metrica e quadratura del miniapp non coincidono sempre col paper |
| Banda sinusoidale IMR | Algoritmo descritto; la \(\eta\) del §6.1 non è il default di tid 5 | Va reimplementata. Il trasporto di \(\eta\) (advection vs interpolazione) cambia il transitorio |
| Poisson hr e anello analitico | Algoritmo in MFEM (`TMOPHRSolver`, tid 9–11 per l’anello) | Il Poisson del §4.1 non è un comando dell’header. Formula da ricontrollare sul PDF |
| Poisson / elastici 2026 | “implemented using MFEM”; non è un sample run del miniapp | Formule e tabelle bastano per reimplementare. Pesi e filtro sono tarati, non fissati da una ricetta unica |
| Triple Point, gas impact, ball impact in BLAST | BLAST chiuso | Il problema idrodinamico è citato; la catena ALE+remap no. Laghos è aperto e copre solo lo snapshot lagrangiano 3D del §4.5 hr |
| Knupp 2012 / Mesquite | Mesquite è il codice citato | I casi MATLAB 2007 e le mesh del talk 2010 non risultano pacchettizzati col paper |

## Cosa è già stato usato per confrontare stimatori

- hr 2022, §4.1: confronta \(h\), \(r\) e \(hr\) sulla stessa soluzione esatta, in norma dell’energia. Il target di size e aspect ratio viene dal gradiente esatto (oracolo), non da uno stimatore. Lo stimatore entra solo nella decisione di raffinare.
- hr 2022, §3.6: dice che quello stimatore TMOP si può sostituire. Non lo fa.
- Computers & Fluids 2020: nomina “known error estimates” come ingresso ammesso. I numeri pubblicati usano volume fraction e interfacce. L’errore riportato è la diffusione dell’interfaccia rispetto al lagrangiano, non l’errore di uno stimatore.
- Eng. Comput. 2026, §8: unico confronto pubblicato di più misure di errore dentro un’ottimizzazione TMOP. Le tre misure sono a posteriori o energetiche (variazione locale, funzionale di carico, continuità del gradiente alla Zienkiewicz–Zhu). Non ci sono stimatori a priori, né un confronto GCEE-style “stesso \(W\), stimatori diversi”.

Uno stimatore a priori, in questa letteratura, sostituisce la costruzione di \(\zeta\) e \(\rho\) (e, se serve, di \(\theta\)). Non sostituisce \(\mu\).

## Shortlist per confrontare uno stimatore a priori

Sei punti, in ordine di implementazione. I primi quattro bastano per un confronto pulito; gli ultimi due coprono anisotropia e una feature non liscia senza aprire BLAST.

1. **Poisson a fronte circolare, ricetta 2026, §8.1.** \(\Omega=[0,1]^d\), \(u=\arctan(\beta(r-0.7))\), centro \((-0.05,-0.05)\), \(\beta=20\), mesh quad da \(4\times 4\), ordini 1–3, errore \(L^2\) contro la mesh uniforme. Formula non ambigua, soluzione esatta, stesso fronte già usato dal gruppo TMOP. GCEE costruisce la size (e l’aspect ratio) al posto della misura di variazione locale o del funzionale di carico. Il fattore “circa 10” sull’errore è il riferimento pubblicato, con la riserva che i pesi sono tarati.

2. **Stesso Poisson, protocollo hr 2022, §4.1, solo ramo \(r\).** Norma dell’energia, target (15) shape+size, metrica \(\mu_9\), partenza \(4\times 4\). Qui il confronto naturale è: gradiente esatto contro GCEE contro, se si vuole un a posteriori già nella famiglia, la continuità di gradiente del 2026. Tenere \(\alpha=200\) e \(\beta=20\) come due sharpness, non come lo stesso test. Ricontrollare l’equazione (13) sul PDF prima di fissare \(\alpha\).

3. **Anello analitico, hr §4.2 e `mesh-optimizer -tid 9` / `-tid 10`.** Nessuna PDE. Serve a verificare che il campo di size di GCEE entra in \(W\) e che il funzionale scende, prima di misurare un errore di soluzione. Tid 9 è sola size; tid 10 aggiunge l’aspect ratio angolare. Aperto, un comando.

4. **Banda sinusoidale, IMR §6.1.** \(16\times 16\), ordine 3, tre target (size \(\mu_7\), shape \(\mu_{58}\), shape+size \(\mu_9\)) con riduzioni pubblicate \(47\%\), \(22\%\), \(65\%\). È il test storico di r-adattività su campo discreto. GCEE può sostituire \(\eta\), a patto di dichiarare che i numeri IMR restano il riferimento del campo analitico, non dello stimatore. La \(\eta\) va codificata: non è `size_indicator`.

5. **Shock obliquo allineato, 2026 §8.2.** Stesso \(\theta=\pi y(1-y)\cos(2\pi x)\) di `HessianCoefficient` (mid 14 e 36) e del tid 8. Tabella 1: errore \(L^2\) da 0.0034 a 0.00025, poi risalita se il peso dell’errore schiaccia la qualità. È il banco giusto se GCEE produce un target anisotropo e si vuole vedere il trade-off con \(\mu\), non solo la size.

6. **Interfaccia del miniapp, `-tid 6`, oppure size discreta `-tid 5`.** Tid 6 (`material_indicator_2d`, metrica \(\mu_7\)) è una feature non liscia, ruotata, con size e aspect ratio dal gradiente dell’indicatore: il surrogato aperto del gas impact, senza BLAST. Tid 5 (`size_indicator` su `square01.mesh`, \(\mu_{94}\)) è il caso minimo di forma ideale più size assegnata, utile come harness: si sostituisce solo la funzione indicatrice con l’output di GCEE.

Da lasciare fuori dal primo confronto: pala, ICF, sfera pizzicata, untangling (`jagged`, `cube-holes-inv`). Misurano forma e validità, non uno stimatore. Triple Point, gas impact completo e ball impact richiedono BLAST. Lo snapshot Laghos del §4.5 è riproducibile a metà (Laghos è aperto, il trasferimento e il codice ALE no) e va dopo i sei punti sopra.
