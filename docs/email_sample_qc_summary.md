Subject: Cinque Terre WGS – sample-level QC summary

Ciao,

abbiamo completato il QC preliminare a livello di individuo sui 50 campioni WGS delle Cinque Terre, lavorando sul VCF finale post-VQSR e limitando le valutazioni agli autosomi 1–22.

Abbiamo usato come indicatori principali profondità media, missingness e tasso di eterozigosità, con numero di varianti, singleton burden e Ts/Tv come diagnostica secondaria. I campioni flaggati sono stati poi controllati cromosoma per cromosoma e con i campi di genotipo GT/DP/GQ/AD.

Quattro campioni mostrano anomalie riproducibili e genome-wide:

- TSBC6060: eccesso di eterozigosità su 22/22 autosomi, associato a forte squilibrio dell'allele balance nei genotipi eterozigoti.
- TSBC6389: stesso pattern di eccesso di eterozigosità su 22/22 autosomi e allele balance anomalo, con anche depth/GQ più bassi rispetto al resto della coorte.
- TSBD8047: eccesso di eterozigosità su 22/22 autosomi con allele balance anomalo.
- TSBD8199: quadro diverso e chiaramente tecnico, con profondità molto bassa, missingness elevata, GQ ridotto e profilo di eterozigosità anomalo su quasi tutti gli autosomi.

Per i primi tre il pattern è compatibile con un artefatto tecnico a livello di campione; non lo definirei però come contaminazione certa senza un test dedicato sui BAM/CRAM. Per TSBD8199 la classificazione più appropriata è low-coverage/low-quality sample.

Sulla base di questi risultati abbiamo deciso di escludere i quattro campioni prima del QC a livello di variante. La decisione non deriva da una soglia automatica né dal voler raggiungere un numero prefissato di campioni: è il risultato della convergenza di più segnali QC indipendenti e genome-wide.

Sto preparando anche due figure sintetiche per mostrare visivamente gli outlier: depth vs missingness e heterozygosity vs allelic imbalance.

A presto
