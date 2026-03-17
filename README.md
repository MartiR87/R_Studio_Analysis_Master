# Analisi dati clinici e molecolari in R

Questa repository contiene lo script R utilizzato per la tesina Master di II livello in Bioinformatica e Data Science, dal titolo:

“MSI-H e varianti passenger nei tumori del colon-retto, dell’endometrio e dell’ovaio: analisi di una coorte di pazienti oncologici per una possibile stratificazione terapeutica”.

Lo script riguarda l’analisi di dati clinici e molecolari di pazienti con tumore del colon-retto, ovaio ed endometrio, e include tutti i passaggi necessari per:

- la pulizia e normalizzazione dei dati clinici;
- selezione dei pazienti con tumori di interesse (ovaio, colon-retto ed endometrio)
- l’analisi della distribuzione dei pazienti in base al valore di MSI (Microsatellite Instability);
- il filtraggio dei pazienti con varianti hotspot BRCA1/2;
- l’analisi delle varianti nei geni MMR;
- la generazione di grafici e tabelle riassuntive per pazienti MSI-H (MSI>20) e MSI-L (MSI<20);
- l’identificazione di varianti novel in regioni omopolimeriche nei geni BRCA1/2
- Approfondimento per Endometrio e confronto tra MMR-mutati Vs WT


## Contenuto

- **Script_R_Comandi.R**: script completo in R che esegue:
  
  - Pulizia e normalizzazione dei dati clinici
  - Selezione dei pazienti con i tre tumori studiati (ovaio, colon-retto, endometrio)
  - Analisi della distribuzione dei pazienti in base all'MSI (Microsatellite Instability)
  - Filtraggio pazienti con varianti hotspot in BRCA1/2 già note in letteratura scientifica
  - Analisi delle varianti nei geni MMR nei pazienti con MSI-H (MSI>20)
  - Grafici e tabelle riassuntive per pazienti MSI-H e MSI-L
  - Analisi delle varianti novel in regioni omopolimeriche nei geni BRCA1/2


## Dati

I dati utilizzati nello script includono:

- **data_clinical_sample**: informazioni cliniche dei pazienti
- **data_mutations_extended**: dati mutazionali NGS filtrati
- **data_clinical_patient**: informazioni generali dei pazienti

> Nota: i dataset non sono inclusi in questa repository per motivi di privacy e non sono pubblicamente distribuibili. L’accesso è stato fornito alla commissione tramite condivisione privata.

## Istruzioni per l’uso

1. Clonare la repository sul proprio computer:
   ```bash
   git clone https://github.com/MartiR87/R_Studio_Analysis_Master.git
   
2. Aprire lo script in RStudio:
Assicurarsi di avere installato i pacchetti R necessari: install.packages(c("dplyr", "ggplot2", "scales", "ggsignif", "dunn.test", "tidytext", "readr"))

4. Eseguire i comandi passo passo per generare grafici, tabelle e analisi.
   
6. Link
Repository GitHub: https://github.com/MartiR87/R_Studio_Analysis_Master

Licenza

Questa repository è destinata a uso accademico e non contiene dati sensibili dei pazienti.
