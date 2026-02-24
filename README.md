# Analisi dati clinici e molecolari in R

Questa repository contiene lo script R utilizzato per la tesina di Master di II livello in Bioinformatica e Data Science, dal titolo:

“MSI-H e varianti passenger nei tumori del colon-retto, dell’endometrio e dell’ovaio: analisi di una coorte di pazienti oncologici per una possibile stratificazione terapeutica”.

Lo script riguarda l’analisi di dati clinici e molecolari di pazienti con tumore del colon-retto, ovaio ed endometrio, e include tutti i passaggi necessari per:

- la pulizia e normalizzazione dei dati clinici;
- l’analisi della distribuzione di MSI (Microsatellite Instability);
- il filtraggio dei pazienti con varianti BRCA1/2;
- l’analisi delle mutazioni nei geni MMR;
- la generazione di grafici e tabelle riassuntive per pazienti MSI-H e MSI-L;
- l’identificazione di varianti novel e note nei geni BRCA1/2.


## Contenuto

- **Script_R_Comandi.R**: script completo in R che esegue:
  - Pulizia e normalizzazione dei dati clinici
  - Analisi della distribuzione MSI (Microsatellite Instability)
  - Filtraggio pazienti con varianti BRCA1/2
  - Analisi delle mutazioni nei geni MMR
  - Grafici e tabelle riassuntive per pazienti MSI-H e MSI-L
  - Analisi delle varianti novel e note nei geni BRCA1/2

## Dati

I dati utilizzati nello script includono:

- **data_clinical_sample**: informazioni cliniche dei pazienti
- **data_mutations_extended**: dati mutazionali NGS filtrati

> Nota: i dati reali non sono inclusi in questa repository per motivi di privacy.

## Istruzioni per l’uso

1. Clonare la repository sul proprio computer:
   ```bash
   git clone https://github.com/MartiR87/R_Studio_Analysis_Master.git
2. Aprire lo script in RStudio:
Assicurarsi di avere installato i pacchetti R necessari: install.packages(c("dplyr", "ggplot2", "scales", "ggsignif", "dunn.test", "tidytext", "readr"))
3. Eseguire i comandi passo passo per generare grafici, tabelle e analisi.
4. Link
Repository GitHub: https://github.com/MartiR87/R_Studio_Analysis_Master

Licenza

Questa repository è destinata a uso accademico e non contiene dati sensibili dei pazienti.
