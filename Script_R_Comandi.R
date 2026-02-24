# ============================================================
# ANALISI DESCRITTIVA DEI DATI
# ============================================================
# Obiettivo: filtrare i pazienti con tumore del colon-retto,
# endometrio e ovaio e con MSI elevato (MSI > 20), 
# visualizzare distribuzione dei tumori e boxplot MSI-H
# ============================================================

# ---------------------------- LIBRERIE ------------------------ 
library(dplyr)       # manipolazione dati
library(ggplot2)     # grafici
library(scales)      # formattazione percentuali
library(ggsignif)    # barre di significatività

# ---------------------------- PARAMETRI ------------------------ 
tumori_interesse <- c("Colorectum", "Ovary", "Endometrium")

# ---------------------------- 1. PULIZIA DATI ------------------------ 
# Mantieni solo valori numerici nella colonna MSI (V3)
data_clean <- data_clinical_sample[grepl("^[0-9.]+$", data_clinical_sample$V3), ]
data_clean$V3 <- as.numeric(data_clean$V3)

# ---------------------------- 2. GRAFICO DISTRIBUZIONE TUMORI ------------------------ 
# Calcola conteggi e percentuali dei tumori selezionati
totali_per_tumore <- data_clean %>%
  filter(V14 %in% tumori_interesse) %>%
  count(Tumore = V14) %>%
  mutate(
    perc = n / sum(n),
    label = paste0(round(perc * 100, 1), "%")
  )

# Grafico a barre della distribuzione dei tumori
ggplot(totali_per_tumore, aes(x = Tumore, y = perc, fill = Tumore)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  labs(
    title = "Distribuzione dei tumori nel dataset completo",
    x = "Tipo di tumore",
    y = "Percentuale sul totale"
  ) +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# ---------------------------- 3. FILTRO MSI-H ------------------------ 
# Seleziona pazienti con MSI > 20 e tumori di interesse
msi_alti <- data_clean %>%
  filter(
    V14 %in% tumori_interesse,
    V3 > 20,
    !is.na(V1)
  ) %>%
  select(Paziente_ID = V1, Tumore = V14, MSI = V3)

# ---------------------------- 4. NORMALIZZAZIONE MSI-H PER TUMORE ------------------------ 
# Conta pazienti MSI-H per tumore
msi_per_tumore <- msi_alti %>% count(Tumore)

# Calcola percentuale normalizzata rispetto al totale dei pazienti
msi_normalizzato <- msi_per_tumore %>%
  left_join(totali_per_tumore %>% select(Tumore, n_Totale = n), by = "Tumore") %>%
  mutate(
    perc_MSI_H = n / n_Totale,
    label = paste0(round(perc_MSI_H * 100, 1), "%")
  )

# Grafico percentuale MSI-H normalizzata
ggplot(msi_normalizzato, aes(x = Tumore, y = perc_MSI_H, fill = Tumore)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), vjust = -0.7, size = 3) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, max(msi_normalizzato$perc_MSI_H)*1.2)) +
  theme_minimal() +
  labs(
    title = "Percentuale di pazienti MSI-H per tumore",
    x = "Tipo di tumore",
    y = "Percentuale MSI-H"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

# ---------------------------- 5. BOXPLOT MSI-H CON SIGNIFICATIVITA' ------------------------ 
# Definizione delle coppie da confrontare
comparisons <- list(
  c("Colorectum", "Endometrium"),
  c("Colorectum", "Ovary"),
  c("Endometrium", "Ovary")
)

# Boxplot pulito con jitter e barre di significatività Wilcoxon
ggplot(msi_alti, aes(x = Tumore, y = MSI, color = Tumore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("skyblue", "lightgreen", "salmon")) +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.1,
    textsize = 4
  ) +
  labs(
    title = "Distribuzione MSI per tipo di tumore (MSI > 20) con significativitÃ ",
    x = "Tipo di tumore",
    y = "MSI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )



# ============================================================
# SECONDA PARTE: ANALISI DELLA NORMALITA'
# ============================================================
# Obiettivo: verificare se la distribuzione dei valori di MSI
# nei diversi tumori segue una distribuzione normale,
# e confrontare le medie tra gruppi tramite test statistici.
# ============================================================

# ---------------------------- LIBRERIE AGGIUNTIVE ------------------------ 
# Assicurati di aver già caricato:
library(dplyr)
library(ggplot2)
library(ggsignif)
library(dunn.test)   # per test post-hoc di Dunn

# ---------------------------- 1. ISTOGRAMMI E DENSITA' ------------------------ 
# Visualizzazione della distribuzione MSI-H per tumore
ggplot(msi_alti, aes(x = MSI, fill = Tumore)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", alpha = 0.5) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Tumore, scales = "free") +
  labs(
    title = "Distribuzione MSI-H per tipo di tumore",
    x = "MSI",
    y = "DensitÃ "
  ) +
  theme_minimal()

# ---------------------------- 2. TEST SHAPIRO-WILK ------------------------ 
# Verifica formalmente la normalità per ciascun tumore
shapiro_results <- by(msi_alti$MSI, msi_alti$Tumore, shapiro.test)
shapiro_results
# Nota: p-value < 0.05 indica deviazione dalla normalità

# ---------------------------- 3. Q-Q PLOT ------------------------ 
# Q-Q plot dei valori MSI-H per ciascun tumore
ggplot(msi_alti, aes(sample = MSI)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~ Tumore, scales = "free") +
  labs(
    title = "Q-Q plot dei valori MSI-H per tipo di tumore",
    x = "Quantili teorici",
    y = "Quantili osservati"
  ) +
  theme_minimal()

# ---------------------------- 4. TEST KRUSKAL-WALLIS ------------------------ 
# Test non parametrico tra gruppi tumorali (MSI-H)
kruskal_results <- kruskal.test(MSI ~ Tumore, data = msi_alti)
kruskal_results
# Nota: p-value < 0.05 indica differenze significative tra almeno due gruppi

# ---------------------------- 5. TEST POST-HOC DI DUNN ------------------------ 
# Con correzione Bonferroni per confronti multipli
dunn_results <- dunn.test(msi_alti$MSI, msi_alti$Tumore, method = "bonferroni")
dunn_results
# I risultati mostrano quali coppie di tumori differiscono significativamente



# ============================================================
# TERZA PARTE: PAZIENTI CON VARIANTI BRCA1/2 E ANALISI MSI
# ============================================================
# Obiettivo: filtrare pazienti con mutazioni BRCA1/2 nei tumori
# di interesse e analizzare la distribuzione MSI (alto/basso)
# ============================================================

# ---------------------------- 1. PAZIENTI CON VARIANTI BRCA1/2 ------------------------ 
# Seleziona solo i pazienti con mutazioni BRCA1 o BRCA2
brca_patients <- data_mutations_extended %>%
  filter(V2 %in% c("BRCA1","BRCA2")) %>%
  select(Paziente_ID = V1)  # salva solo l'ID del paziente

# ---------------------------- 2. UNIONE CON DATI CLINICI ------------------------ 
# Recupera tumore e MSI dei pazienti BRCA1/2
brca_clinical <- data_clinical_sample %>%
  filter(V14 %in% c("Colorectum", "Ovary", "Endometrium")) %>%  # solo i 3 tumori
  filter(V1 %in% brca_patients$Paziente_ID) %>%                  # solo pazienti BRCA1/2
  filter(!is.na(V3)) %>%                                         # rimuove MSI NA
  select(Paziente_ID = V1, Tumore = V14, MSI = V3)

# ---------------------------- 3. CLASSIFICAZIONE MSI ------------------------ 
# Definisce MSI alto (>20) o basso (<20)
brca_clinical <- brca_clinical %>%
  mutate(MSI_status = ifelse(MSI > 20, "MSI > 20", "MSI < 20"))

# ---------------------------- 4. CONTEGGIO PAZIENTI ------------------------ 
# Numero di pazienti per tumore e stato MSI
brca_count <- brca_clinical %>%
  group_by(Tumore, MSI_status) %>%
  summarise(Pazienti = n(), .groups = "drop")

# ---------------------------- 5. NORMALIZZAZIONE ------------------------ 
# Calcola la percentuale di pazienti per tumore
brca_normalizzato <- brca_count %>%
  group_by(Tumore) %>%
  mutate(Percentuale = round(Pazienti / sum(Pazienti) * 100, 1)) %>%
  ungroup()

# ---------------------------- 6. GRAFICO A BARRE ------------------------ 
# Visualizza numero e percentuale di pazienti BRCA1/2 per tumore e MSI
ggplot(brca_normalizzato, aes(x = Tumore, y = Pazienti, fill = MSI_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0(Percentuale, "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  labs(
    title = "Distribuzione normalizzata di pazienti BRCA1/2 per tumore e MSI",
    x = "Tipo di tumore",
    y = "Numero di pazienti",
    fill = "MSI status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )



# ============================================================
# QUARTA PARTE: ANALISI VARIANTI BRCA E GENI MMR
# ============================================================
# Obiettivo: analizzare le varianti ricorrenti nei geni BRCA1/2
# e mutazioni nei geni MMR nei pazienti MSI-H e MSI-L
# ============================================================

# ----------------------------
# 1. PARAMETRI PRINCIPALI
# ----------------------------
tumori_interesse <- c("Colorectum", "Ovary", "Endometrium")
varianti_target <- c("rs80357569", "rs80357522", "rs80359306",
                     "rs80359479", "rs80359507", "rs397507419")
geni_MMR <- c("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM")

# ----------------------------
# 2. Dataset clinico pulito con MSI numerico
# ----------------------------
data_clean <- data_clinical_sample %>%
  filter(V14 %in% tumori_interesse,
         grepl("^[0-9.]+$", V3)) %>%
  mutate(
    Paziente_ID = V1,
    Tumore = V14,
    MSI = as.numeric(V3),
    Stato_MSI = ifelse(MSI > 20, "MSI-H", "MSS")
  )

# ----------------------------
# 3. Pazienti MSI-H con varianti BRCA1/2 target
# ----------------------------
msi_alti_varianti <- data_clean %>%
  filter(Stato_MSI == "MSI-H") %>%
  inner_join(
    data_mutations_extended %>%
      filter(V2 %in% c("BRCA1","BRCA2"),
             V15 %in% varianti_target) %>%
      select(Paziente_ID = V1, Gene = V2, Variante = V15),
    by = "Paziente_ID"
  ) %>%
  select(Paziente_ID, Tumore, MSI, Gene, Variante) %>%
  arrange(Tumore, MSI)

# ----------------------------
# 4. Pazienti MSI-L con varianti BRCA1/2 target
# ----------------------------
msi_bassi_varianti <- data_clean %>%
  filter(Stato_MSI == "MSS") %>%
  inner_join(
    data_mutations_extended %>%
      filter(V2 %in% c("BRCA1","BRCA2"),
             V15 %in% varianti_target) %>%
      select(Paziente_ID = V1, Gene = V2, Variante = V15),
    by = "Paziente_ID"
  ) %>%
  select(Paziente_ID, Tumore, MSI, Gene, Variante) %>%
  arrange(Tumore, MSI)

# ----------------------------
# 5. Conteggio e normalizzazione pazienti MSI-H per tumore, gene e variante
# ----------------------------
msi_brca_normalizzato <- msi_alti_varianti %>%
  group_by(Tumore, Gene, Variante) %>%
  summarise(N_pazienti = n_distinct(Paziente_ID), .groups = "drop") %>%
  group_by(Tumore) %>%
  mutate(Percentuale = round(N_pazienti / sum(N_pazienti) * 100, 1)) %>%
  ungroup()

# ----------------------------
# 6. Grafico MSI-H
# ----------------------------
ggplot(msi_brca_normalizzato, aes(x = Variante, y = N_pazienti, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = paste0(Percentuale, "%")),
            position = position_dodge(width = 0.8),
            vjust = 0.5, hjust = 0.5, angle = 90,
            color = "black", size = 3) +
  facet_wrap(~ Tumore) +
  labs(
    title = "Pazienti MSI > 20 con varianti BRCA1/2",
    x = "Variante",
    y = "Numero di pazienti",
    fill = "Gene"
  ) +
  scale_fill_manual(values = c("BRCA1" = "#E69F00", "BRCA2" = "#56B4E9")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------
# 7. Conteggio e normalizzazione pazienti MSI-L
# ----------------------------
msi_basso_normalizzato <- msi_bassi_varianti %>%
  group_by(Tumore, Gene, Variante) %>%
  summarise(N_pazienti = n_distinct(Paziente_ID), .groups = "drop") %>%
  group_by(Tumore) %>%
  mutate(Percentuale = round(N_pazienti / sum(N_pazienti) * 100, 1)) %>%
  ungroup()

# ----------------------------
# 8. Grafico MSI-L
# ----------------------------
ggplot(msi_basso_normalizzato, aes(x = Variante, y = N_pazienti, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = paste0(Percentuale, "%")),
            position = position_dodge(width = 0.8),
            vjust = 0.5, hjust = 0.5, angle = 90,
            color = "black", size = 3) +
  facet_wrap(~ Tumore) +
  labs(
    title = "Pazienti MSI < 20 con varianti BRCA1/2",
    x = "Variante",
    y = "Numero di pazienti",
    fill = "Gene"
  ) +
  scale_fill_manual(values = c("BRCA1" = "#E69F00", "BRCA2" = "#56B4E9")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------
# 9. Analisi mutazioni geni MMR (MLH1, MSH2, MSH6, PMS2, EPCAM)
# ----------------------------
mutazioni_MMR <- data_mutations_extended %>%
  filter(V2 %in% geni_MMR) %>%
  rename(Paziente_ID = V1, Gene = V2)

mutazioni_MMR_con_tumore <- mutazioni_MMR %>%
  inner_join(data_clean %>% select(Paziente_ID, Tumore, Stato_MSI),
             by = "Paziente_ID")

conteggio_MMR <- mutazioni_MMR_con_tumore %>%
  group_by(Tumore, Gene, Stato_MSI) %>%
  summarise(Pazienti_con_mutazione = n_distinct(Paziente_ID),
            .groups = "drop") %>%
  left_join(
    data_clean %>%
      group_by(Tumore) %>%
      summarise(Tot_Pazienti = n_distinct(Paziente_ID), .groups = "drop"),
    by = "Tumore"
  ) %>%
  mutate(Percentuale = round(100 * Pazienti_con_mutazione / Tot_Pazienti, 1)) %>%
  arrange(Tumore, Stato_MSI, desc(Pazienti_con_mutazione))

# ----------------------------
# 10. Visualizzazione grafica mutazioni MMR
# ----------------------------
library(tidytext)  # per reorder_within

ggplot(conteggio_MMR,
       aes(x = reorder_within(Gene, Percentuale, Tumore),
           y = Percentuale,
           fill = Stato_MSI)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = paste0(Percentuale, "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 4, angle = 90) +
  scale_x_reordered() +
  facet_wrap(~ Tumore, scales = "free_x") +
  labs(title = "Percentuale di pazienti con mutazioni nei geni MMR",
       subtitle = "Normalizzate sul totale dei pazienti per tumore",
       x = "Gene MMR",
       y = "Percentuale di pazienti (%)",
       fill = "Stato MSI") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# ============================================================
# QUINTA PARTE: Analisi omopolimeri in BRCA1/2 (MSI > 20)
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)  # per eventuale esportazione CSV

# ---------------------------
# 1. Pulizia dati clinici
# ---------------------------
clinical_selected <- data_clinical_sample %>%
  select(V1, V3, V14) %>%
  rename(
    SAMPLE_ID = V1,
    MSI = V3,
    NEOPLASM = V14
  ) %>%
  filter(NEOPLASM %in% c("Colorectum", "Ovary", "Endometrium")) %>%
  mutate(MSI = as.numeric(MSI))

# Calcolo totale pazienti per tumore (per normalizzazione)
totali_per_tumore <- clinical_selected %>%
  group_by(NEOPLASM) %>%
  summarise(Totale_pazienti = n_distinct(SAMPLE_ID), .groups = "drop")

# ---------------------------
# 2. Pulizia dati mutazioni e unione con clinici
# ---------------------------
mutations_selected <- data_mutations_extended %>%
  select(V1, V2, V5, V6, V7, V8, V11, V15, V35, V36) %>%
  rename(
    SAMPLE_ID = V1,
    Hugo_Symbol = V2,
    NCBI_Build = V5,
    Chromosome = V6,
    Start_Position = V7,
    End_Position = V8,
    Variant_Type = V11,
    dbSNP_RS = V15,
    HGVSc = V35,
    HGVSp = V36
  )

merged_data <- left_join(clinical_selected, mutations_selected, by = "SAMPLE_ID")

# ---------------------------
# 3. Filtri principali: MSI > 20, BRCA1/2, tipo variante INS/DEL
# ---------------------------
filtered_final <- merged_data %>%
  filter(
    MSI > 20,
    Hugo_Symbol %in% c("BRCA1", "BRCA2"),
    Variant_Type %in% c("INS", "DEL")
  )
cat("âœ… Varianti filtrate:", nrow(filtered_final), "\n")

# ---------------------------
# 4. Rimozione varianti note già analizzate
# ---------------------------
varianti_da_escludere <- c("rs80357569", "rs80357522", "rs80359306",
                           "rs80359479", "rs80359507", "rs397507419")

dataset_filtrato <- filtered_final %>%
  filter(!(trimws(dbSNP_RS) %in% varianti_da_escludere))
cat("âœ… Righe totali nel dataset filtrato (senza varianti note):", nrow(dataset_filtrato), "\n")

# ---------------------------
# 5. Separazione varianti novel e con rsID valido
# ---------------------------
varianti_novel <- dataset_filtrato %>% filter(!grepl("^rs", dbSNP_RS))
varianti_con_rs <- dataset_filtrato %>% filter(grepl("^rs", dbSNP_RS))
cat("Varianti con rsID:", nrow(varianti_con_rs), "\n")
cat("Varianti novel:", nrow(varianti_novel), "\n")

# ---------------------------
# 6. Unione varianti uniche (rsID unici + novel)
# ---------------------------
varianti_con_rs_uniche <- varianti_con_rs %>% distinct(dbSNP_RS, .keep_all = TRUE)
varianti_uniche_finali <- bind_rows(varianti_con_rs_uniche, varianti_novel)

# ---------------------------
# 7. Calcolo conteggi e normalizzazione per tumore
# ---------------------------
grafico_data <- varianti_uniche_finali %>%
  group_by(NEOPLASM) %>%
  summarise(Conteggio_varianti = n_distinct(SAMPLE_ID), .groups = "drop") %>%
  left_join(totali_per_tumore, by = "NEOPLASM") %>%
  mutate(Percentuale = round(100 * Conteggio_varianti / Totale_pazienti, 1))

# ---------------------------
# 8. Grafico barre normalizzato con percentuali
# ---------------------------
ggplot(grafico_data, aes(x = reorder(NEOPLASM, -Percentuale), y = Percentuale)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentuale, "%")), vjust = -0.5, size = 5) +
  labs(
    title = "Distribuzione normalizzata delle varianti per tipo di tumore (MSI > 20)",
    x = "Tipo di tumore",
    y = "Percentuale di pazienti con varianti"
  ) +
  theme_minimal(base_size = 14)

# ---------------------------
# 9. Elenco pazienti con varianti novel
# ---------------------------
pazienti_varianti_novel <- varianti_novel %>%
  select(SAMPLE_ID, NEOPLASM, Hugo_Symbol, Chromosome, Start_Position, HGVSc, HGVSp) %>%
  arrange(SAMPLE_ID)

print(pazienti_varianti_novel)