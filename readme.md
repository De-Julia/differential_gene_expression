Autor: Julia Denkewicz
# Projekt zrealizowany w ramach zaliczenia I z Bioinformatyki roślin w semestrze zimowym 2022/2023 UPWr

Informacje ogólne:
**Skrypt next-flow wykonujący analizę porównania ekspresji genów**

4 procesy: index, mapping, quantification, DESeq2

Możliwość zmiany filtrów w pliku „script_deseq2.R”

Plik wynikowy z rozszezeniem csv zapisywany w folderze "results"

* Do poprawnego działania skryptu niezbędne jest dostarczenie plików wsadowych ("22.fa", "22.gtf") w folderze "refs" oraz folder "reads" z plikami genów HBR i UHR w formacie .fq 

Wywołanie skryptu musi następować z katalogu, w którym znajduje się nadrzędny folder z plikami wsadowymi o nazwie "sample_data"
        > nextflow run skrypt.nf

# Wymagania sprzętowe:
>>>>
N E X T F L O W
      version 21.10.6 build 5660
      created 21-12-2021 16:55 UTC (17:55 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
>>>>
RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS
Mozilla/5.0 (Macintosh; Intel Mac OS X 13_2_0) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36

>>>>
featureCounts v2.0.1 (quay.io/biocontainers/subread:2.0.1--hed695b0_0)
https://quay.io/repository/biocontainers/subread?tab=tags&tag=2.0.1--hed695b0_0
Wymagane dodanie ścieżki do pliku wykonywalnego FeatureCounts do zmiennej środowiskowej PATH (do pliku .bash_profile dodać następujący wpis: 
				”export PATH=$PATH:/path/to/featureCounts”


>>>>
DESeq2 v1.38.0 (quay.io/biocontainers/bioconductor-deseq2:1.38.0--r42hc247a5b_0)
https://quay.io/repository/biocontainers/bioconductor-deseq2?tab=tags&tag=1.38.0--r42hc247a5b_0
