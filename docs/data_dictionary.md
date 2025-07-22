# Data Dictionary: Seroepidemiology Dataset

This document describes the variables used in the formatted dataset files:
- [`data/DireDawa_SeroEpi24.csv`](../data/DireDawa_SeroEpi24.csv)
- [`data/AddisAbaba_Seroepi24.csv`](../data/AddisAbaba_Seroepi24.csv)


| Variable       | Description |
|----------------|-------------|
| `HH_ID`        | Household ID code (unique within site) |
| `D1.AGEX`      | Age in years |
| `D1.SEXX`      | Sex of participant (`M` = male, `F` = female) |
| `DengM`        | Dengue IgM positivity (1 = pos, 0 = neg) |
| `ZIKM`         | Zika IgM positivity  (1 = pos, 0 = neg) |
| `CHIKM`        | Chikungunya IgM positivity  (1 = pos, 0 = neg) |
| `DengG`        | Dengue IgG positivity  (1 = pos, 0 = neg) |
| `ZIKG`         | Zika IgG positivity  (1 = pos, 0 = neg) |
| `CHIKG`        | Chikungunya IgG positivity  (1 = pos, 0 = neg) |
| `RDT.DENM`     | Raw categorical RDT result for Dengue IgM ("NR" = not reactive; "R" = reactive)|
| `RDT.DENIgM`   | Raw photometric intensity for Dengue IgM |
| `RDT.ZIKM`     | Raw categorical RDT result for Zika IgM  ("NR" = not reactive; "R" = reactive)|
| `RDT.ZIKIgM`   | Raw photometric intensity for Zika IgM |
| `RDT.CHKM`     | Raw categorical RDT result for Chikungunya IgM  ("NR" = not reactive; "R" = reactive)|
| `RDT.CHKIgM`   | Raw photometric intensity for Chikungunya IgM |
| `RDT.DENG`     | Raw categorical RDT result for Dengue IgG  ("NR" = not reactive; "R" = reactive)|
| `RDT.DENIgG`   | Raw photometric intensity for Dengue IgG |
| `RDT.ZIKG`     | Raw categorical RDT result for Zika IgG  ("NR" = not reactive; "R" = reactive)|
| `RDT.ZIKIgG`   | Raw photometric intensity for Zika IgG |
| `RDT.CHKG`     | Raw categorical RDT result for Chikungunya IgG  ("NR" = not reactive; "R" = reactive)|
| `RDT.CHKIgG`   | Raw photometric intensity for Chikungunya IgG |
