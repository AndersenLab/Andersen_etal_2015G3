Data
====

This directory contains all of the raw cytometer data from *insert citation here*. These data were collected from a Union Biometrica BIOSORT large-particle flow cytometer in August of 2011. The data represent 2 experimental condition, control and paraquat. The data are organized using the following scheme:

```
Data
├── 20110811_RIAILs0a
│   ├── ctrlstrains.Rds
│   ├── pqstrains.Rds
│   ├── reports
│   │   └── All report files
│   ├── score
│   │   └── Score data files
│   └── setup
│       └── Setup data files
├── 20110812_RIAILs0b
│   ├── ctrlstrains.Rds
│   ├── pqstrains.Rds
│   ├── reports
│   │   └── All report files
│   ├── score
│   │   └── Score data files
│   └── setup
│       └── Setup data files
├── 20110818_RIAILs0c
│   ├── ctrlstrains.Rds
│   ├── pqstrains.Rds
│   ├── reports
│   │   └── All report files
│   ├── score
│   │   └── Score data files
│   └── setup
│       └── Setup data files
├── 20110819_RIAILs0d
│   ├── ctrlstrains.Rds
│   ├── pqstrains.Rds
│   ├── reports
│   │   └── All report files
│   ├── score
│   │   └── Score data files
│   └── setup
│       └── Setup data files
├── README.md
└── RIAILs0_ProcessedPhenotypes.csv
```

The names of each data file can be broken apart to reveal information about the assay. For example the first data folder can be used to identify:

+ **20110811** - The date of the assay
+ **RIAILs** - The type of assay (using Recombinant Inbred Advanced Intercross Lines or RIAILs)
+ **0** - The round number, this was the first of set of high-throughput assays conducted
+ **a** - The assay number, this was the first day on which plates were run through the cytometer