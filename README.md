Companion Github Repository for the paper ["The molecular impact of cigarette smoking resembles aging across tissues"](https://www.biorxiv.org/content/10.1101/2024.03.14.585016v1 )

# Getting the code
```
git clone https://github.com/Mele-Lab/2024_GTEx_Smoking.git
```

# Software dependencies

| Package | Version | Packages | Version |
| -------- | ------- | ------- |------- |
| edgeR | 4.0.16 | dplyr | 1.1.4 |
| optparse | 1.7.5 | limma | 3.58.1 |
| missMethyl | 1.38.0 | ggplot2 | 3.5.0 |
| reshape2 | 1.4.4 | caret | 6.0-94 |
| hier.part | 1.0-6 | sandwich | 3.1-0 |
| multcomp | 1.4-25 | lmtest | 0.9-40 |


# Running the code

Code to run differential gene expression analysis
```
/analysis/Expression_models.sh
```

Code to run differential splicing analysis
```
/analysis/Splicing_models.sh
```

Code to run differential methylation analysis
```
/analysis/Methylation_models.sh
```
Code to run enrichement analsyis on gene expresion results 

```
/analysis/Enrichement_analysis.sh
```

All expected output is available as supplementary tables in the paper
