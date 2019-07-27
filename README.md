# OrzekMasterThesis

# Installation
If you want to install OrzekMasterThesis from GitHub, use the following commands in R:

```{r, eval=FALSE}
if(!require(devtools))install.packages("devtools")
devtools::install_github("jhorzek/OrzekMasterThesis")
```

# Running the simulations:
## Simulation 1:

```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(regsem)
library(tcltk)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
simulation1(sampleSize = sampleSize, seed = seed, wd = wd, total_repetitions = total_repetitions, zeroThresh = .001)
```

## Simulation 2:
```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(regsem)
library(tcltk)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
gamma = 0
simulation2(sampleSize = sampleSize, seed = seed, wd = wd, total_repetitions = total_repetitions, gamma = gamma)
```

## Simulation 3:
### 10 percent non-zero:
```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(tcltk)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
crossEffect = .2
autoEffect = .5
simulation3_10(sampleSize, seed, wd, total_repetitions, crossEffect, autoEffect)
```

### 25 percent non-zero:
```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(tcltk)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
crossEffect = .2
autoEffect = .5
simulation3_25(sampleSize = sampleSize, seed = seed, wd = wd, total_repetitions = total_repetitions, crossEffect = crossEffect, autoEffect = autoEffect)
```

### 50 percent non-zero:
```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(tcltk)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
crossEffect = .2
autoEffect = .5
simulation3_50(sampleSize = sampleSize, seed = seed, wd = wd, total_repetitions = total_repetitions, crossEffect = crossEffect, autoEffect = autoEffect)
```

## Simulation 4:
```{r, eval=FALSE}
library(OrzekMasterThesis)
library(OpenMx)
library(caret)
library(laremm)
library(tcltk)
library(ctsem)
sampleSize = 200 # set the sample size
seed = 12341324 # set a seed
wd = "~/Desktop" # set a working directory. The outcome file will be saved here
total_repetitions = 1000 # number of repetitions
simulation4(sampleSize = sampleSize, seed = seed, wd = wd, total_repetitions = total_repetitions, Scale = T)
```
