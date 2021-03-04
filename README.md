# cibn
__Causal independence Bayesian networks__

`cibn` is an R package implementing elicitation, estimation and inference functionalities for Bayesian networks under the causal independence assumption, i.e., non-interacting parent variables, called _causal independence Bayesian networks_. They allow three exaustive types of variables (graded, double-graded and multi-valued nominal variables) and admit the Causal Independence Decomposition (CID), which increases efficiency of elicitation, estimation and inference. Causal interactions can be added upon need.
The reference paper is:

A. Magrini (2021). Efficient decomposition of Bayesian networks with non-graded variables. _International Journal of Statistics and Probability_, 10(2): 52-67. https://www.doi.org/10.5539/ijsp.v10n2p52


R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `cibn` package. R can be downloaded from https://www.r-project.org/.

To install the `cibn` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/cibn")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

In this example, we create a simple Bayesian network to infer risk attitude of bank customers.
The variables in the network are:
- `Age`: age in years, double-graded variable with sample space: `(18_30, 31_50, 51_)`;
- `Edu`: education level, double-graded variable with sample space: `(primary_or_less, secondary, tertiary)`;
- `Marital`: marital status, graded variable with sample space: `(single, convivent)`;
- `Parent`: parentship, graded variable with sample space: `(no, yes)`;
- `Risk`: risk attitude, double-graded variable with sample space: `(low, normal, high)`;
- `Portf`: type of portfolio, double-graded variable with sample space: `(money_market, mixed, stock_market)`;
- `Life`: life insurance, multi-valued nominal variable with sample space: `(long_term, short_term, none)`.

The hypothesized Directed Acyclic Graph (DAG) is shown ![here](https://github.com/alessandromagrini/cibn/blob/main/bin/bankrisk_dag.pdf).
Also, we hypothesize a causal interaction between marital status and parentship in determining risk attitude.

Here is the model code with prior parameter values.
Probability distributions are automatically normalized, thus the user can specify them using counts (null counts are interpreted as a uniform distribution).
```
bankrisk_code <- '

  variable Age {
    type = DGRAD
    states = (18_30, 31_50, 51_)
    parents = ()
    description = <Age>
    }

  model Age {
    omitted = (1,3,2)
    }

  variable Edu {
    type = DGRAD
    states = (primary_or_less, secondary, tertiary)
    parents = ()
    description = <Education level>
    }

  model Edu {
    omitted = (1,7,5)
    }

  variable Marital {
    type = GRAD
    states = (single, convivent)
    parents = (Age)
    description = <Single or convivent>
    }

  model Marital {
    omitted = (2,3)
    Age:18_30 = (3,1)
    Age:51_ = (1,3)
    }

  variable Parent {
    type = GRAD
    states = (no, yes)
    parents = (Age)
    description = <Parentship>
    }

  model Parent {
    omitted = (3,1)
    Age:18_30 = (4,1)
    Age:51_ = (2,1)
    }

  variable Risk {
    type = DGRAD
    states = (low, normal, high)
    parents = (Age, Edu, Marital_Parent)
    description = <Risk attitude>
    }

  interaction Marital_Parent {
    from = (Marital, Parent)
    description = <Interaction between marital status and parentship>
    }

  model Risk {
    omitted = (1,3,2)
    Age:18_30 = (2,1,5)
    Age:51_ = (5,1,2)
    Edu:primary_or_less = (3,2,1)
    Edu:tertiary = (1,3,4)
    Marital_Parent:convivent+no = (1,1,2)
    Marital_Parent:single+yes = (2,1,1)
    Marital_Parent:convivent+yes = (3,1,0)
    }

  variable Portf {
    type = DGRAD
    states = (money_market, mixed, stock_market)
    parents = (Risk)
    description = <Type of portfolio>
    }

  model Portf {
    omitted = (3,5,2)
    Risk:low = (3,1,0)
    Risk:high = (0,1,3)
    }

  variable Life {
    type = NOM
    states = (long_term, short_term, none)
    parents = (Risk)
    description = <Life insurance>
    }

  model Life {
    omitted = (2,4,1)
    Risk:low = (0,1,2)
    Risk:high = (2,1,0)
    }
  '
```
Creation of the network
```
bankrisk_bn <- new.cibn(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code, maximal=FALSE)  ## disable maximal CID
```
Graphic of the DAG
```
plot(bankrisk_bn, attrs=list(edge=list(arrowsize=0.5)))
plot(bankrisk_bn, attrs=list(edge=list(arrowsize=0.5)), full=TRUE)  ## full DAG
```
Exact inference (interface to gRain package)
```
# see the sample spaces
getStates(bankrisk_bn)
# risk attitude for a customer aged 31-50 years and with a mixed portfolio
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="mixed"))
# risk attitude for a customer aged 31-50 years and with a money market portfolio
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="money_market"))
# risk attitude for a customer aged 31-50 years and with a stock market portfolio
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="stock_market"))
```
Random sample from the network
```
bankrisk_sam <- sample.cibn(bankrisk_bn, nsam=100)
head(bankrisk_sam)
summary(bankrisk_sam)
