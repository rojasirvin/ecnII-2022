---
title: "Datos en panel, GMM y endogeneidad"
summary: " "
weight: 2
type: book
toc: false
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=T,
                      warning = F,
                      message = F)

rm(list = ls())
options(scipen=999) # Prevenir notación científica

library(tidyverse)
library(janitor)
library(sandwich)
library(clubSandwich)
library(lfe)
library(AER)
library(gmm)
library(stargazer)
```
El propósito de esta nota es visibilizar las decisiones que están implícitas en el uso de