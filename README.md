# Balancing-oysters
Simulations and evidence to support a 'temporal balancing selection' hypothesis for larval oyster development
---
title: "Exploring how genomic architecture may present as balancing selection"
author: "Evan Durland"
date: "1/24/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##Introduction
This is a walkthrough of a simulation method to explore how contrasting genomic architectures in oysters may create changes in pooled minor allele frequency (MAF) over larval development.  

Take a scenario when you have a single genetic locus with two alleles (C/G) and it is linked to two genes with developmental importance to oysters:

....geneA....(_C_).....geneB....

....geneA....(*G*).....geneB....

Now assume that geneA becomes mutated (*x*), with detrimental effects early in development (between time 1 and 2) and it origininally showed up with the 'C' allele of the marker.  The population now has three genetic variants:

....geneA....(*C*).....geneB....

....geneA*x*....(*C*).....geneB....

....geneA....(*G*).....geneB....

in a pool, changes in the frequency of the 'C' allele reflect the death of individuals where 'C' is associated with 'geneA*x*'. Now assume a second mutation shows up for geneB which kills larvae late in development and is associated with the 'G' allele.  The population has four variants:

....geneA....(*C*).....geneB....

....geneA*x*....(*C*).....geneB....

....geneA....(*G*).....geneB....

....geneA....(*G*).....geneB*x*....

In this scenario, you can have a relative decrease in allele '*C*' early on, due to some of those mutations and then it would actually rise in frequency again as '*G*' allele individuals with the 'geneB*x*' variants die off.  




```{r}

```
