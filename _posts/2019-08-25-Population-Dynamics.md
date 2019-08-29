---
title: "Genetic Dynamics in Populations"
date: 2019-08-25
tags: [biology, evolution, excel]
layout: splash
header:
  overlay_image: "/images/abstracthead.jpg"
excerpt: "Modeling gene frequency in populations"
mathjax: true
---

Population genetics is the study of genetic composition in a group and its change over time. In this post I'll cover the conventional methods of modeling population genetics, starting with a simple model and refining it by adding terms representing the effects of additional parameters. I've chosen to use Excel for this purpose since it's easier to visualize the processes, but the same principles can easily be utilized in R or Python as well.

Since this post is intended to be self-contained, I'll introduce/recap some biology info before building up to the model itself.

All organisms have DNA containing many genes, individual regions that code for some functional product or trait. Each gene contains loci (plural of locus), sites at which 1 of 2+ alternate alleles can exist. Only one allele can be present within each strand of DNA, and correspondingly, double-stranded DNA organisms can have 2 total alleles.

Alleles are usually notated by their dominance, the interaction they have with other alleles. Dominant alleles have their contents expressed if at least one copy is present, while recessive alleles are only expressed if two copies are present. I'll denote dominant alleles with an 'A' and recessive alleles with 'a'. A genotype refers to the set of alleles an organism has, with homozygote indicating that both alleles are the same while heterozygote denoting that they're different.

Finally, I'll use $$p$$ and $$q$$ as the percentage of dominant and recessive alleles respectively. With the preliminary material out of the way we can begin building up our model.

#### Hardy-Weinberg Equilibrium

Using initial conditions for $$p$$, $$q$$, and population, we can calculate the starting population and its genotype breakdown. Since $$p$$ and $$q$$ are the percentage of alleles of the total, we can use them as probabilities and treat each allele as an independent event: an allele has probability $$p$$ to be dominant and $$q$$ to be recessive, where $$p = 1 - q$$. From this we can see that there is $$p^2$$ probability to be AA, and $$q^2$$ probability to be aa. Since there are two orientations for heterozygotes (Aa and aA), the respective probability is $$2pq$$.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene1.JPG" alt="Genotype Counts">

We can then use these probabilities, and multiply them by the initial condition to find the number of individuals with each genotype at $$t = 0$$.

We can also work backwards from genotype counts to calculate $$p$$ and $$q$$ by finding the frequency of each allele. $$p = (AA + .5Aa)/Population$$, and $$q = (aa +.5Aa)/Population$$.

To incorporate birth and death rates, we can apply them to ending populations to determine starting population in the next period. Side note: this process has autoregressive properties, the output of the time series is related to the output at t - 1. Birth and death rates can be arithmetic $$f(t+1) = f(t) +- h$$ or geometric $$f(t+1) = f(t) /* u$$. Arithmetic rates are level independent, while geometric rates can change based on the current level. Since population births and deaths depend on the population size, it makes sense to use geometric rates for this. We'll incorporate this feature into the model itself, but I'll use rate inputs of births/deaths per 1000.

Example: 100 birth rate, 90 death rate, population of 2000. The birth rate results in a change in population of $$100*(2000/1000) = 200$$, while the death rate results in a change in population of $$-90*(2000/1000) = 180$$.

We can combine these into a net term: $$(BR - DR) * (Population/1000)$$. Adding this term to the previous total population gives the total population for the start of the next period. Using the previous period's $$p$$ and $$q$$ probabilities we can calculate counts for each genotype using the new population and the independent probabilities of alleles. By repeating this process we can project the population counts for any number of periods.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene2.JPG" alt="HW Equilibrium">

By playing around with the birth and death rate inputs, we can observe that the population growth or decay varies, eventually becoming 0 or growing in an exponential fashion. But the allele frequency remains constant through the prediction period. We can also adjust the starting allele frequencies, and observe the same effect. $$p'$$ and $$q'$$ remain the same as the initial conditions in this case. This special case is called [Hardy Weinberg Equilibrium](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2475721/), as described in simple algebraic terms by geneticists G. H Hardy and Wilhelm Weinberg. This equilibrium predicts that given the following assumptions, allele frequencies will remain the same regardless of initial conditions.

1. No selection
 - Individuals with different genotypes have the same chance of survival and reproduction
2. No mutation
 - Individual's genotypes remain permanent until death
3. No migration
 - The population of interest doesn't gain or lose individuals
4. Large populations
 - This assumption is required due to the non-uniform nature of birth and death rates in reality. All individuals born or death can be distributed in any way across genotype counts. In small populations this can have the effect of large variations in actual net change per genotype count. In an application of the Law of large numbers, as the population size increases, the expected change per genotype count approaches the rate times the percent of the total population the genotype count represents.
5. Random mating
 - Random mating ensures that variation is maintained. Without this assumption, we could imagine a case where certain genotypes don't mate with others, which would change $$p$$ and $$q$$.

This model was a departure from previous rationale that more common alleles would invariably spread and become fixed. Hardy Weinberg Equilibrium represents the null model of population dynamics, and can be refined to incorporate other input parameters.

#### Simple Selection
Simple selection incorporates a term for percentage survival for each genotype. Different genotypes can vary in their survivability due to physiological effects. Additionally, even though heterozygotes express the same traits as dominant homozygotes, they can have differing survival rates due to underlying effects of allele interaction:

  - An example of heterozygotes having differential survival rates is the allele for human sickle cell anemia. In areas with high malaria-carrying mosquito populations, heterozygotes have higher survival rates than recessive or dominant homozygotes. The sickle cell allele impacts the shape of red blood cells in the body, which impacts their ability to clot and carry nutrients. The dominant wild allele codes for a functional red blood cell, while the recessive mutant allele codes for a misshapen cell (in the shape of a hook or sickle). Dominant homozygotes have only normal red blood cells, and are susceptible to malaria. Recessive homozygotes have only sickled red blood cells, and are susceptible to the physiological effects of sickle cell anemia, but are resistant to contracting malaria. Heterozygotes have some red blood cells that are normal and some that are sickled, and do not experience the full negative effects of sickle cell anemia, but also maintain the resistance to malaria. As a result, heterozygotes have higher survival rates than homozygotes.

We can incorporate this easily by multiplying it by the output of the HW model. Multiplying the genotype counts by the corresponding survival rate gives the post birth-death population counts (in this example I've treated the genetic drift birth and death terms as "normal" effects due to life expectancy of the individuals in the population, and the selection rates as deaths attributed to external factors such as predation).

<img src="{{ site.url }}{{ site.baseurl }}/images/gene4.JPG" alt="Simple Selection">

#### Complex Selection
The parameters for simple selection can be difficult to estimate, making the previous model difficult to use in practice. The selection against genotypes is often immeasurable: genotype information requires DNA analysis, which is often difficult or impossible to perform until the individual has died. As a result we can use complex selection parameters, in which selection against alleles is estimated instead of selection against genotypes. Selection against alleles can be estimated by performing test crosses and observing lethality due to differences at the loci of interest.

Replacing the survival rate term in the previous model with $$(1 - selection)$$ where selection represents the strength of adverse selection for the relevant allele, scaled from 0 to 1. Heterozygotes are selected against by the dominant selection term, since they express the dominant phenotype.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene5.JPG" alt="Complex Selection">

By playing around with these selection parameters, we can see that the allele frequencies can move to one extreme or remain in balance at a level corresponding to the initial conditions.

#### Genetic Drift Without Death
In all the previous models, we have applied the net change consistently across the population of interest. Given a population of 1000 AA, 1000 Aa, 1000 aa with birth rate of 100 per 1000, these models would allocate 100 births per genotype. Obviously this expectation seems rational based on the  average outcome, but we should only expect this to be the actual outcome in large populations over long time periods. In shorter time frames and in smaller populations the randomness of the rates dominates the outcome. For example, births and deaths could be distributed unevenly across genotype counts, causing the allele frequencies to change even without selection. The randomness of the net rate can force one allele to increase until it becomes "fixed" (frequency of 1).

We can model this fairly easily by using the *RAND()* function in excel. *RAND()* returns a number from a uniform distribution between 0 and 1, which can be combined with a probability term to represent the change to birth. We want the birth rates for each genotype to sum to 1, where each birth rate represents the proportion of the total period births we can attribute to the genotype. The number of births each genotype count receives should be proportional to the previous period count which can be modeled as shown below.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene6.JPG" alt="Stochastic Birth">

By dividing each birth term by the sum of current birth terms, we can arrive at the stochastic proportion of births attributed to each genotype. This stochastic term is nonuniformly distributed between 0 and 1. The distribution is no longer uniform like the *RAND()* function due to the level dependent rate terms we introduced, and the shape of this distribution changes with initial conditions and over time.

By adding the stochastic term times the birth rate times the previous period's ending total population divided by 1000, we can find the stochastic birth change per genotype. Adding this term to the previous total population times the complex selection term times the corresponding allele probabilities gives us the starting genotype population for the next period.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene7.JPG" alt="Genetic Drift" type = center>

##### Genetic Drift with Death
The previous model ignored the stochastic death term, which will now be introduced. The same concept as the birth term can be used to model the death term in stochastic terms.

The death and birth terms must be modeled using separate *RAND()* functions, since using the same *RAND()* outcome would cause the birth and death to become dependent. We can rationalize that the number of births and deaths that a genotype group experiences are independent events if they occur at the same moment (otherwise we would have to adjust the population based on the event that occurred first before calculating the level-dependent rate of the second event). Modeling using this method allows any combination of birth and death rate distribution across genotype counts.

Adding this random death term to the previous formula requires a bit more work to implement. Since the death terms can result in the population dropping to or below zero, we have to include a few simple *IF()* statements to terminate the model in case this happens.

The intuition is fairly simple: If the count is less than or equal to 0, output 0, else output the count. The stochastic terms can also be changed to be 0 if the corresponding genotype count is 0; there's no need to calculate the births or deaths if the population is all dead. Lastly, to get the 0's to carry through the model we just wrap the previous model term with an *IF()* statement that produces a 0 if the $$t - 1$$ population count was 0. The corresponding formulas are shown below.

<img src="{{ site.url }}{{ site.baseurl }}/images/gene8.JPG" alt="Full Genetic Drift">

##### Recap
I've built up to a simple model for modeling and simulating population genetics, including stochastic birth and death terms, simple and complex selection terms, and initial frequency parameters. This model requires nontrivial assumptions about the underlying behavior of the genes in question:

1. 2 alleles per locus. This model could be expanded to include more than 2 alleles, though it would dramatically increase the complexity. Each additional allele introduces more possible genotypes, and terms to the calculation of allele frequency. Ex: 2 alleles produce 3 genotypes AA Aa aa, while 3 alleles (A, B, and C for clarity) produce 6 genotypes AA BB CC AB AC BC. The number of possible genotypes for any number of alleles equals $$n(n+1)/2$$.

2. Random segregation of alleles. The model assumes that an allele at one position is independent of the opposing position. In reality this isn't always true: some alleles exhibit some degree of linkage or are likely to occur together.

3. Death and birth timing. The model assumes that death and birth occur at the same time, and occur instantaneously. In reality birth and death occur randomly in continuous time. The issue of discrete time is a major limitation of the model developed here.

The workbook I've made is available [*here*](https://github.com/arun-krishnaraj/arun-krishnaraj.github.io/blob/master/uploads/Population%20Dynamics%20v%201.1.xlsx), play around with the input parameters and observe how the model predictions respond.

###### A note on discrete time
We've treated the population as stationary during the periods between the modeled periods. If we treat the time intervals as days, then our rates should correspondingly be daily births/deaths per 1000. By reducing the time interval we can approach a continuous process. I'll leave the transition to continuous time for a later post, in which I'll move away from Excel for development of a more robust model.
