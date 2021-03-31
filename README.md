# Package "rondo"

Analytical calculation of groundwater flow in a single aquifer under circular shaped polders.

Based on the publication "Berekening van de stationaire grondwaterstroming van hogere naar lagere (ronde) gebieden".
C.H. van Immerzeel, Internationale Agrarische Hogeschool Larenstein Velp (1990).


## Installation

You can install the released version of menyanthes from with:

`install_github("KeesVanImmerzeel/rondo")`

Then load the package with:

`library("rondo")` 

## Functions in this package
- `rd_init()`: Initialise a list (y) with parameters used for the calculations of the groundwater heads and fluxes.
- `rd_phi()`: Calculate the head at radius x.
- `rd_q()`: Calculate the lateral discharge at radius x.
- `rd_seep()`: Calculate the seepage intensity at radius x.
- `av_seep()`: Calculate the average seepage intensity between radius r1 and r2.

## Get help

To get help on the functions in this package type a question mark before the function name, like `?rd_init`

## Example usage

`kD <- c(1000, 1000, 1000)`

`c <- c(1000, 2000, 3000)`

`r <- c(1000, 2000)`

`h <- c(10, 9, 8)`

`y <- rd_init(kD, c, r, h)`

`x <- seq(0, 3000, by = 100)`

`plot(x, rd_phi(x,y))`

# References

1. *Steady flow of ground water towards wells*
Comm. voor Hydrologisch Onderzoek T.N.O. Verslagen en mededelingen No. 10 (1964).

2. *Wegzijging en kwel* 
De grondwaterstroming van hogere naar lagere gebieden. L.F. Ernst, I.C.W.-rapport nr. 7 (1983).
