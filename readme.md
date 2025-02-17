# Simulation model
This model can be applied in studies of tree species with different characteristics, from tropical and temperate forests, to assess population persistence in restoration sites. 
 
## author
Patricia Sanae Sujii

## description
This is an individual-based model that allows to evaluate the effect of
different levels of initial genetic diversity on the short and
mid-terms’ population viability.

## input
The input necessary is a tabular delimited table with individual
location (x, y), genotypes, and age. There is no limitation for
number of individuals and the size of the restoration area. See example file (datain.txt).

## output
The outpu will be:
* files with individual location (x, y), genotypes, age, plant parents and survival probability, in two formats (this program's format and FSTAT's format).
* summary genetic statistics
* example output files listed here are named as: number of mother trees_number of pollen donors_area size. Example: 1m_1f_5ha.zip are the result files of simulations of a restoraion area created with seeds from one mother tree, one pollen donor in an area of 5ha.
