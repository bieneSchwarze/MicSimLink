This repository includes a generic microsimulation model and the associated R software micSimLink for generating family structures and generational relationships. 
*micSimLink* allows to create and develop family biographies from family members for feasible population shares over calendar time. 
The crucial characteristics are sex, children, and living arrangement or partnership status. 
Further characteristics such as nationality or education level can be added if required. 
Continuous time (semi-)markov multi state models form the backbone of the model. 
For the simulation to run, an initial population and transition rates have to be provided for the characteristics (i.e., states) considered. 
The simulation itself works by applying the concept of competing risks and drawing random variables of the possible waiting times in different states. 
Construction of family linkages occurs from a female-centered perspective. 
That is, women’s fertility, partnership, dissolution and death events cause related events for linked men. 
This approach ensures consistency in the linked life histories of women and men. 
For the linkage of women and men, micSimLink includes generic matching rules that can be easily determined by the person using the approach and the software. 
Newborns are automatically linked to their mothers and, if present, also to the man to whom the mother is linked at the time of birth (who is then defined as the father).

The software itself extends the portfolio of the R package "MicSim" (https://CRAN.R-project.org/package=MicSim), which is already available for free in the R library CRAN. 
In order to deal with large and huge populations, the MicSim package includes features for running simulations in parallel on different cores. 
This repository includes a short live illustration of the approach and the software (in the file *test_MicSimLink_Exp1.R*).
