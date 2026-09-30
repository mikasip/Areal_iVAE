## Identifiable Variational Autoencoders for Blind Source Separation of Multivariate Areal Spatio-Temporal Data
This repository contains the code and supplementary materials for reproducing the results presented in manuscript "Identifiable Variational Autoencoders for Blind Source Separation of Multivariate Areal Spatio-Temporal Data".

## Repository Structure
```
├── helpers/                                # Helper files for generating data and auxiliary variables
├── simulations/                            # Files for running simulations and plotting the results
├── analysis_tycho.R                        # A file to reproduce the results of the case study
├── simulation_figures.R                    # Code to reproduce some of the figures in the simulation study
├── tycho_aggregated.RData                  # The preprocessed dataset
└── README.md                               # This file
```

## Dependencies

The main methods used in the analysis are in R package ```NonlinearBSS```, which can be installed by running:

```
devtools::install_github("mikasip/NonlinearBSS")
```

## Data

The case study data are publicly available from Project Tycho (https://www.tycho.pitt.edu). The data was constructed from 6 separate datasets, one for each of the considered diseases [1-6].

## References

[1] Van Panhuis W., Cross A., Burke D., Counts of Measles reported in UNITED STATES OF AMERICA: 1888-2002: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.14189004

[2] Van Panhuis W., Cross A., Burke D., Counts of Diphtheria reported in UNITED STATES OF AMERICA: 1888-1981: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.397428000

[3] Van Panhuis W., Cross A., Burke D., Counts of Pertussis reported in UNITED STATES OF AMERICA: 1888-2017: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.27836007

[4] Van Panhuis W., Cross A., Burke D., Counts of Scarlet fever reported in UNITED STATES OF AMERICA: 1888-1969: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.30242009

[5] Van Panhuis W., Cross A., Burke D., Counts of Smallpox reported in UNITED STATES OF AMERICA: 1888-1952: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.67924001

[6] Van Panhuis W., Cross A., Burke D., Counts of Typhoid fever reported in UNITED STATES OF AMERICA: 1888-2005: (version 2.0, April 1, 2018): Project Tycho data release, DOI: 10.25337/T7/ptycho.v2.0/US.4834000 


