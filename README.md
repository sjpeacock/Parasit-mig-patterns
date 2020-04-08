# Parasit-mig-patterns
This repo contains code for reproducing simulations described in the paper:

- Peacock, SJ, M Krkosek, MA Lewis, PK Molnar. A unifying framework for the transient parasite dynamics of migratory hosts. PNAS #2019-08777. Submitted 23-May-2019, accepted 19-Mar-2020.

The host-macroparasite model we simulated is fully described in [Peacock et al. (2018)](https://www.sciencedirect.com/science/article/pii/S0040580916301034?via%3Dihub). We simplified that model to focus on the conditions (i.e., parameters) that lead to migratory escape, recovery, culling, and stalling.

Interested readers are encouraged to email me <stephanie.j.peacock@gmail.com> with any questions.

## Files

### Simulation code

The following files run simulations over a range of parameter values (using parallel computation) and save a .rds file of output that is used in plotting. The full results of each simulation are not stored (as that requires too much memory), but each simulation is summarized by time until parasite burdens peak, the time until parasite burdens decline below initial, peak mean parasite burden, and (for culling and stalling) the proportion of the host population that is still migrating when parasite burdens decline to initial. These summary statistics are stored as a list and saved as .rds files to be called by the `figures.R` file below. Note that the parameter names in the code correspond to the original parameterization in [Peacock et al. 2018](https://www.sciencedirect.com/science/article/pii/S0040580916301034?via%3Dihub), which differ slightly from the notation of the PNAS paper.

* `culling_time2decline.R`: Simulates the host-parasite dynamics over a range of transmission rates and rate of parasite-induced mortality, capturing migratory culling.
* `escape_time2decline.R`: Simulates the host-parasite dynamics over a range of transmission rates, host migration speed (c), and natural parasite mortality (mu_P), capturing migratory escape and recovery.
* `stalling_time2decline.R`: Simualtes the host-parasite dynamics over a range of transmission rates and rate of parasite-induced stopping, capturing migratory stalling.

### Results files
- `figures.R`: Reads .rds files produced by the simulation code and plots figures found in the manuscript.
