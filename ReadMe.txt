Boddepalli et al., "Integrative Modelling of Innate Immune Response Dynamics during Virus Infection"

This repository contains MATLAB codes of the Viral-Immune Model and different analyses that were done. 

Cite the article upon use of the codes.

The function of different files is listed below.

Core Model Files

ODEs.m – Defines the system of ODEs (model equations).
Param.xlsx – Contains information of parameter values.
Param_vector.m – Generates a parameter vector used in model simulations (unless specified otherwise).
RunThis.m – Runs the model simulation and plots system dynamics (corresponds to Figure 1b in the manuscript).

Bifurcation_plotting.m – Varies immune activation levels and computes threshold values of viral antagonism, generating the bifurcation front (Figure 2b).

Latin_Hypercube_Sampling folder contains the following set of codes:
-ParamSets_LHS.m – Generates 25,000 parameter sets using Latin Hypercube Sampling and saves them as bounds.mat.
-RunLHS.m – Runs simulations across generated parameter sets and stores steady-state values in session.mat.
-LHS_plotting.m – Analyzes viral load distributions and parameter correlations (SI Figure S3, Figure 2c).

Run_BPSA: Run_BPSA.m – Performs bifurcation point sensitivity analysis (Figure 2d and SI Figure S4). 

Run_DiffRegimes: Simulates system behavior under different viral antagonism regimes and parameter perturbations (0.5× and 2×) for the following parameters- k_{t, ISG RNA}, mu_{ISG RNA}, k_{r, V}, k_{c, V}, k_{t, V} (Figure 2e, 2f and SI Figures S5-S6).

Clustering_Run: Clusters immune response dynamics and visualizes patterns (Figure 3 and SI Figures S8-10).

IFN_SequentialDosing: Simulates sequential IFN dosing in a JAK–STAT signaling model (Figure 4).

IFN_effects:  Simulates pre-, co-, and post-infection IFN treatment scenarios for different viral antagonism strengths (Figure 5a and 5c).

IFN_Antiviral_Combn_Run: This script simulates the combined effect of IFN pretreatment and antiviral by varying parameter (e.g. k_c) and IFN doses (Figure 5c and SI Figure S12). 







