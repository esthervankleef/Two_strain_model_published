# Why susceptible bacteria are resistant to  hospital infection control

This release contains additional sensitivity analysis showing the impact of changing the degree of bacterial inference and mixing.

This model simulates the transmission dynamics of two competing strains, i.e. a hospital-adapted resistant strain and a community-adapted sensitive strain, in the hospital and its catchment area.

It has been argued that since infection control measures, such as hand hygiene, should affect resistant and sensitive strains equally, observed discordant changes must have largely resulted from other factors, such as changes in antibiotic use. This model is used to test the validity of this reasoning.

model_functions.R contains all functions including the model, plotting functions and functions to calculate R0.

run_model.R contains the code to run the model and generate output. Parameters listed are those at baseline.

plot_output.R contains the code to plot the model output.

run_model_cdi.R contains the code to run the model and generate output for the Lancet ID correspondence to Dingle et al 2017 (DOI: http://dx.doi.org/10.1016/S1473-3099(16)30514-X). Parameters listed are those used to generate the plot.

plot_model_cdi.R contains the code to generate the figure in the Lancet ID correspondence.
