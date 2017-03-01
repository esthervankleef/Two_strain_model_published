# Two strain model 

This model simulates the transmission dynamics of two competing strains, i.e. a hospital-adapted resistant strain and a community-adapted sensitive strain. 

It has been argued that since infection control measures, such as hand hygiene, should affect resistant and sensitive strains equally, 
observed discordant changes must have largely resulted from other factors, such as changes in antibiotic use. 
This model is used to test the validity of this reasoning.

Model_functions.R contains all functions including the model, plotting functions and functions to calculate R0.

Run_model.R contains the code to run the model and generate output

plot_output.R contains the code to plot the model output. 

Reference: Van Kleef E, Luangasanatip, N, Bonten MJM, Cooper BS. Why susceptible bacteria are resistant to  hospital infection control. Wellcome Open Res 2017
