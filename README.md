# AI Clinician: reinforcement learning in intensive care


Code for a reinforcement learning model applied to the management of intravenous fluids and vasopressors in patients with sepsis in intensive care.

Related to publication: https://www.nature.com/articles/s41591-018-0213-5

Author: Dr Matthieu Komorowski, Imperial College London, 2015-2019 - m.komorowski14@imperial.ac.uk

The 2 datasets used in the research are:
- MIMIC-III : https://mimic.physionet.org/
- eICU-RI: not publicly available in full, subset available here: http://eicu-crd.mit.edu/

Cohort definition: all adult patients fulfilling the sepsis-3 definition: http://jamanetwork.com/journals/jama/fullarticle/2492881


This repository contains:
I. the Jupyter notebook to perform data extraction in MIMIC-III
II. the Matlab code to identify the cohort of patients with sepsis in MIMIC-III
III. the Matlab code to re-create the MIMIC-III dataset
IV. the Matlab code (AIClinician_core_111118.m) to:
1.	build 500 different discrete state and action MDP models from the MIMIC-III training dataset;
2.	select the best policy from off-policy evaluation on the MIMIC-III validation set;
3.	test this optimal policy on the eICU-RI dataset;
4.	compute the main results and key figures.
V. Additional files:
1. The unique identifiers for these patients in both datasets are provided (patientIDs_MIMIC3.csv and patientIDs_eRI.csv). Note: you'll need to add 200,000 to all the patient identifiers in patientIDs_MIMIC3 to match the numbering found in the initiail database.
2. A detailed desciption of the datasets (Dataset description Komorowski 111118.xlsx). 


External files and toolboxes used:

- MDP toolbox (some functions modified):  https://uk.mathworks.com/matlabcentral/fileexchange/25786-markov-decision-processes--mdp--toolbox
- Fastknnsearch : https://uk.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit
