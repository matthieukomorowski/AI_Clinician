# AI Clinician: reinforcement learning in intensive care


Code for a reinforcement learning model applied to the management of intravenous fluids and vasopressors in patients with sepsis in intensive care.

Related to publication: https://www.nature.com/articles/s41591-018-0213-5

Author: Dr Matthieu Komorowski, Imperial College London, 2016-2018 - m.komorowski14@imperial.ac.uk

The 2 datasets used in the research are:
- MIMIC-III : https://mimic.physionet.org/
- eICU-RI: not publicly available in full, subset available here: http://eicu-crd.mit.edu/

Cohort definition: all adult patients fulfilling the sepsis-3 definition: http://jamanetwork.com/journals/jama/fullarticle/2492881
The unique identifiers for these patients in both datasets are provided, along with a detailed desciption of the datasets.

This repository contains Matlab code to:
1.	build 500 different discrete state and action MDP models from the MIMIC-III training dataset;
2.	select the best policy from off-policy evaluation on the MIMIC-III validation set;
3.	test this optimal policy on the eICU-RI dataset;
4.	compute the main results and key figures.

External files and toolboxes used:

- MDP toolbox (some functions modified):  https://uk.mathworks.com/matlabcentral/fileexchange/25786-markov-decision-processes--mdp--toolbox
- Fastknnsearch : https://uk.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit
