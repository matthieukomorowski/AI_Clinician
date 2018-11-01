# AI Clinician: reinforcement learning in intensive care


Code for a reinforcement learning model applied to the management of intravenous fluids and vasopressors in patients with sepsis in intensive care.


Author: Dr Matthieu Komorowski, Imperial College London, 2016-2018 - m.komorowski14@imperial.ac.uk


This repository contains PostgreSQL and Matlab code to:

1.	define cohorts of patients fulfilling the sepsis-3 definition in two databases: MIMIC-III (https://mimic.physionet.org/) and eICU-RI (not publicly available in full, subset available here: http://eicu-crd.mit.edu/)
2.	extract the data of interest from both databases
3.	build 500 different discrete state and action MDP models from the MIMIC-III training dataset
4.	select the best policy from off-policy evaluation
5.	test this optimal policy on the eICU-RI dataset
6.	compute the main results and key figures

External files and toolboxes used:

MDP toolbox (some functions modified):  https://uk.mathworks.com/matlabcentral/fileexchange/25786-markov-decision-processes--mdp--toolbox

Fastknnsearch : https://uk.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit

References for the sepsis-3 definition: Singer, JAMA 2016 http://jamanetwork.com/journals/jama/fullarticle/2492881
