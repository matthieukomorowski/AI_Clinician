function [ bootql,bootwis] = offpolicy_multiple_eval_010518( qldata3,physpol, gamma,do_ql,iter_ql,iter_wis)
%performs all off-policy algos in one run. bootstraps to generate CIs.

if do_ql==1  %to save time, do offpol q_learning or not
[bootql]=offpolicy_eval_tdlearning( qldata3,physpol, gamma ,iter_ql);
else
bootql=55;  %gives an approximate value
end

fprintf('   Mean value of physicians'' policy by TD Learning : %f \n',nanmean(bootql));
   
[ bootwis ] = offpolicy_eval_wis( qldata3,gamma ,iter_wis);

fprintf('   Mean value of AI policy by WIS : %f \n',nanmean(bootwis));

end

