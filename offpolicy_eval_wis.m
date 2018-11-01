function [ bootwis ,c,individual_trial_estimators] = offpolicy_eval_wis( qldata3,gamma ,num_iter)
% WIS estimator of AI policy
% Thanks to Omer Gottesman (Harvard) for his assistance

bootwis=cell(num_iter,1);
p=unique(qldata3(:,8));
prop=25000/numel(p); %25000 patients of the samples are used
prop=min([prop 0.75]);  %max possible value is 0.75 (75% of the samples are used)


for jj=1:num_iter
    
ii=floor(rand(size(p,1),1)+prop); % prop% of the samples are used
j=ismember(qldata3(:,8),p(ii==1));
q=qldata3(j==1,:);
fence_posts=find(q(:,1)==1);
num_of_trials=size(fence_posts,1);
individual_trial_estimators = NaN(num_of_trials,1);
rho_array=NaN(num_of_trials,1);
 c=0;  %count of matching pairs pi_e + pi_b
     
for i=1:num_of_trials-1
      rho=1;
      for t=fence_posts(i):fence_posts(i+1)-2 %stops at -2
            rho=rho*   q(t,6)/q(t,5);  
      end
      if rho>0
             c=c+1;
      end
      rho_array(i)=rho;    
end

ii=isinf(rho_array)|isnan(rho_array);  %some rhos are INF
normalization=nansum(rho_array(~ii));

 for i=1:num_of_trials-1

        current_trial_estimator = 0;
        rho = 1;
        discount = 1/gamma   ;
            
            for t=fence_posts(i):fence_posts(i+1)-2 %stops at -2 otherwise ratio zeroed
               
                   rho=rho*   q(t,6)/q(t,5);        
                   discount =discount* gamma;
                   current_trial_estimator =current_trial_estimator+ discount * q(t+1,4);


            end
            
      individual_trial_estimators(i) =  current_trial_estimator*rho;
     
 end
 
    bootwis(jj) = {nansum( individual_trial_estimators(~ii) )/normalization};
    
end

individual_trial_estimators = individual_trial_estimators(~ii)./rho_array(~ii);  %for the last iteration only!

bootwis=cell2mat(bootwis);

end

