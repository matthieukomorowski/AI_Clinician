function [ Q, sumQ] = OffpolicyQlearning150816( qldata3 , gamma, alpha, numtraces)
% OFF POLICY Q LEARNING

%initialisation of variables
sumQ=zeros(numtraces,1);  %record sum of Q after each iteration
nact=numel(unique(qldata3(:,3)))-1;   %nr of actions
ncl=numel(unique(qldata3(:,2)));
Q=zeros (ncl, nact);  
maxavgQ=1;
modu=100;
listi=find(qldata3(:,1)==1);   %position of 1st step of each episodes in dataset
nrepi=numel(listi);  %nr of episodes in the dataset
jj=1;

 for j=1:numtraces
    
    
    i=listi(floor(rand()*(nrepi-2))+1);  %pick one episode randomly (not the last one!)
    trace = [];
    
    while qldata3(i+1,1)~=1 
    S1=qldata3(i+1,2);
    a1=qldata3(i+1,3);
    r1=qldata3(i+1,4);
     step = [ r1, S1, a1 ];
     trace = [trace ; step];
    i=i+1;
    end

    tracelength = length(trace(:,1));
    return_t = trace(tracelength,1); % get last reward as return for penultimate state and action.
    
    for t=tracelength-1:-1:1       %Step through time-steps in reverse order
        s = trace(t,2); % get state index from trace at time t
        a = trace(t,3); % get action index
        Q(s,a) = (1-alpha)*Q(s,a) + alpha*return_t; % update Q.
        return_t = return_t*gamma + trace(t,1); % return for time t-1 in terms of return and reward at t
    end
    
     sumQ(jj,1)=sum(sum(Q));
     jj=jj+1;
     
 if mod(j,500*modu)==0  %check if can stop iterating (when no more improvement is seen)
%      sumQ(jj,1)=sum(sum(Q));
%      jj=jj+1;
     s=mean(sumQ(j-49999:j));
     d=(s-maxavgQ)/maxavgQ;
     if abs(d)<0.001
         break   %exit routine
     end
     maxavgQ=s;
 end
 

 end

 sumQ(jj:end)=[];
 
 
end

