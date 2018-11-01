function [V,mean_discrepancy]=mdp_eval_policy_TD_0(P,R,discount,policy,N)

% mdp_eval_policy_TD_0  Evaluation of the value function, using the TD(0) algorithm 
%
% Arguments
% -------------------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount  = discount rate in ]0; 1[
%   policy(S) = optimal policy
%   N(optional) = number of iterations to execute, default value: 10000.
%              It is an integer greater than the default value. 
% Evaluation --------------------------------------------------------------
%   V(S)   = optimal value function.
%   mean_discrepancy(N/100) = vector of V discrepancy mean over 100 iterations
%             Then the length of this vector for the default value of N is 100.

% MDPtoolbox: Markov Decision Processes Toolbox
% Copyright (C) 2009  INRA
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright notice, 
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%      this list of conditions and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the <ORGANIZATION> nor the names of its contributors 
%      may be used to endorse or promote products derived from this software 
%      without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.


% Arguments checking
if iscell(P); S = size(P{1},1); else S = size(P,1); end;
if iscell(P); A = length(P); else A = size(P,3); end
if (nargin < 4 || nargin > 5)    
    disp('------------------------------------------')
    disp('The number of arguments must be 4 or 5')
    disp('------------------------------------------')       
elseif (discount <= 0 || discount >= 1)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0,1[')
    disp('--------------------------------------------------------')   
elseif size(policy,1)~=S || any(mod(policy,1)) || any(policy<1) || any(policy>A)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: policy must be a (1xS) vector with integer from 1 to A')
    disp('--------------------------------------------------------')
elseif (nargin == 5) && (N < 10000)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: N must be upper than 10000')
    disp('--------------------------------------------------------') 
else
    
    % initialization of optional argument
    if nargin < 5; N = 10000; end;
        
    
    % Initializations
    s=randi([1,S]);      % Initial state choice
    V=zeros(S,1);
    mean_discrepancy=[]; % vector of mean V variations
    discrepancy=[];      % vector of V variation

    for n=1:N

        % Reinitialisation of trajectories every 100 transitions
	if (mod(n,100)==0); s = randi([1,S]); end;
        
        % Select an action: here action of the policy
        a = policy(s);
    
        % Simulate next state s_new 
        p_s_new = rand(1);
        p = 0; 
        s_new = 0;
        while ((p <= p_s_new) && (s_new < S)) 
            s_new = s_new+1;   
            if iscell(P)
                p = p+P{a}(s,s_new);
            else
                p = p+P(s,s_new,a);
            end
        end; 
    
        % Update V
        if iscell(R)
            r = R{a}(s,s_new);
        else
            if ndims(R) == 3
                r = R(s,s_new,a);
            else
                r = R(s,a);
            end
        end
        delta = r + discount*V(s_new)-V(s);
        dV = (1/sqrt(n+2))*delta;
        V(s) = V(s)+dV;
    
        % Update current state
        s = s_new;
     
        % Memorize mean V variations on each trajectory (100 transitions) 
        discrepancy(mod(n,100)+1) = abs(dV);  
        if (length(discrepancy)==100)     
            mean_discrepancy = [mean_discrepancy mean(discrepancy)];
            discrepancy = [];
        end;   
    end;

end;
