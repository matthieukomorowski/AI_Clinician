function [V, policy, iter, cpu_time, Q] = mdp_policy_iteration_with_Q(P, R, discount, policy0, max_iter, eval_type)

% max_iter=2000;
% eval_type=1;

% mdp_policy_iteration  Resolution of discounted MDP 
%                       with policy iteration algorithm
% Arguments ---------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA) = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount = discount rate, in ]0, 1[
%   policy0(S) = starting policy, optional 
%   max_iter = maximum number of iteration to be done, upper than 0, 
%              optional (default 1000)
%   eval_type = type of function used to evaluate policy: 
%              0 for mdp_eval_policy_matrix, else mdp_eval_policy_iterative
%              optional (default 0)
% Evaluation --------------------------------------------------------------
%   V(S)   = value function 
%   policy(S) = optimal policy
%   iter     = number of done iterations
%   cpu_time = used CPU time
%--------------------------------------------------------------------------
% In verbose mode, at each iteration, displays the number 
% of differents actions between policy n-1 and n

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


cpu_time = cputime;

global mdp_VERBOSE;

% check of arguments
if iscell(P); S = size(P{1},1); else S = size(P,1); end;
if iscell(P); A = length(P); else A = size(P,3); end

if discount <= 0 || discount >= 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1[')
    disp('--------------------------------------------------------')
elseif nargin > 3 && (size(policy0,1)~=S || any(mod(policy0,1)) || any(policy0<1)|| any(policy0>A) )
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: policy0 must a (Sx1) vector with integer from 1 to A')
    disp('--------------------------------------------------------')
elseif nargin > 4 && max_iter <= 0
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
    disp('--------------------------------------------------------')
else
    
    PR = mdp_computePR(P,R);

    % initialization of optional arguments
    if nargin < 6; eval_type = 0; end;
    if nargin < 5; max_iter = 1000; end;
    if nargin < 4;
        % initialization of policy: 
        % the one wich maximizes the expected immediate reward
        [nil, Q, policy0] = mdp_bellman_operator_with_Q(P,PR,discount,zeros(S,1));
    end;
    
%      if mdp_VERBOSE; disp('  Iteration  Number_of_different_actions'); end;
%         disp('  Iteration  Number_of_different_actions');

    iter = 0;
    policy = policy0;
    is_done = false;
    while ~is_done
        iter = iter + 1;
        if  (eval_type==0)   
            V = mdp_eval_policy_matrix(P,PR,discount,policy);         
        else
            V = mdp_eval_policy_iterative(P,PR,discount,policy);
        end;
        [nil, Q, policy_next] = mdp_bellman_operator_with_Q(P,PR,discount,V);
        
	n_different = sum(policy_next ~= policy);
%          if mdp_VERBOSE; disp(['       ' num2str(iter) '                 '  num2str(n_different)]); end;
%      disp(['       ' num2str(iter) '                 '  num2str(n_different)]);

        if all(policy_next==policy) || iter == max_iter || (iter > 20 && n_different <=5) || iter > 50
            is_done = true; 
        else
            policy = policy_next;
        end;
    end;
    
end; 

cpu_time = cputime - cpu_time;
