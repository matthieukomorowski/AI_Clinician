function [V, policy, iter, cpu_time] = mdp_policy_iteration_modified(P, R, discount, epsilon, max_iter)


% mdp_policy_iteration_modified    Resolution of discounted MDP  
%                                  with modified policy iteration algorithm
% Arguments -------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA) = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount = discount rate in ]0, 1]
%              beware to check conditions of convergence for discount = 1.
%   epsilon  = epsilon-optimal policy search, upper than 0,
%              optional (default : 0.01)
%   max_iter = maximum number of iteration to be done in the inner loop,
%              upper than 0, optional (default: 10)
% Evaluation -------------------------------------------------------------
%   V(S)     = value function
%   policy(S)= epsilon-optimal policy
%   iter     = number of main iterations
%   cpu_time = used CPU time
%--------------------------------------------------------------------------
% In verbose mode, at each iteration, displays the variation of V

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
if discount <= 0 || discount > 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
    disp('--------------------------------------------------------')
elseif nargin > 4 && epsilon <= 0
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: epsilon must be upper than 0')
    disp('--------------------------------------------------------')
elseif nargin > 5 && max_iter <= 0
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
    disp('--------------------------------------------------------')
else

    if discount == 1  
        disp('-------------------------------------------------------')
        disp('MDP Toolbox WARNING: check conditions of convergence.')
        disp('With no discount, convergence is not always assumed.')
        disp('--------------------------------------------------------')
    end;
    
    if iscell(P); S = size(P{1},1); else S = size(P,1); end;
    
    PR = mdp_computePR(P,R);

    % initialization of optional arguments
    if nargin < 5; max_iter = 10; end;
    if nargin < 4; epsilon = 0.01; end;

    % computation of threshold of variation for V for an epsilon-optimal policy
    if discount ~= 1
        thresh = epsilon * (1-discount)/discount;
    else 
        thresh = epsilon;
    end;

    if discount == 1
        V = zeros(S,1);
    else
        V = 1/(1-discount)*min(min(PR))*ones(S,1);
    end;     

%    if mdp_VERBOSE; disp('  Iteration  V_variation'); end;
        disp('  Iteration  V_variation');
        
    iter = 0;
    is_done = false;
    while ~is_done

        iter = iter + 1;
        
        [Vnext, policy] = mdp_bellman_operator(P,PR,discount,V);
        %[Ppolicy, PRpolicy] = mdp_computePpolicyPRpolicy(P, PR, policy);
        
        variation = mdp_span(Vnext - V);
%         if mdp_VERBOSE; 
%              disp(['      ' num2str(iter,'%5i') '         ' num2str(variation)]); 
%         end;
          
             disp(['      ' num2str(iter,'%5i') '         ' num2str(variation)]); 
        

        V=Vnext;
        if variation < thresh
            is_done = true; 
        else
	    is_verbose = false;
            if mdp_VERBOSE; mdp_VERBOSE = 0; is_verbose = true; end;
            V = mdp_eval_policy_iterative(P, PR, discount, policy, V, epsilon, max_iter);
            if is_verbose; mdp_VERBOSE = 1; end;
        end;
    end;
end;

cpu_time = cputime - cpu_time;
