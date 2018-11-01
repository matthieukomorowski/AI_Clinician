function Vpolicy = mdp_eval_policy_iterative(P, R, discount, policy, V0, epsilon, max_iter)


% mdp_eval_policy_iterative   Policy evaluation using iteration. 
% Arguments -------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%               P could be an array with 3 dimensions or 
%               a cell array (1xS), each cell containing a matrix possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount  = discount rate in ]0; 1[
%   policy(S) = a policy
%   V0(S)     = starting value function, optional (default : zeros(S,1))
%   epsilon   = epsilon-optimal policy search, upper than 0,
%               optional (default : 0.0001)
%   max_iter  = maximum number of iteration to be done, upper than 0, 
%               optional (default : 10000)
% Evaluation -------------------------------------------------------------
%   Vpolicy(S) = value function, associated to a specific policy
%--------------------------------------------------------------------------
% In verbose mode, at each iteration, displays the condition which stopped iterations:
% epsilon-optimum value function found or maximum number of iterations reached.

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


global mdp_VERBOSE;
 
% check of arguments
if iscell(P); S = size(P{1},1); else S = size(P,1); end;
if iscell(P); A = length(P); else A = size(P,3); end
if discount <= 0 || discount >= 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0,1[')
    disp('--------------------------------------------------------')
elseif size(policy,1)~=S || any(mod(policy,1)) || any(policy<1) || any(policy>A)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: policy must be a (1xS) vector with integer from 1 to A')
    disp('--------------------------------------------------------')
elseif nargin > 4  && size(V0,1) ~= S
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: V0 must have the same dimension as P')
    disp('--------------------------------------------------------')
elseif (nargin > 5) && (epsilon <= 0)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: epsilon must be upper than 0')
    disp('--------------------------------------------------------')
elseif (nargin > 6) && (max_iter <= 0)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: max_iter must be upper than 0')
    disp('--------------------------------------------------------')
else

    % initialization of optional arguments 
    if nargin < 5; V0 = zeros(S,1); end;
    if nargin < 6; epsilon = 0.0001; end;
    if nargin < 7; max_iter = 10000; end; 


    [Ppolicy, PRpolicy] = mdp_computePpolicyPRpolicy(P, R, policy);
    if mdp_VERBOSE; disp('  Iteration    V_variation'); end;
    iter = 0;
    Vpolicy = V0;
    is_done = false; 
    while ~is_done
        iter = iter + 1;
        Vprev = Vpolicy;     
        Vpolicy = PRpolicy + discount * Ppolicy * Vprev;
        variation = max(abs(Vpolicy - Vprev));
        if mdp_VERBOSE;disp(['      ' num2str(iter,'%5i') '         ' num2str(variation)]); end;
        if variation < ((1-discount)/discount)*epsilon % to ensure |Vn - Vpolicy| < epsilon
            is_done = true; 
            if mdp_VERBOSE; disp('MDP Toolbox: iterations stopped, epsilon-optimal value function'); end;
        elseif iter == max_iter
            is_done = true; 
            if mdp_VERBOSE; disp('MDP Toolbox: iterations stopped by maximum number of iteration condition'); end;
        end;
    end;
    
end;


