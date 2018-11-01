function Vpolicy = mdp_eval_policy_matrix(P, R, discount, policy)


% mdp_eval_policy_matrix    Evaluation of the value function of a policy
% Arguments -------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA) = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount = discount rate in ]0; 1[
%   policy(S) = a policy
% Evaluation -------------------------------------------------------------
%   Vpolicy(S) = value function of the policy

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


% check of arguments
if iscell(P); S = size(P{1},1); else S = size(P,1); end;
if iscell(P); A = length(P); else A = size(P,3); end
if discount <= 0 || discount >= 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1[')
    disp('--------------------------------------------------------')
elseif size(policy,1)~=S || any(mod(policy,1)) || any(policy<1) || any(policy>A)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: policy must be a (1xS) vector with integer from 1 to A')
    disp('--------------------------------------------------------')
else
        
    [Ppolicy, PRpolicy] = mdp_computePpolicyPRpolicy(P, R, policy);
 
    % V = PR + gPV  => (I-gP)V = PR  => V = inv(I-gP)*PR
    Vpolicy = (speye(S,S) - discount*Ppolicy) \ PRpolicy;

end
