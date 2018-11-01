function [Ppolicy, PRpolicy] = mdp_computePpolicyPRpolicy(P, R, policy)


% mdp_computePpolicyPRpolicy  Computes the transition matrix and the reward matrix for a policy
% Arguments --------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   policy(S) = a policy
% Evaluation -------------------------------------------------------------
%   Ppolicy(SxS)  = transition matrix for policy
%   PRpolicy(S)   = reward matrix for policy

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


if iscell(P); A = length(P); else A = size(P,3); end

for a=1:A % avoid looping over S
    
    ind = find(policy == a); % the rows that use action a
    if ~isempty(ind)
        if iscell(P)
            Ppolicy(ind,:) = P{a}(ind,:);
        else
            Ppolicy(ind,:) = P(ind,:,a);
        end
        PR = mdp_computePR(P,R);
        PRpolicy(ind,1) = PR(ind,a);
    end
end

if issparse(PR)
    PRpolicy = sparse(PRpolicy);
end
   
