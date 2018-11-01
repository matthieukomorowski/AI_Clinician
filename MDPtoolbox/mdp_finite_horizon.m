function [V, policy, cpu_time] = mdp_finite_horizon(P, R, discount, N, h)


% mdp_finite_horizon   Reolution of finite-horizon MDP with backwards induction 
% Arguments -------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA) = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount = discount factor, in ]0, 1]
%   N        = number of periods, upper than 0
%   h(S)     = terminal reward, optional (default [0; 0; ... 0] )
% Evaluation -------------------------------------------------------------
%   V(S,N+1)     = optimal value function
%                  V(:,n) = optimal value function at stage n
%                         with stage in 1, ..., N
%                         V(:,N+1) = value function for terminal stage 
%   policy(S,N)  = optimal policy
%                  policy(:,n) = optimal policy at stage n
%                         with stage in 1, ...,N
%                         policy(:,N) = policy for stage N
%   cpu_time = used CPU time
%--------------------------------------------------------------------------
% In verbose mode, displays the current stage and policy transpose.

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
if N < 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: N must be upper than 0')
    disp('--------------------------------------------------------')
elseif discount <= 0 || discount > 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
    disp('--------------------------------------------------------')
else
    
    if iscell(P); S = size(P{1},1); else S = size(P,1); end;
    
    V = zeros(S,N+1);
    if nargin == 5; V(:,N+1) = h; end;
    
    PR = mdp_computePR(P,R);
    
    for n=0:N-1
        [W,X]=mdp_bellman_operator(P,PR,discount,V(:,N-n+1));
        V(:,N-n)=W; 
        policy(:,N-n) = X;
        if mdp_VERBOSE
            disp(['stage:' num2str(N-n) '      policy transpose : ' num2str(policy(:,N-n)')]); 
        end;
    end
    
end;

cpu_time = cputime - cpu_time;

