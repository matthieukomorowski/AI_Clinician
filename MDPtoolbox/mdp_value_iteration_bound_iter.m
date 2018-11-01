function [max_iter, cpu_time] = mdp_value_iteration_bound_iter(P, R, discount, epsilon, V0)


% mdp_value_iteration_bound_iter  Computes a bound for the number of iterations 
%                                 for the value iteration algorithm
%                                 to find an epsilon-optimal policy 
%                                 with use of span for the stopping criterion
% Arguments --------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount  = discount rate in ]0; 1[
%   epsilon   = |V - V*| < epsilon,  upper than 0,
%               optional (default : 0.01)
%   V0(S)     = starting value function, 
%               optional (default : zeros(S,1))
% Evaluation -------------------------------------------------------------
%   max_iter  = bound of the number of iterations for the value iteration algorithm
%               to find an epsilon-optimal policy with use of span for the stopping criterion
%   cpu_time  = used CPU time

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

% check of arguments
if iscell(P); S = size(P{1},1); else S = size(P,1); end;
if discount <= 0 || discount >= 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1[')
    disp('--------------------------------------------------------')
elseif nargin > 3 && (epsilon < 0)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: epsilon must be upper than 0')
    disp('--------------------------------------------------------')
elseif nargin > 4 && size(V0,1) ~= S
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: V0 must have the same dimension as P')
    disp('--------------------------------------------------------')
else
    
    if iscell(P); A = length(P); else A = size(P,3); end
    
    PR = mdp_computePR(P,R);
    
    % set default values
    if nargin < 5; V0 = zeros(S,1); end;
    if nargin < 4; epsilon = 0.01; end;
    
    % See Markov Decision Processes, M. L. Puterman, 
    % Wiley-Interscience Publication, 1994 
    % p 202, Theorem 6.6.6
    % k =    max     [1 - S min[ P(j|s,a), p(j|s',a')] ]
    %     s,a,s',a'       j
    if iscell(P)
        for ss=1:S;
            PP = [];
            for tt=1:A,
                PP = [PP;P{tt}(:,ss)];
            end;
            h(ss) = min(min(PP)); 
        end;
    else
        for ss=1:S; h(ss) = min(min(P(:,ss,:))); end;
    end
    k = 1 - sum(h);
    V1 = mdp_bellman_operator(P,PR,discount,V0);
    % p 201, Proposition 6.6.5
    max_iter = log ( (epsilon*(1-discount)/discount) / mdp_span(V1-V0) ) / log(discount*k);
end;

max_iter = ceil(max_iter);

cpu_time = cputime - cpu_time;
