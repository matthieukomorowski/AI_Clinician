function [policy, iter, cpu_time] = mdp_value_iterationGS(P, R, discount, epsilon, max_iter, V0)


% mdp_value_iterationGS   Resolution of discounted MDP with value iteration Gauss-Seidel algorithm 
% Arguments --------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount  = discount rate in ]0; 1]
%               beware to check conditions of convergence for discount = 1.
%   epsilon   = epsilon-optimal policy search, upper than 0,
%               optional (default : 0.01)
%   max_iter  = maximum number of iteration to be done, upper than 0, 
%               optional (default : computed)
%   V0(S)     = starting value function, optional (default : zeros(S,1))
% Evaluation -------------------------------------------------------------
%   policy(S) = epsilon-optimal policy
%   iter      = number of done iterations
%   cpu_time  = used CPU time
%--------------------------------------------------------------------------
% In verbose mode, at each iteration, displays the variation of V
% and the condition which stopped iterations: epsilon-optimum policy found
% or maximum number of iterations reached.

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
if discount <= 0 || discount > 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
    disp('--------------------------------------------------------')
elseif nargin > 3 && (epsilon <= 0)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: epsilon must be upper than 0')
    disp('--------------------------------------------------------')
elseif nargin > 4 && max_iter <= 0
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
    disp('--------------------------------------------------------')
elseif nargin > 5 && size(V0,1) ~= S
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: V0 must have the same dimension as P')
    disp('--------------------------------------------------------')    
else

     if discount == 1  
         disp('--------------------------------------------------------')
         disp('MDP Toolbox WARNING: check conditions of convergence.')
         disp('With no discount, convergence is not always assumed.')
         disp('--------------------------------------------------------')
     end;
     
     if iscell(P); A = length(P); else A = size(P,3); end
     
     PR = mdp_computePR(P,R);
     
     % initialization of optional arguments
     if nargin < 6; V0 = zeros(S,1); end;
     if nargin < 4; epsilon = 0.01; end;
     % compute a bound for the number of iterations
     if discount ~= 1
        computed_max_iter = mdp_value_iteration_bound_iter(P, R, discount, epsilon, V0);
     end;
     if nargin < 5
         if discount ~= 1
             max_iter = computed_max_iter;
         else
             max_iter = 1000;
         end;
     else
         if discount ~= 1 && max_iter > computed_max_iter
             disp(['MDP Toolbox WARNING: max_iter is bounded by ' num2str(computed_max_iter,'%12.1f') ])
             max_iter = computed_max_iter;
         end;
     end;
     
     % computation of threshold of variation for V for an epsilon-optimal policy
     if discount ~= 1
         thresh = epsilon * (1-discount)/discount;
     else 
         thresh = epsilon;
     end;
     
     iter = 0;
     V = V0;
     is_done = false;
     if mdp_VERBOSE; disp('  Iteration  V_variation'); end;
     while ~is_done
         iter = iter + 1;
         Vprev = V;
         
         for s = 1:S
             for a = 1:A
                 if iscell(P)
                     Q(a) =  PR(s,a)  +  discount * P{a}(s,:) * V; 
                 else
                     Q(a) =  PR(s,a)  +  discount * P(s,:,a) * V; 
                 end              
             end;
             V(s) = max(Q);
         end;
         
         variation = mdp_span(V - Vprev);
         if mdp_VERBOSE; 
             disp(['      ' num2str(iter,'%5i') '         ' num2str(variation)]); 
         end;
         if variation < thresh 
              is_done = true; 
              if mdp_VERBOSE 
                  disp('MDP Toolbox : iterations stopped, epsilon-optimal policy found')
              end;
         elseif iter == max_iter
              is_done = true; 
              if mdp_VERBOSE 
                  disp('MDP Toolbox : iterations stopped by maximum number of iteration condition')
              end;
         end;
    end;
    
    for s = 1:S
        for a = 1:A
            if iscell(P)
                Q(a) =  PR(s,a) + P{a}(s,:) * discount * V; 
            else
                Q(a) =  PR(s,a) + P(s,:,a) * discount * V; 
            end            
        end;
        [V(s) , policy(s,1)] = max(Q);
    end;

end;

cpu_time = cputime - cpu_time;
