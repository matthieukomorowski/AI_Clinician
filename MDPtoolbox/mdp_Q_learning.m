function [Q, V, policy, mean_discrepancy] = mdp_Q_learning(P, R, discount, N)

% mdp_Q_learning   Evaluation of the matrix Q, using the Q learning algorithm 
%
% Arguments
% -------------------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS)
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount  = discount rate in ]0; 1[
%   N(optional) = number of iterations to execute, default value: 10000.
%                 It is an integer greater than the default value. 
% Evaluation --------------------------------------------------------------
%   Q(SxA) = learned Q matrix 
%   V(S)   = learned value function.
%   policy(S) = learned optimal policy.
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


% check of arguments
if (discount <= 0 || discount >= 1)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0,1[')
    disp('--------------------------------------------------------')   
elseif (nargin >= 4) && (N < 10000)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: N must be upper than 10000')
    disp('--------------------------------------------------------') 
else

    % initialization of optional arguments
    if (nargin < 4); N=10000; end;      
    
    % Find number of states and actions
    if iscell(P)
        S = size(P{1},1);
        A = length(P);
    else
        S = size(P,1);
        A = size(R,2); 
    end;
    
    % Initialisations
    Q = zeros(S,A);
    dQ = zeros(S,A);
    mean_discrepancy = [];
    discrepancy = [];

    % Initial state choice
    s = randi([1,S]);

    h = waitbar(0,'Initializing waitbar...');
    
    for n=1:N

           
           
        % Reinitialisation of trajectories every 100 transitions
        if (mod(n,100)==0); s = randi([1,S]); waitbar(n/N,h,n/N*100); end;
        
        % Action choice : greedy with increasing probability
        % probability 1-(1/log(n+2)) can be changed
        
        
        
        
        
        
        pn = rand(1);
        if (pn < (1-(1/log(n+2))))
          [nil,a] = max(Q(s,:));
        else
          a = randi([1,A]);
        end;
 
        % Simulating next state s_new and reward associated to <s,s_new,a> 
        p_s_new = rand(1);
        p = 0; 
        s_new = 0;
        while ((p < p_s_new) && (s_new < S)) 
            s_new = s_new+1;
            if iscell(P)
                p = p + P{a}(s,s_new);
            else   
                p = p + P(s,s_new,a);
            end;
        end; 
        if iscell(R)
            r = R{a}(s,s_new); 
        elseif ndims(R) == 3
            r = R(s,s_new,a); 
        else
            r = R(s,a); 
        end;

        % Updating the value of Q   
        % Decaying update coefficient (1/sqrt(n+2)) can be changed
        delta = r + discount*max(Q(s_new,:)) - Q(s,a);
        dQ = (1/sqrt(n+2))*delta;
        Q(s,a) = Q(s,a) + dQ;
    
        % Current state is updated
        s = s_new;
 
        % Computing and saving maximal values of the Q variation  
        discrepancy(mod(n,100)+1) = abs(dQ);  
    
        % Computing means all over maximal Q variations values  
        if (length(discrepancy) == 100)     
           mean_discrepancy = [ mean_discrepancy mean(discrepancy)];
           discrepancy = [];
        end;   
    
    end;

    %compute the value function and the policy
    [V, policy] = max(Q,[],2);        
close(h);
end;
