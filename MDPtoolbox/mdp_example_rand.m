function [P, R] = mdp_example_rand (S, A, is_sparse, mask)


% mdp_example_rand    Generate a random Markov Decision Process
% Arguments -------------------------------------------------------------
%   S = number of states (> 0)
%   A = number of actions (> 0)
%   is_sparse = false to have matrices in plain format, true to have sparse matrices
%               optional (default false).
%   mask(SxS) = matrix with 0 and 1 (0 indicates a place for a zero probability), 
%               optional (default, ones(S,S) )
% Evaluation -------------------------------------------------------------
%   P(SxSxA) = transition probability matrix 
%   R(SxSxA) = reward matrix

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



% arguments checking
if S < 1 || A < 1
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: The number of states S ')
    disp('and the number of actions A must be upper than 1.')
    disp('--------------------------------------------------------')
elseif nargin >= 4 && ( size(mask,1) ~= S || size(mask,2) ~= S )
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: mask must be a SxS matrix') 
    disp('--------------------------------------------------------')
else
    
    % initialization of optional arguments
    if nargin < 3; is_sparse = false; end;
    if nargin < 4; mask = ones(S,S); end;
    
    if is_sparse
        % definition of transition matrix : square stochastic matrix
        P = {};
        for a=1:A
            PP = sparse(mask .* rand(S));
            for s=1:S
                PP(s,:) = PP(s,:) / sum( PP(s,:) );
            end;
            P{a} = PP;
        end;
        
        % definition of reward matrix (values between -1 and +1)
        R = {};
        for a=1:A
	    R{a} = sparse(mask .* ( 2*rand(S) - ones(S,S) ));
        end;
    else
        % definition of transition matrix : square stochastic matrix
        for a=1:A
            P(:,:,a) = mask .* rand(S);
            for s=1:S
                P(s,:,a) = P(s,:,a) / sum( P(s,:,a) );
            end;
        end;
        
        % definition of reward matrix (values between -1 and +1)
        for a=1:A
	    R(:,:,a) = mask .* ( 2*rand(S) - ones(S,S) );
        end;
    end
    
end;



     
