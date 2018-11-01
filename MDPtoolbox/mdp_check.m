function error_msg = mdp_check(P , R)

% mdp_check    Check if the matrices P and R define a Markov Decision Process
% Let S = number of states, A = number of actions
%      The transition matrix P must be on the form P(SxSxA) 
%      and P(:,:,a) must be stochastic
%      The reward matrix R must be on the form (SxSxA) or (SxA)
% Arguments ---------------------------------------------------------------
%   P(SxSxA) = transition matrix
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
% Evaluation --------------------------------------------------------------
%   error_msg = error message or '' if correct

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


is_error_detected = false;
error_msg = '';

% check of P
if iscell(P)
    s1 = size(P{1},1);
    s2 = size(P{1},2);
    a1 = length(P);
else
    [s1 s2 a1] = size(P);
end

if s1 < 1 || a1 < 1 || s1 ~= s2 
    error_msg = 'MDP Toolbox ERROR: The transition matrix must be on the form P(S,S,A) with S : number of states greater than 0 and A : number of action greater than 0';
    is_error_detected = true;
end;

if ~is_error_detected
    a = 1;
    while a <= a1
        if iscell(P)
            error_msg = mdp_check_square_stochastic ( P{a} );
        else 
            error_msg = mdp_check_square_stochastic ( P(:,:,a) );
        end
        if isempty(error_msg)
	    a = a + 1;
        else
	    a = a1 + 1;
        end
    end;
end;


% check of R
if ~is_error_detected
    if iscell(R)
        s3 = size(R{1},1);
        s4 = size(R{1},2);
        a2 = length(R);
    elseif issparse(R)
        s3 = size(R,1);
        s4 = s3;
        a2 = size(R,2);
    elseif ndims(R) == 3
        [s3 s4 a2] = size(R);
    else
        [s3 a2] = size(R);
        s4 = s3;
    end
    if s3 < 1 || a2 < 1 || s3 ~= s4
        error_msg = 'MDP Toolbox ERROR: The reward matrix R must be an array (S,S,A) or (SxA) with S : number of states greater than 0 and A : number of actions greater than 0';
        is_error_detected = true;
    end;
end;

if ~is_error_detected
     if s1~=s3 || a1 ~= a2
        error_msg = 'MDP Toolbox ERROR: Incompatibility between P and R dimensions';
	    is_error_detected = true;
     end;
end;


