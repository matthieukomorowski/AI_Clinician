function [P, R] = mdp_example_forest (S, r1, r2, p)


% mdp_example_forest    Generate a Markov Decision Process example based on 
%                       a simple forest management 
%                       (see the related documentation for more detail)
% Arguments -------------------------------------------------------------
%   S = number of states (> 0), optional (default 3)
%   r1 = reward when forest is in the oldest state and action Wait is performed,
%        optional (default 4)
%   r2 = reward when forest is in the oldest state and action Cut is performed, 
%        optional (default 2)
%   p = probability of wild fire occurence, in ]0, 1[, optional (default 0.1)
% Evaluation -------------------------------------------------------------
%   P(SxSxA) = transition probability matrix 
%   R(SxA) = reward matrix

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
if nargin >= 1 && S <= 1
  disp('----------------------------------------------------------')
  disp('MDP Toolbox ERROR: Number of states S must be upper than 1')
  disp('----------------------------------------------------------')
elseif nargin >= 2 && r1 <= 0    
  disp('-----------------------------------------------------------')
  disp('MDP Toolbox ERROR: The reward value r1 must be upper than 0')
  disp('-----------------------------------------------------------')
elseif nargin >= 3 && r2 <=0   
  disp('-----------------------------------------------------------')
  disp('MDP Toolbox ERROR: The reward value r2 must be upper than 0')
  disp('-----------------------------------------------------------')
elseif nargin >= 4 &&( p < 0 || p > 1  )  
  disp('--------------------------------------------------------')
  disp('MDP Toolbox ERROR: Probability p must be in [0; 1]')
  disp('--------------------------------------------------------') 
else

  % initialization of optional arguments
  if nargin < 4; p=0.1; end;
  if nargin < 3; r2=2; end;
  if nargin < 2; r1=4; end;
  if nargin < 1; S=3; end;


  % Definition of Transition matrix P(:,:,1) associated to action Wait (action 1) and
  % P(:,:,2) associated to action Cut (action 2)
  %             | p 1-p 0.......0  |                  | 1 0..........0 |
  %             | .  0 1-p 0....0  |                  | . .          . |
  %  P(:,:,1) = | .  .  0  .       |  and P(:,:,2) =  | . .          . |
  %             | .  .        .    |                  | . .          . |
  %             | .  .         1-p |                  | . .          . |
  %             | p  0  0....0 1-p |                  | 1 0..........0 |        
  P1=zeros(S,S)+(1-p)*diag(ones(S-1,1),1);
  P1(:,1)=p;
  P1(S,S)=1-p;
  P2=zeros(S,S);
  P2(:,1)=1;
  P=cat(3,P1,P2);

  % Definition of Reward matrix R1 associated to action Wait and 
  % R2 associated to action Cut
  %           | 0  |                   | 0  |
  %           | .  |                   | 1  |
  %  R(:,1) = | .  |  and     R(:,2) = | .  |	
  %           | .  |                   | .  |
  %           | 0  |                   | 1  |                   
  %           | r1 |                   | r2 |
  R1=zeros(S,1);
  R1(S)=r1;
  R2=ones(S,1);
  R2(1)=0;
  R2(S)=r2;
  R=[R1 R2];

end;
