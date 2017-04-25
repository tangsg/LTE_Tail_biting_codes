%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Viterbi Encoder
%              -----------------                  
%   Author: Dama Sreekanth        
%                                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Specifications
% 1. 1/3 code rate
% 2. memory size 6
% 3. polynomials [1011011] & [1111001] [1110101]

function [p1,p2,p3] = encode_tail(ip1)

N = length(ip1);

%ip1 = rand(1,N)>0.5;
%ip2 =ip1+zeros(1,N);
ip = [circshift(ip1,[0 6]) ip1(N-5:N)];

% g = [1 0 1 1 0 1 1;1 1 1 1 0 0 1;1 1 1 0 1 0 1];
   % convolutional coding, rate - 1/3, generator polynomial - [133,171,165] octal
   cip1 = mod(conv(ip,[1 0 1 1 0 1 1]),2);
   cip2 = mod(conv(ip,[1 1 1 1 0 0 1]),2);
   cip3 = mod(conv(ip,[1 1 1 0 1 0 1]),2);
   p1 = cip1(7:N+6);
   p2 = cip2(7:N+6);
   p3 = cip3(7:N+6);
   