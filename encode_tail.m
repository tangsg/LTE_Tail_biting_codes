% /*
% Tail biting convolutional encoder
%     Copyright (C) 2017  sreekanth
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  * Author: sreekanth dama
%  * Contact: sreekanth@iith.ac.in
%  **/

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
   