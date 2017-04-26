% /*
% ML viterbi algorithm for tail biting codes
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
function [ipHat_v ] = viterbi_ML(y)

% Viterbi decoding
 LUT = [zeros(64,32) ones(64,32)];
ref_main = [0 0 0;1 1 1;1 0 0;0 1 1;0 0 1;1 1 0;1 0 1;0 1 0]; 
ref_1 = [ref_main;double(xor(ref_main,kron([1 1 0],ones(size(ref_main,1),1))))];
ref_2 = [ref_1;double(xor(ref_1,kron([1 1 1],ones(size(ref_1,1),1))))];
ref_3 = [ref_2;double(xor(ref_2,kron([0 1 1],ones(size(ref_2,1),1))))];
ref_last1 = [ref_3;double(xor(ref_3,kron([1 1 1],ones(size(ref_3,1),1))))];

state_ML = zeros(1,64);
r =zeros(1,3,size(y,1)); % taking 3 coded bits
rv = zeros(128,3,size(y,1));
euclidianDist1 = zeros(128,1,size(y,1));
%ref_soft=1/sqrt(2)*(ones(size(ref_last))-2*ref_last);
   state_trellis = zeros(64,length(y)/3+1);
   path_surv = zeros(64,length(y)/3,64);
   for tt=1:64  %run viterbi 64 times
   for ii = 1:size(y,2)/3
       for kk = 1:size(y,1)   % 16 copies of the message
      r(1,:,kk) = y(kk,3*ii-2:3*ii); % taking 3 coded bits
      ref_soft=1/sqrt(2)*(ones(size(ref_last1))-2*ref_last1);
      
      % computing the Euclidian distance distance between ip coded sequence with [00x;x01;x10;x11]
      rv(:,:,kk) = kron(ones(128,1),r(1,:,kk));
      euclidianDist1(:,:,kk) = sum((rv(:,:,kk).*ref_soft(:,:,kk)),2);
      euclidianDist = sum(euclidianDist1,3);
        end
            if (ii<=6)      %for first level 1
       
       for jj = 1:32
           % state 1-32
          bm = euclidianDist(2*jj-mod(ceil(tt/(2^(ii-1))),2))-euclidianDist(64+2*jj-mod(ceil(tt/(2^(ii-1))),2))+state_trellis(2*jj-mod(ceil(tt/(2^(ii-1))),2),ii);
          state_trellis(jj,ii+1) = bm;
          path_surv(jj,ii,tt) = 2*jj-mod(ceil(tt/(2^(ii-1))),2);
       end
       
       for jj = 1:32    
           % state 33-64
           bm1 = euclidianDist(64+2*jj-mod(ceil(tt/(2^(ii-1))),2))-euclidianDist(2*jj-mod(ceil(tt/(2^(ii-1))),2))+state_trellis(2*jj-mod(ceil(tt/(2^(ii-1))),2),ii);
          state_trellis(32+jj,ii+1)= bm1;
          path_surv(32+jj,ii,tt) = 2*jj-mod(ceil(tt/(2^(ii-1))),2);
       end
                
            else 
      for jj = 1:32
                    % state 1-32
          bm1 = euclidianDist(2*jj-1)-euclidianDist(64+2*jj-1)+state_trellis(2*jj-1,ii);
          bm2 = euclidianDist(2*jj)-euclidianDist(64+2*jj)+state_trellis(2*jj,ii);
          [state_trellis(jj,ii+1), idx] = max([bm1,bm2]);
          path_surv(jj,ii,tt) = idx+2*(jj-1);
      end
                    % states 33-64
      for jj = 1:32
          bm1 = euclidianDist(64+2*jj-1)-euclidianDist(2*jj-1)+state_trellis(2*jj-1,ii);
          bm2 = euclidianDist(64+2*jj)-euclidianDist(2*jj)+state_trellis(2*jj,ii);
          [state_trellis(32+jj,ii+1), idx] = max([bm1,bm2]);
          path_surv(32+jj,ii,tt) = idx+2*(jj-1);
      end
            end
   end
   state_ML(tt) = state_trellis(tt,length(y)/3+1);
   end
    % trace back unit
    [~, currState1] = max(state_ML);
   
    currState = currState1;
   ipHat_v = zeros(1,length(y)/3);
   for jj = length(y)/3:-1:1
      prevState =  path_surv(currState,jj,currState1); 
      ipHat_v(jj) = int32(LUT(prevState,currState));
      currState = prevState;
   end
