% /*
% Test encoder and decoder
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
%%test
clear all
N=10^2;
x = rand(1,N)>0.5;
 x = x+0;
[p1,p2,p3] = encode_tail(x);

y = [p1;p2;p3];
y = y(:).';
y = 1-2*y;
z = viterbi_ML(y);
 length(find(x-z))
 
 z = viterbi_wava(y,4);
 length(find(x-z))