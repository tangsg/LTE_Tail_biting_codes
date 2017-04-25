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