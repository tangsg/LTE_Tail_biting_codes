% /*
% BER performance of tail biting codes
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

clear all
%coder -build max_log_map.prj;
%matlabpool 24;
iterations = [4];
ttt=200;
%Ebnodb =[0.75:.25:2];%
%Ebnodb =[6:2:8];%
Ebnodb =[1:1:4];
Ecnodb =Ebnodb-10*log10(3); %[1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.25 2.5 2.75 3 3.5 4 4.5 5];
N=100; 
ppp=zeros(length(iterations),length(Ecnodb));
ppp_fer=zeros(length(iterations),length(Ecnodb));

ppp1=zeros(1,length(Ecnodb));
ppp_fer1=zeros(1,length(Ecnodb));
tic
for iter=1:length(iterations)
%tic
iter
it = iterations(iter);
for SNR=1:length(Ecnodb)
%x = [1 0 1 0 1 1 1 0 0 0];
SNR
ppp(iter,SNR)=0;
for kkk=1:ttt
 
 x = rand(1,N)>0.5;
 x = x+0;
% encoding
[xk,zk,zk1]=encode_tail(x);

% add noise

data = reshape([xk;zk;zk1],1,[]);
mod_code = (1-2*data);
% mod_code = QPSK_mod(data).';
tx = mod_code+((10^(-Ecnodb(SNR)/20))*(randn(1,length(mod_code))+1i*randn(1,length(mod_code))))/sqrt(2);
% rx = demod_QPSK_soft(tx);
rx = real(tx);
[z]= viterbi_wava(rx,it);
errr = length(find(x-z));
ppp(iter,SNR)=ppp(iter,SNR)+errr;
ppp_fer(iter,SNR)=ppp_fer(iter,SNR)+ceil(errr/N);
%% ML
z1 = viterbi_ML(rx);
errr = length(find(x-z1));
ppp1(SNR)=ppp1(SNR)+errr;
ppp_fer1(SNR)=ppp_fer1(SNR)+ceil(errr/N);
end
end
ppp_fer(iter,:) = ppp_fer(iter,:)/(ttt);
ppp(iter,:) = ppp(iter,:)/(N*ttt);
%toc
end
toc
time_taken=toc;
%save maxlog_1000
%semilogy(Ebnodb,ppp(1,:),Ebnodb,ppp(2,:),Ebnodb,ppp(3,:),0:1:10,erfc(10.^((0:1:10)/20))/2);%,Ebnodb,ppp(2,:));%,Ebnodb,ppp(3,:));%,Ebnodb,ppp(4,:),Ebnodb,ppp(5,:),Ebnodb,erfc(10.^(Ebnodb/20))/2);
%semilogy(Ebnodb,ppp(1,:));
semilogy(Ebnodb,ppp(1,:),Ebnodb,ppp1/(N*ttt),0:1:10,erfc(10.^((0:1:10)/20))/2);
hold all
%semilogy(Ebnodb,ppp_fer(1,:),Ebnodb,ppp_fer(2,:),Ebnodb,ppp_fer(3,:),Ebnodb,ppp_fer(4,:),Ebnodb,ppp_fer(5,:),Ebnodb,erfc(10.^(Ebnodb/20))/2);
% semilogy(Ebnodb,ppp_fer(1,:),Ebnodb,ppp_fer(2,:),Ebnodb,ppp_fer(3,:),Ebnodb,erfc(10.^(Ebnodb/20))/2);
% grid on
H = legend('1','4','8','uncoded theory');%,'iter 10','iter 20','uncoded');
% hold on

