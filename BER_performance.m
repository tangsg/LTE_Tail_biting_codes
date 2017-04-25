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

% semilogy(Ebnodb,ppp(1,:),Ebnodb,ppp(2,:),Ebnodb,ppp(3,:),Ebnodb,ppp(4,:),Ebnodb,ppp(5,:),Ebnodb,ppp(6,:),Ebnodb,erfc(10.^(Ebnodb/20))/2);
%  H = legend('iter 1','iter 2','iter 3','iter 5','iter 10','iter 18','uncoded');
% %grid on
% hold on
% H = legend('iter 1');






%% Trash

% first encoder
% y1 = turbo_enco(x);
% %interleaver
% i = 0:N-1;
% f1=3;f2=10;
% p = mod(f1*i+f2*(i.^2),N)+1;
% x1 = x(p);
% 
% % second encoder
% y2 = turbo_enco(x1);
% 
% code = reshape([x;y1;y2],1,[]);
% %% viterbi decoding
% z = turbo_viterbi(1-2*y1);
% 
% length(find(z-x))


%% Turbo decoder

%z1 = turbo_decoder(x,y1,y2);
