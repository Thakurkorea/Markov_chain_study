%% ���-�ð迭�� ���� ���縶���� �� �����ϱ�
% ����Ȯ��(TRANS)�� ���Ȯ��(EMIS)�� �˰� �ִ� ��쿡 
% ���� ��ǵ��� �ð迭�� �����ϴ� ���:
% 2���� ������°� �ְ� 6���� �������� ����� �ִ� ���:
tic
clc;
clear 
close all;
%% Read the sequence file

input_data=xlsread('Hmm_llto.csv');

%% Get the individual sequence 
no_of_states=2;
init_TRANS=[0.86 0.14;...
    0.04 0.96];
% init_EMIS=zz;
%% Fixed init_EMIS
init_EMIS=[0.8 0.2
    0.3 0.7]
for i=2:8;
    seq1=input_data(:,2);
%% Training the sequence
[estTR,estE] = hmmtrain(seq1',init_TRANS,init_EMIS)

%% Generating the Sequence
[seq,states] = hmmgenerate(length(seq1'),estTR,estE);


%% ������ ����� ��Ȯ�� ����
likelystates = hmmviterbi(seq, estTR,estE);
score_1=sum(states==likelystates)/length(seq1)*100; % states-�ð迭�� ������ ����Ѵ�. 
%score_2=sum(states==likelystates2)/100; % ������ǽ������ΰ�� ������ �ſ� ����. 

writematrix([estTR,estE],'2_state_all_11.xlsx','WriteMode','append')
% writematrix(score_1,'Score_2states_10000','WriteMode','append')
%  writematrix([init_TRANS,estTR,estE],'Spring_all_10000_3states.xls','WriteMode','append')
%  writematrix(score_1,'Score_3states_10000','WriteMode','append')
% writematrix(seq1,'gensequences','WriteMode','append')
%writematrix([init_TRANS,estTR,estE],'Spring_all.xls',)

%xlswrite('Spring_all.xlsx',[estTR,estE])

%histc(states,1)  % no of state repeat

end 

toc
% n=10;
% klist=2:n;%the number of clusters you want to try
% x=seq1
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(x,myfunc,'CalinskiHarabasz','klist',klist)
% classes=kmeans(x,eva.OptimalK);