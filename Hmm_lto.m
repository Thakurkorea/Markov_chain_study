%% 사건-시계열에 대한 히든마코프 모델 구성하기
% 전이확률(TRANS)과 출력확률(EMIS)을 알고 있는 경우에 
% 관측 사건들의 시계열을 생성하는 방법:
% 2개의 히든상태가 있고 6개의 관측가능 사건이 있는 경우:
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


%% 추정된 행렬의 정확도 점수
likelystates = hmmviterbi(seq, estTR,estE);
score_1=sum(states==likelystates)/length(seq1)*100; % states-시계열로 점수를 계산한다. 
%score_2=sum(states==likelystates2)/100; % 랜덤사건시퀀스인경우 점수가 매우 낮다. 

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