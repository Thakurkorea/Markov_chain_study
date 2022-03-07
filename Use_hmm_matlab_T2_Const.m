%% 사건-시계열에 대한 히든마코프 모델 구성하기
% 전이확률(TRANS)과 출력확률(EMIS)을 알고 있는 경우에 
% 관측 사건들의 시계열을 생성하는 방법:
% 2개의 히든상태가 있고 6개의 관측가능 사건이 있는 경우:
tic
clc;
clear 
close all;
%% Read the sequence file

input_data=xlsread('HMM_individuals.xlsx');

%% Get the individual sequence 
Number_of_combination=18;
itr= 2; % iteration for each cases
N_max=Number_of_combination*itr;
n=8; % number of Events
no_of_states=2;
for k=1:N_max;
    if k<N_max/Number_of_combination+1
        sequece=input_data(:,1);  %Spring all
    elseif N_max/Number_of_combination <k && 2*N_max/Number_of_combination+1>k
        sequece=input_data(:,3);  % Sp1905
    elseif 2*N_max/Number_of_combination <k && 3*N_max/Number_of_combination+1>k
        sequece=input_data(:,5);  % Sp1904
    elseif 3*N_max/Number_of_combination <k && 4*N_max/Number_of_combination+1>k
    sequece=input_data(:,7);  % Sp1903
    elseif 4*N_max/Number_of_combination <k && 5*N_max/Number_of_combination+1>k
    sequece=input_data(:,9);  % Sp44469
    elseif 5*N_max/Number_of_combination <k && 6*N_max/Number_of_combination+1>k
    sequece=input_data(:,11);  % Sp44467
    elseif 6*N_max/Number_of_combination <k && 7*N_max/Number_of_combination+1>k
    sequece=input_data(:,13);  % win44467
    elseif 7*N_max/Number_of_combination <k && 8*N_max/Number_of_combination+1>k
    sequece=input_data(:,15);  % Sum_all
    elseif 8*N_max/Number_of_combination <k && 9*N_max/Number_of_combination+1>k
    sequece=input_data(:,17);  % Sum1906
    elseif 9*N_max/Number_of_combination <k && 10*N_max/Number_of_combination+1>k
    sequece=input_data(:,19);  % Sum1903
    elseif 10*N_max/Number_of_combination <k && 11*N_max/Number_of_combination+1>k
    sequece=input_data(:,21);  % sum44467
    elseif 11*N_max/Number_of_combination <k && 12*N_max/Number_of_combination+1>k
    sequece=input_data(:,23);  % Aut_all
    elseif 12*N_max/Number_of_combination <k && 13*N_max/Number_of_combination+1>k
    sequece=input_data(:,25);  % Aut_44467
    elseif 13*N_max/Number_of_combination <k && 14*N_max/Number_of_combination+1>k
    sequece=input_data(:,27);  % Aut_44468
    elseif 14*N_max/Number_of_combination <k && 15*N_max/Number_of_combination+1>k
    sequece=input_data(:,29);  % Aut_1906
    elseif 15*N_max/Number_of_combination <k && 16*N_max/Number_of_combination+1>k
    sequece=input_data(:,31);  % Male_all_Spring
    elseif 16*N_max/Number_of_combination <k && 17*N_max/Number_of_combination+1>k
    sequece=input_data(:,33);  % Male_summer
    else
    sequece=input_data(:,35);  % Mall_all_Autumn
    end
seq1=rmmissing(sequece);
% multiple component
%% TRANS와 EMIS  

%        for i=no_of_states
%             for j=no_of_states
%                 x1 = randi([1, 1000], [i,j]);
%                 if i==2
%                 y1=sum(x1(1,:));
%                 y2=sum(x1(2,:));
%                 z=[[x1(1,:)/y1];[x1(2,:)/y2]]; % 2states
%                 elseif i==3
%                 y1=sum(x1(1,:));
%                 y2=sum(x1(2,:));
%                 y3=sum(x1(3,:));
%                 z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3]];  %3states
%                 elseif i==4
%                 y1=sum(x1(1,:));
%                 y2=sum(x1(2,:));
%                 y3=sum(x1(3,:));
%                 y4=sum(x1(4,:));
%                 z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4]];  %4states 
%                 elseif i==5
%                 y1=sum(x1(1,:));
%                 y2=sum(x1(2,:));
%                 y3=sum(x1(3,:));
%                 y4=sum(x1(4,:));
%                 y5=sum(x1(5,:));      
%                 z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4];[x1(5,:)/y5]];% state 5
%                 else
%                     error("Check the state for initial TPM")
%                 end
%             end
%         end
%              
%             
% init_TRANS=z;
init_TRANS=[0.005073481	 0.994926519
0.204262072	0.795737928];
%%%%%%%%%%%%%%%
%% Initial Emission matrix
%%%%%%%%%%%%%%%
% % init_EMIS = [ones(no_of_states,n)/n];
%         for i=no_of_states
%             for j= n
%                 x1 = randi([1, 1000], [i,j]);
%                 if i==2
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                     zz=[[x1(1,:)/y1];[x1(2,:)/y2]]; % 2states
%                 elseif i==3
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                     zz=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3]];  %3states
%                 elseif  i==4
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                      y4=sum(x1(4,:));
%                     zz=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4]];  %4states 
%                 elseif i==5
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                      y4=sum(x1(4,:));
%                      y5=sum(x1(5,:));
%                      zz=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4];[x1(5,:)/y5]];
%                 else
%                     error('Check the state number')
%                 end
%             end
%         end
%                       
% init_EMIS=zz;
%% Fixed init_EMIS
init_EMIS=[0.710615341	0.067822004	0.123909435	0.011049632	0	0	0.072048375	0.014555212
0.708434134	0.067690476	0.12377901	0.011042002	0	0	0.074528758	0.01452562];


%% Training the sequence
[estTR,estE] = hmmtrain(seq1',init_TRANS,init_EMIS);

%% Generating the Sequence
[seq,states] = hmmgenerate(length(seq1'),estTR,estE);


%% 추정된 행렬의 정확도 점수
likelystates = hmmviterbi(seq, estTR,estE);
score_1=sum(states==likelystates)/length(seq1)*100; % states-시계열로 점수를 계산한다. 
%score_2=sum(states==likelystates2)/100; % 랜덤사건시퀀스인경우 점수가 매우 낮다. 

writematrix([estTR,estE],'2s_Splowestval_spring_all.xlsx','WriteMode','append')
% writematrix([estTR,estE],'2_state_all_10000.xlsx','WriteMode','append')
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