% trans = [0.36,0.64;
%       0.63,0.37];
% emis = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;
%    1/10, 1/10, 1/10, 1/10, 1/10, 1/2];
trans = [0.272221	0.727791
0.675706	0.32435]
emis =[0.700203	0.065392	0.128901	0.010434	0.080405	0.014605
0.714036	0.064215	0.132276	0.013901	0.065886	0.009182]


% seq1 = hmmgenerate(100,trans,emis);
% % seq2 = hmmgenerate(200,trans,emis);
% % seqs = {seq1,seq2};
% [estTR,estE] = hmmtrain(seq1,trans,emis)
%%
input_data=xlsread('HMM_individuals.xlsx');
Number_of_combination=18;
itr= 1; % iteration for each cases
N_max=Number_of_combination*itr;
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
n=8; % number of Events

    no_of_states=2;
        for i=no_of_states
            for j=no_of_states
                x1 = randi([1, 1000], [i,j]);
                y1=sum(x1(1,:));
                y2=sum(x1(2,:));
                %y3=sum(x1(3,:));
                %y4=sum(x1(4,:));
                %y5=sum(x1(5,:));
                
                z=[[x1(1,:)/y1];[x1(2,:)/y2]]; % 2states
                %z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3]];  %3states
                %z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4]];  %4states        
                %z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4];[x1(5,:)/y5]];%
                init_TRANS=round(z,2);
            end
        end
%%%%%%%%%%%%%%%
%% Initial Emission matrix
%%%%%%%%%%%%%%%
init_EMIS = [ones(no_of_states,n)/n];

%% Training the sequence
[estTR,estE] = hmmtrain(seq1',init_TRANS,init_EMIS)

%% Generating the Sequence
[seq,states] = hmmgenerate(length(seq1'),estTR,estE);


%% 추정된 행렬의 정확도 점수
likelystates = hmmviterbi(seq, estTR,estE);
score_1=sum(states==likelystates)/length(seq1)*100; % states-시계열로 점수를 계산한다. 
%score_2=sum(states==likelystates2)/100; % 랜덤사건시퀀스인경우 점수가 매우 낮다. 

% writematrix([init_TRANS,estTR,estE],'Spring_all_1000_2states.xls','WriteMode','append')
% writematrix(score_1,'Score_2states_100','WriteMode','append')
% %  writematrix([init_TRANS,estTR,estE],'Spring_all_1000_3states.xls','WriteMode','append')
% %  writematrix(score_1,'Score_3states_1000','WriteMode','append')
% writematrix(seq1,'gensequences','WriteMode','append')
%writematrix([init_TRANS,estTR,estE],'Spring_all.xls',)

%xlswrite('Spring_all.xlsx',[estTR,estE])


end 
