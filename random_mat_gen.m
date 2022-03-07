
clc
clear all
n=10;
for k=1:3
no_of_states=2;
        for i=no_of_states
            for j=no_of_states
                x1 = randi([1, 1000], [i,j]);
                if i==2
                y1=sum(x1(1,:));
                y2=sum(x1(2,:));
                z=[[x1(1,:)/y1];[x1(2,:)/y2]]; % 2states
                elseif i==3
                y1=sum(x1(1,:));
                y2=sum(x1(2,:));
                y3=sum(x1(3,:));
                z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3]];  %3states
                elseif i==4
                y1=sum(x1(1,:));
                y2=sum(x1(2,:));
                y3=sum(x1(3,:));
                y4=sum(x1(4,:));
                z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4]];  %4states 
                elseif i==5
                y1=sum(x1(1,:));
                y2=sum(x1(2,:));
                y3=sum(x1(3,:));
                y4=sum(x1(4,:));
                y5=sum(x1(5,:));      
                z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4];[x1(5,:)/y5]];% state 5
                else
                    error("Check the state for initial TPM")
                end
            end
        end
             
       
        
            
        init_TRANS=round(z,2)
end

% no_of_states=3;
%         for i=no_of_states
%             for j= n
%                 x1 = randi([1, 1000], [i,j]);
%                 if i==2
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                     z=[[x1(1,:)/y1];[x1(2,:)/y2]]; % 2states
%                 elseif i==3
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                     z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3]];  %3states
%                 elseif  i==4
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                      y4=sum(x1(4,:));
%                     z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4]];  %4states 
%                 elseif i==5
%                     y1=sum(x1(1,:));
%                     y2=sum(x1(2,:));
%                      y3=sum(x1(3,:));
%                      y4=sum(x1(4,:));
%                      y5=sum(x1(5,:));
%                      z=[[x1(1,:)/y1];[x1(2,:)/y2];[x1(3,:)/y3];[x1(4,:)/y4];[x1(5,:)/y5]];
%                 else
%                     error('Check the state number')
%                 end
%             end
%         end
%                       
% init_TRANS=round(z,2);
%     
% init_TRANS
% % 0.0200    0.1600    0.1600    0.1300    0.1500    0.1100    0.1500    0.1200
% %     0.0700    0.2300    0.2300    0.0300    0.1900    0.0100    0.1600    0.0900
% %     0.1100    0.0300    0.1000    0.0900    0.2000    0.1800    0.1600    0.1400
