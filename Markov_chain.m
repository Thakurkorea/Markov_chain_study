% % function fp = first_passage(P,E,x0,T)
clc
clear all

% transition_probabilities = [0.67	0.07	0.13	0.13
% 0.21	0.47	0.16	0.16
% 0.23	0.14	0.5	0.14
% 0.09	0.45	0.36	0.09]; 
% b='Autumn M44468';
transition_probabilities=[0.239	0.611	0.000	0.064	0.007	0.079
0.076	0.798	0.002	0.070	0.002	0.053
0.000	0.667	0.000	0.333	0.000	0.000
0.035	0.371	0.002	0.559	0.000	0.033
0.143	0.500	0.071	0.071	0.143	0.071
0.048	0.454	0.000	0.048	0.021	0.430];
% starting_value = 1; 
% chain_length = 1000;
  mc = dtmc(transition_probabilities,'StateNames',["U1" ,"U2", "U3", "U3" ,"U5",  "U6"]) % 1904P event
%  mc = dtmc(transition_probabilities,'StateNames',["P1" ,"P2", "P3" ,"P4", "P7", "P8"]) % 1904P event
% mc = dtmc(transition_probabilities,'StateNames',["H1"	"H2" "H3" "H4" "H5" "H6" "H7"	"H8" "H13" "H14"]) % 1904H
% mc = dtmc(transition_probabilities,'StateNames',["H1"	"H2" "H5" "H6" "H13" "H14"]) %1905Hevent
% mc = dtmc(transition_probabilities,'StateNames',["H1"	"H2" "H3" "H4" "H5" "H6" "H8" "H13" "H14" "H15" "H16"]) % 44467H
% mc = dtmc(transition_probabilities,'StateNames',["H1"	"H2" "H5" "H6" ]) % 44469H
% mc = dtmc(transition_probabilities,'StateNames',["H1"	"H2" "H3" "H4" "H5" "H6" "H7" "H8"  "H13" "H14" "H15" "H16"]) % 4indvH
% mc = dtmc(transition_probabilities,'StateNames',["U1" "U2"  "U5"  "U6"]) % 4indvH
figure;
g=digraph(mc.P);
% h=graphplot(mc,'ColorEdges',true);
h=graphplot(mc,'ColorEdges',true);
% h.EdgeLabel=round(g.Edges.Weight,2);
% h.EdgeColor=matrix
colormap('winter')
h.ArrowSize=6;
h.ArrowPosition=0.7;
h.EdgeFontSize=6;
% layout(h,'layered');
% layout(h);
% h.XData=h.XData
% h.YData=h.YData*2
h.EdgeFontWeight= 'bold'
% colorbar off
% axis equal
axis off
% axis tight
set(gcf,'color','w');
% figure;
%  graphplot(mc,'ColorNodes',true);
% % 
% mc1 = subchain(mc,1)
% % mc2 = subchain(mc,6);
[x1,t1] = asymptotics(mc)   %stationary state

for i=1:length(mc.StateNames);
% figure
ht(:,i)=hittime(mc,1,'Graph',false);
%  hp = hitprob(mc,i,'Graph',true)
axis tight
set(gcf,'color','w');
hold off
end
heat=heatmap(mc.StateNames,mc.StateNames,round(ht,2))
axp = struct(heat);
axp.Axes.XAxisLocation = 'top';
colormap(lines(9))
set(gcf,'color','w');
% 
% 
% %% markov sequence generate based on mc
% % X = simulate( mc , 450); 
% % figure;
% % simplot(mc,X);
% % figure;
% % simplot(mc,X,'Type','transition');
% %% plotting the transition plot without value
% % figure;
% % imagesc(mc.P);
% % colormap(jet);
% % axis square;
% % colorbar;
% 
% % %% find the simulated transition matrix
% % m = max(X);
% % n = numel(X);
% % y = zeros(m,1);
% % pp = zeros(m,m);
% % for k=1:n-1
% %     y(X(k)) = y(X(k)) + 1;
% %     pp(X(k),X(k+1)) = pp(X(k),X(k+1)) + 1;
% % end
% % pp = bsxfun(@rdivide,pp,y); pp(isnan(pp)) = 0;
% 
% %% heatmap
% figure
% xvalues=mc.StateNames;
% % xvalues={"V1" "V2" "V3" "V4"  "V6"  "V10"  "V11"   "V12"};
% yvalues=xvalues;
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'Colormap',copper,'GridVisible','off')
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'colormap',parula(10),'GridVisible','off')
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'colormap',jet,'GridVisible','off')
% h=heatmap(xvalues,yvalues,transition_probabilities)
% h.FontSize=12
% 
% h.ColorScaling = 'scaledrows';
% % colorbar off
% h.Title = b;
% set(gcf,'color','w');
% %%
% % figure;
% % distplot(mc,X,'Type','histogram','FrameRate',0.5);
% %     chain = zeros(1,chain_length);
% %     chain(1)=starting_value;
% %     for i=2:chain_length
% %         this_step_distribution = transition_probabilities(chain(i-1),:);
% %         cumulative_distribution = cumsum(this_step_distribution);
% %         r = rand();
% %         chain(i) = find(cumulative_distribution>r,1);
% %     end
% figure
% % hh=heatmap(xvalues,yvalues,ht, 'colormap',parula(10),'GridVisible','off')
% hh=heatmap(xvalues,yvalues,ht)
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'colormap',jet,'GridVisible','off')
% hh.ColorScaling = 'scaledrows';
% hh.Title = 'First hitting time for different habitats';
% hh.FontSize=12
% set(gcf,'color','w');
% 
% %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculation the centrality using adjecency matrix
% figure
% A=mc.P;
% b=A+A';
% G = graph(b,'omitselfloops')
% p = plot(G,'Layout','force','EdgeAlpha',0.5,'NodeColor','r');
% deg_ranks = (centrality(G,'degree','Importance',G.Edges.Weight)); % transpose for horrizontal row
% rank=deg_ranks'
% % deg_ranks = (centrality(G,'degree'))'
% % ucc = centrality(G,'closeness')
% 
% % Use discretize to place the nodes into 6 equally-spaced bins based on their centrality scores.
% edges = linspace(min(deg_ranks),max(deg_ranks),6);
% bins = discretize(deg_ranks,edges);
% p.MarkerSize = bins;
% p.NodeCData = deg_ranks;
% colormap jet
% colorbar
% % C1 = (centrality(digraph(mc.P), 'authorities'))'
% % C2 = (centrality(digraph(mc.P), 'hubs'))'
% C3 = (centrality(digraph(mc.P), 'pagerank'))'
% C4 = (centrality(digraph(mc.P), 'incloseness'))'
%  C5 = (centrality(digraph(mc.P), 'outcloseness'))'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


% figure
% 
% PP=[0.03	0.82	0.01	0.12	0	0.02
% 0.14	0.49	0	0.31	0	0.05
% 0.09	0.72	0	0.1	0	0.08
% 0.02	0.56	0	0.15	0	0.26
% 0.05	0.8	0	0.08	0	0.06
% 0.12	0.8	0	0	0.02	0.06
% 0.09	0.75	0	0.13	0	0.04];
% 
% xvalues=mc.StateNames;
% % xvalues={"U1" "U2" "U3" "U4"   "U7"   "U8"};
% yvalues={'SPF1903'	'SPM1905'	'SPM44467'	'SMRM44467'	'AUTM44467'	'AUTM44468'	'WINM44467'	};
% yvalues={'#1903F'	'#1905M'	'#44467 M'	'#44467M'	'#44467  M'	'#44468M'	'#44467  M'	};
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'Colormap',copper,'GridVisible','off')
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'colormap',parula(10),'GridVisible','off')
% % h=heatmap(xvalues,yvalues,transition_probabilities, 'colormap',jet,'GridVisible','off')
% h=heatmap(xvalues,yvalues,PP)
% h.FontSize=12
% axp = struct(h);
% axp.Axes.XAxisLocation = 'top';
% h.ColorScaling = 'scaledrows';
% colorbar off
% colormap('cool')
% % h.Title = b;
% set(gcf,'color','w');
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centrality with reaf field data without Adjacency matrix 
s=g.Edges.EndNodes(:,1);
t=g.Edges.EndNodes(:,2);
Wt=g.Edges.Weight;
GG=digraph(s,t,[])
% gp_ranks=(centrality(GG,'pagerank'))'  % pagerank transpose for row
% wbc = centrality(GG,'incloseness')'/sum(centrality(GG,'incloseness')')
% outclo_rank = centrality(GG,'outcloseness')'/sum(centrality(GG,'outcloseness')')
inclo_rank = centrality(GG,'incloseness','Cost',Wt)'/ sum(centrality(GG,'incloseness','Cost',Wt)')
outclo_rank = centrality(GG,'outcloseness','Cost',Wt)'/ sum(centrality(GG,'outcloseness','Cost',Wt)')
hub_rank = centrality(GG,'hubs','Importance',Wt)'
aut_rank = centrality(GG,'authorities','Importance',Wt)'
Pgrank = centrality(GG,'pagerank','Importance',Wt)'
between = centrality(GG,'betweenness','Cost',Wt)'/sum(centrality(GG,'betweenness','Cost',Wt)')
indegree = centrality(GG,'indegree','Importance',Wt)'/sum(centrality(GG,'indegree','Importance',Wt)')
% 
% %%%%%%%%%%%%%%correlation between ranks 
% ranks=round( [x1,
%     inclo_rank, 
%     outclo_rank,
%     hub_rank,
%     aut_rank,
%     Pgrank],2)
% rank_names={['Markov'], ['Incloseness'], ['Outcloseness'],'Hubs', 'Authorities', 'PageRank'};
% aaa=size(ranks);
% figure
% heatmap(mc.StateNames,rank_names,ranks)
% colormap('summer'); % Choose jet or any other color scheme
% colorbar % 
% set(gcf,'color','w')
% 
% 
%  imagesc(ranks)
% set(gca, 'XTick', 1:aaa(2)); % center x-axis ticks on bins
% set(gca, 'YTick', 1:aaa(1)); % center y-axis ticks on bins
% set(gca, 'XTickLabel',mc.StateNames ); % set x-axis labels
% set(gca, 'YTickLabel',rank_names); % set y-axis labels
% title('Node importance with different measures', 'FontSize', 10); % set title
% colormap('winter'); % Choose jet or any other color scheme
% colorbar % 
% set(gcf,'color','w')
% figure 
% cor= corr(ranks)
% imagesc(cor)
% set(gca, 'XTick', 1:aaa(2)); % center x-axis ticks on bins
% set(gca, 'YTick', 1:aaa(1)); % center y-axis ticks on bins
% set(gca, 'XTickLabel',mc.StateNames ); % set x-axis labels
% set(gca, 'YTickLabel',rank_names); % set y-axis labels
% title('Your Title', 'FontSize', 10); % set title
% colormap('jet'); % Choose jet or any other color scheme
% set(gcf,'color','w')
% 
% %%%%%%%%%%%%% rank plot %%%%%%%%%%%%%%%%
% 
% indeg_rank = centrality(GG,'betweenness','Cost',Wt)'/sum(centrality(GG,'betweenness','Cost',Wt)')
% % ranks=round( [x1,    inclo_rank,     outclo_rank,    hub_rank,    aut_rank,    Pgrank],2)
% indeg_rank=aut_rank
% %ranking the indivi
% rank=zeros(size(indeg_rank))
% for i =1:length(indeg_rank)
% a=find(indeg_rank==max(indeg_rank))
% indeg_rank(a)=0
% rank(a)=i
% end

%##################################################################
% % % generating the sequence and testing the sequences 
% a=[407, 439, 529 675 378 253 171 236 135 2050 3223];
% % 
% rng(1); % For reproducibility
% sequences = simulate(mc,253,'X0',[0 0 0 0 0 1000]);
% figure;
% simplot(mc,sequences);

% pp=[0.29	0.64	0	0.07	0	0
% 0.05	0.85	0.01	0.08	0	0.01
% 0	1	0	0	0	0
% 0	0.59	0	0.37	0	0.04
% 0	0	0	0	0	0
% 0	1	0	0	0	0];
% pp(:)





