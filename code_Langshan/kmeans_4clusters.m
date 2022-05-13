clear;
clc;

load 'data_DHS.txt';
X = data_DHS(:,:);

rng(1,'simdTwister')
[idx,C,sumd] = kmeans(X,4)

figure
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',16)
hold on
plot(X(idx==2,1),X(idx==2,2),'m.','MarkerSize',16)
hold on
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',16)
hold on
plot(X(idx==4,1),X(idx==4,2),'b.','MarkerSize',16)

hold on
% % plot(C(:,1),C(:,2),'kx')
plot(C(:,1),C(:,2),'k^','MarkerSize',8)
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster Centroid')

h=gca;
xlim ([0, 15]);
ylim ([0, 10]);
set(h,'xtick',1:1:14);
set(h, 'ytick', 0:1:10);
% % grid on;
xlabel('Serial number, a.u.')
ylabel('Depth of half saturation (DHS) (mm)')
box on;