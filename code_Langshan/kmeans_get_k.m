clc;
clear;

load 'data_DHS.txt';
data = data_DHS(:,:);
[n,p]=size(data);

rng(1,'simdTwister')
K=8;
for k=2:K    
[lable,c,sumd,d]=kmeans(data,k,'dist','sqeuclidean');
sse1 = sum(sumd.^2);
D(k,1) = k;
D(k,2) = sse1;
end

figure
plot(D(2:end,1),D(2:end,2))
hold on;
plot(D(2:end,1),D(2:end,2),'or');
xlabel('Number of clusters (k)') 
ylabel('Sum of the squared errors') 
