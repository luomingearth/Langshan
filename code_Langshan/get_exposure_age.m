clear all;
clc;

file1 = '1K.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);

x1 = 0:0.1:16;% mm
F = 3.3/506.0;
t = 0.003;

X0 = [0 0];

modelfun = @(b,x1)((b(2)*exp(-1.0*b(1)*x1).*(exp(-1.0*t*((b(2)*exp(-1.0*b(1)*x1))+F)))+F)./((b(2)*exp(-1.0*b(1)*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,[0.0,0.0]);
ci = nlparci(beta,R,'covar',CovB); % 估计置信区间
sigma_mu=sqrt(CovB(1,1));
sigma_phi=sqrt(CovB(2,1));
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

figure()
h=gca;
xlim ([0, 16]);
ylim ([0, 1.3]);
set(h,'xtick',[0 2 4 6 8 10 12 14 16]);
set(h,'ytick',[0 0.2 0.4 0.6 0.8 1.0 1.2]);
%axis([0 16 0 1.2]);
box on;
p1 = plot(xdata1,ydata1,'c o') 
hold on
p2 = plot(x1,ypred1,'c-','LineWidth',1.5) 
plot(x1,[lower1;upper1],'c--','LineWidth',0.5) 
errorbar(xdata1,ydata1,err1,'c o') 

xlabel('Depth (mm)')
ylabel('Norm. Sens. Corr. IRSL (Ln/Tn)')
mu = beta(1);
phi = beta(2);
mu_1K_sigma = num2str(sigma_mu);
phi_1K_sigma = num2str(sigma_phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LSB1.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
p31 = plot(xdata1,ydata1,'k o','MarkerFaceColor','k') 
errorbar(xdata1,ydata1,err1,'k o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS0_10.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);

modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

p3 = plot(xdata1,ydata1,'r ^','MarkerFaceColor','red')
hold on
p4 = plot(x1,ypred1,'r ','LineWidth',1.5)
plot(x1,[lower1;upper1],'r--','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'r ^')
t_LS0_10 = beta;
t_LS0_10_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS30_40.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);

modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

p5 = plot(xdata1,ydata1,'r s') 
hold on
p6 = plot(x1,ypred1,'r -.','LineWidth',1.5)
plot(x1,[lower1;upper1],'r--','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'r s')
t_LS30_40 = beta;
t_LS30_40_sigma = num2str(sigma_t);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS70_80.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
%figure()
p7 = plot(xdata1,ydata1,'r o','MarkerFaceColor','r') 
hold on
p8 = plot(x1,ypred1,'r :.','LineWidth',1.5)
plot(x1,[lower1;upper1],'r--','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'r o')
t_LS70_80 = beta;
t_LS70_80_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS90_100.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p9 = plot(xdata1,ydata1,'r ^') 
hold on
p10 = plot(x1,ypred1,'r -','LineWidth',3)
plot(x1,[lower1;upper1],'r--','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'r ^')
t_LS90_100 = beta;
t_LS90_100_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS140_150.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB);
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p11 = plot(xdata1,ydata1,'r s','MarkerFaceColor','r') 
hold on
p12 = plot(x1,ypred1,'r -.','LineWidth',3)
plot(x1,[lower1;upper1],'r --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'r s')
t_LS140_150 = beta;
t_LS140_150_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS160_170.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB);
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p13 = plot(xdata1,ydata1,'g ^','MarkerFaceColor','g')
hold on
p14 = plot(x1,ypred1,'g -','LineWidth',1.5)
plot(x1,[lower1;upper1],'g --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'g ^')
t_LS160_170 = beta;
t_LS160_170_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS210_220.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
%figure()
p15 = plot(xdata1,ydata1,'g s')
hold on
p16 = plot(x1,ypred1,'g -.','LineWidth',1.5)
plot(x1,[lower1;upper1],'g --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'g s')
t_LS210_220 = beta;
t_LS210_220_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS290_300.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB);
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p17 = plot(xdata1,ydata1,'g o','MarkerFaceColor','g') 
hold on
p18 = plot(x1,ypred1,'g :','LineWidth',1.5)
plot(x1,[lower1;upper1],'g --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'g o')
t_LS290_300 = beta;
t_LS290_300_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS310_320.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p19 = plot(xdata1,ydata1,'b ^','MarkerFaceColor','b') 
hold on
p20 = plot(x1,ypred1,'b -','LineWidth',1.5)
plot(x1,[lower1;upper1],'b --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'b ^')
t_LS310_320 = beta;
t_LS310_320_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS370_380.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB);
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p21 = plot(xdata1,ydata1,'b s')
hold on
p22 = plot(x1,ypred1,'b -.','LineWidth',1.5)
plot(x1,[lower1;upper1],'b --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'b s')
t_LS370_380 = beta;
t_LS370_380_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS390_400.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);

modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

p23 = plot(xdata1,ydata1,'b o','MarkerFaceColor','b') 
hold on
p24 = plot(x1,ypred1,'b :','LineWidth',1.5)
plot(x1,[lower1;upper1],'b --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'b o')
t_LS390_400 = beta;
t_LS390_400_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS420_430.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

p25 = plot(xdata1,ydata1,'m ^','MarkerFaceColor','m') 
hold on
p26 = plot(x1,ypred1,'m -','LineWidth',1.5)
plot(x1,[lower1;upper1],'m --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'m ^')
t_LS420_430 = beta;
t_LS420_430_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS450_460.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
ci = nlparci(beta,R,'covar',CovB); 
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
p27 = plot(xdata1,ydata1,'m s') 
hold on
p28 = plot(x1,ypred1,'m -.','LineWidth',1.5)
plot(x1,[lower1;upper1],'m --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'m s')
t_LS450_460 = beta;
t_LS450_460_sigma = num2str(sigma_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS590_600.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x1)((phi*exp(-1.0*mu*x1).*(exp(-1.0*b*((phi*exp(-1.0*mu*x1))+F)))+F)./((phi*exp(-1.0*mu*x1))+F));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.0);
sigma_t=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x1,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;
%figure()
p29 = plot(xdata1,ydata1,'m o','MarkerFaceColor','m') 
hold on;
p30 = plot(x1,ypred1,'m :','LineWidth',1.5)
plot(x1,[lower1;upper1],'m --','LineWidth',0.5)
errorbar(xdata1,ydata1,err1,'m o')
t_LS590_600 = beta;
t_LS590_600_sigma = num2str(sigma_t);

h=gca;
xlim ([0, 17]);
ylim ([0, 1.3]);
set(h,'xtick',[0 2 4 6 8 10 12 14 16 18]);
set(h,'ytick',[0 0.2 0.4 0.6 0.8 1.0 1.2]);
legend1 = legend([p31 p1 p3 p5 p7 p9 p11 p13 p15 p17 p19 p21 p23 p25 p27 p29],{'-36\_-46B','1K','0\_10','30\_40','70\_80','90\_100','140\_150','160\_170','210\_220','290\_300','310\_320','370\_380','390\_400','420\_430','450\_460','590\_600'})
legend('boxoff')
hold on;
ah=axes('position',get(gca,'position'),'visible','off');
legend2 = legend(ah, [p2, p4, p6, p8, p10, p12, p14, p16, p18, p20, p22, p24, p26, p28, p30],{'','','','','','','','','','','','','','','',''})

set(legend1,'FontSize',9,'FontWeight','normal')
set(legend2,'FontSize',9,'FontWeight','normal')
legend('boxoff')
