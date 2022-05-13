clear all;
clc;

x = 0:0.1:30;% mm
F = 3.3/506;
u = 1.1044;
phi = 2630.9;

%%
file1 = 'LS0_10.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);
sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS0_10 = beta;
erosion_LS0_10_sigma = num2str(sigma_e);
%%
file1 = 'LS30_40.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS30_40 = beta;
erosion_LS30_40_sigma = num2str(sigma_e);
%%
file1 = 'LS70_80.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS70_80 = beta;
erosion_LS70_80_sigma = num2str(sigma_e);
%%
file1 = 'LS90_100.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS90_100 = beta;
erosion_LS90_100_sigma = num2str(sigma_e);
%%
file1 = 'LS140_150.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS140_150 = beta;
erosion_LS140_150_sigma = num2str(sigma_e);
%%
file1 = 'LS160_170.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS160_170 = beta;
erosion_LS160_170_sigma = num2str(sigma_e);
%%
file1 = 'LS210_220.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS210_220 = beta;
erosion_LS210_220_sigma = num2str(sigma_e);
%%
file1 = 'LS290_300.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS290_300 = beta;
erosion_LS290_300_sigma = num2str(sigma_e);
%%
file1 = 'LS310_320.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);
sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS310_320 = beta;
erosion_LS310_320_sigma = num2str(sigma_e);
%%
file1 = 'LS370_380.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS370_380 = beta;
erosion_LS370_380_sigma = num2str(sigma_e);
%%
file1 = 'LS390_400.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS390_400 = beta;
erosion_LS390_400_sigma = num2str(sigma_e);
%%
file1 = 'LS420_430.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS420_430 = beta;
erosion_LS420_430_sigma = num2str(sigma_e);
%%
file1 = 'LS450_460.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS450_460 = beta;
erosion_LS450_460_sigma = num2str(sigma_e);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'LS590_600.txt';
a= load(file1);
xdata1 = a(:,1);ydata1=a(:,2);
err1 = a(:,3);
modelfun = @(b,x) hypergeom(1,1+(F./(u.*b)),-1.0*(phi*exp(-1.0*u*x)./(u.*b)));
[beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(xdata1,ydata1,modelfun,0.001);

sigma_e=sqrt(CovB);
[ypred1,delta1] = nlpredci(modelfun,x,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower1 = ypred1 - delta1;
upper1 = ypred1 + delta1;

erosion_LS590_600 = beta;
erosion_LS590_600_sigma = num2str(sigma_e);
