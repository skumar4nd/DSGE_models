%% Get beta coefficients for the model across values of ies and rra
% Real consumption growth on real wage uncertainty and inflation

RRA=[-5,-2,-0.5,0.5,2,10,80];
IES=[0.1,0.5,0.8,1.2];
reg1=zeros(length(RRA),length(IES));
reg2=zeros(length(RRA),length(IES));
R=zeros(length(RRA),length(IES));

reg1_3x=zeros(length(RRA),length(IES));
reg2_3x=zeros(length(RRA),length(IES));
reg3_3x=zeros(length(RRA),length(IES));
R_3x=zeros(length(RRA),length(IES));

for kk=1:2
    if kk==1
        % sd as a measure of risk
        exo= xlsread('exog process.xlsx','5_5est','A3:G3');
    elseif kk==2
        % sd_beta as a measure of risk
        exo= xlsread('exog process.xlsx','5_5est','A8:G8');
    end
    for ii=1:length(RRA)
        gamma=RRA(ii)
        for ii=1:length(IES)
            sigma=IES(ii)
            [ebeta,Rsq,ebeta1,Rsq1]=call_sce_exo(gamma,sigma,exo);
            reg1(ii,ii)=ebeta(1); %unc
            reg2(ii,ii)=ebeta(2); % inf
            R(ii,ii)=Rsq;
            reg1_3x(ii,ii)=ebeta1(1); %unc
            reg2_3x(ii,ii)=ebeta1(2); % inf
            reg3_3x(ii,ii)=ebeta1(3); % mean wage growth
            R_3x(ii,ii)=Rsq1;
            filename=['test' num2str(gamma) '_' num2str(sigma) '_' num2str(kk) '.mat'];
            save(filename,'reg1','reg2','R','reg1_3x','reg2_3x','reg3_3x','R_3x','gamma','sigma');
        end
    end
filename1=['full' num2str(kk) '.mat'];    
save(filename1); 

reg1=zeros(length(RRA),length(IES));
reg2=zeros(length(RRA),length(IES));
R=zeros(length(RRA),length(IES));
reg1_3x=zeros(length(RRA),length(IES));
reg2_3x=zeros(length(RRA),length(IES));
reg3_3x=zeros(length(RRA),length(IES));
R_3x=zeros(length(RRA),length(IES));   
end
save('truevalues.mat')

%% Plot regression coefficients
% beta1, coefficient on wage uncertainty across values of risk aversion
hold on
for ii=1:length(IES)
    scatter(RRA(1:end-1),reg1_3x(1:end-1,ii),'filled')
end
xticks([-5 -2 -0.5 0.5 2 10])
xlabel('Risk aversion')
ylabel('Regression coefficient on wage uncertainty')
legend('IES=0.1','IES=0.5','IES=0.8','IES=1.2')
hold off

hold on
for ii=1:length(IES)
    plot(RRA(1:end-1),reg1_3x(1:end-1,ii),'Linewidth',2)
end
xticks([-5 -2 -0.5 0.5 2 10])
xlabel('Risk aversion')
ylabel('Regression coefficient on wage uncertainty')
legend('IES=0.1','IES=0.5','IES=0.8','IES=1.2')
hold off

% beta2, coefficient on inflation expecttaions across values of ies
hold on
for ii=1:length(RRA)-1
    scatter(IES,reg2_3x(ii,:),'filled')
end
xticks([0.1 0.5 0.8 1.2])
xlabel('IES')
ylabel('Regression coefficient on expected inflation')
legend('RRA=-5','RRA=-2','RRA=-0.5','RRA=0.5','RRA=2','RRA=10')
hold off

hold on
for ii=1:length(RRA)-1
    plot(IES,reg2_3x(ii,:),'Linewidth',2)
end
xticks([0.1 0.5 0.8 1.2])
xlabel('IES')
ylabel('Regression coefficient on expected inflation')
legend('RRA=-5','RRA=-2','RRA=-0.5','RRA=0.5','RRA=2','RRA=10')
hold off
