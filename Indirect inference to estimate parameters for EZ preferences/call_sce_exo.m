function [ebeta,Rsq,ebeta1,Rsq1]=call_sce_exo(gamma,sigma,exo)
%% Recursive utility, AR(1) Income growth, AR(1) Stochastic volatility, Debt elastic roi
global M_ options_ oo_
rhos=exo(3);
sg=0.38; %annual
ss=exo(4);
mus=0;
rhog=exo(2);
Rstar=1.04; %Nominal interest rate
rstar=0.04;
mupi=0;
mug=0;
spi=exo(7);
rhopi=exo(6);
isigma=1/sigma;
beta=(1+mug)^(1/sigma) * (1+mupi)/Rstar;
bstar=-0.7442 ; % Uribe and Schmitt (2016) for canada
%bstar=6.6;
cstar=1+mug+(Rstar-(1+mug)*(1+mupi))*bstar;
vstar=cstar*((1-beta)/(1-beta*(1+mug)^(1-1/sigma)))^(1/(1-1/sigma));

save list sigma beta mug rhog bstar Rstar cstar vstar gamma rhos mus ss sg rstar mupi spi rhopi
dynare ez_uysv.mod


% Simulate data and run regression
% draw shock
burnin=500; %periods for convergence
track=1; %period to follow the household
HH=1500; %number of households
Msim=5;
Nexo=3;
shock_mat = zeros(burnin+track,Nexo); % which will go in to simulate data
period=burnin+track;
randn('seed',1235);
mu=[0 0 0];
sig=[sg^2 0 0;0 ss^2 0; 0 0 spi^2];
sig2=chol(sig);
volshock=zeros(period,HH,Nexo,Msim); % idiosyncratic shock for each hh across M simulations
for mm=1:Msim
    for pp=1:period
        volshock(pp,:,:,mm)=repmat(mu,HH,1)+randn(HH,Nexo)*sig2;
    end
end


es_pos=strmatch('es',M_.exo_names,'exact');
eg_pos=strmatch('eg',M_.exo_names,'exact');
epi_pos=strmatch('epi',M_.exo_names,'exact');
ecg_pos=strmatch('Ecg',M_.endo_names,'exact');
s_pos=strmatch('s',M_.endo_names,'exact');
g_pos=strmatch('g',M_.endo_names,'exact');
Eg_pos=strmatch('Eg',M_.endo_names,'exact');
c_pos=strmatch('C',M_.endo_names,'exact');
b_pos=strmatch('B',M_.endo_names,'exact');
pi_pos=strmatch('PI',M_.endo_names,'exact');
exppi_pos=strmatch('expPI',M_.endo_names,'exact');

estbeta=zeros(Msim,2);
estbeta1=zeros(Msim,3);
Rsq=zeros(Msim,1);
Rsq1=zeros(Msim,1);

for mm=1:Msim
    %mm
ecg_g=zeros(1+burnin+track,HH);
g=zeros(1+burnin+track,HH);
Eg=zeros(1+burnin+track,HH);
svol=zeros(1+burnin+track,HH);
c_g=zeros(1+burnin+track,HH);
b_g=zeros(1+burnin+track,HH);
pi_g=zeros(1+burnin+track,HH);
exppi=zeros(1+burnin+track,HH);

initial=oo_.dr.ys;
%initial(strmatch('B',M_.endo_names,'exact'))=0; %begin life with zero assets
%simulate data
for hh=1:HH
   % hh
    shock_mat(:,strmatch('es',M_.exo_names,'exact'))= volshock(:,hh,es_pos,mm);
    shock_mat(:,strmatch('eg',M_.exo_names,'exact'))= volshock(:,hh,es_pos,mm);
    shock_mat(:,strmatch('epi',M_.exo_names,'exact'))= volshock(:,hh,es_pos,mm);
    %sim_data=simult_(stochastic_steady_state',oo_.dr,shock,options_.order)'; %simulate time series for this HH
    sim_data=simult_(M_,options_,oo_.dr.ys,oo_.dr,shock_mat,options_.order)'; %simulate time series for this HH
    %sim_data=simult_(oo_.dr.ys,oo_.dr,shock_mat,options_.order)'; %simulate time series for this HH
    %sim_data=simult_(initial,oo_.dr,shock_mat,options_.order)'; %simulate time series for this HH
    ecg_g(:,hh)=sim_data(:,ecg_pos); %row one is the stochastic steady state
    g(:,hh)=sim_data(:,g_pos);
    Eg(:,hh)=sim_data(:,Eg_pos);
    svol(:,hh)=sim_data(:,s_pos);
    c_g(:,hh)=sim_data(:,c_pos);
    b_g(:,hh)=sim_data(:,b_pos);
    pi_g(:,hh)=sim_data(:,pi_pos);
    exppi(:,hh)=sim_data(:,exppi_pos);
end

%discard burnin period
Ecg_g=ecg_g(1+burnin+1:1+burnin+track,:);
G=g(1+burnin+1:1+burnin+track,:);
EG=Eg(1+burnin+1:1+burnin+track,:);
Galt=G(1,:)';
Svol=svol(1+burnin+1:1+burnin+track,:);
EPI=exppi(1+burnin+1:1+burnin+track,:);
ecg=(Ecg_g+1).*(1+G)-1; %normalize back

% %hpfilter the data
% [ecgt, ecgc]=hpfilter(ecg',1600); %ecg then g
% y=ecgc(:,1); %convert it into percent *100, but scaleing wont matter
% %[svolt,svolc]=hpfilter(sg*exp(Svol'),1600);
% [svolt,svolc]=hpfilter(exp(Svol'),1600); 
% x1=svolc;
% %scatter(x1,y,'filled');
% [epit,epic]=hpfilter(EPI',1600);
% x2=epic;
% [egt,egc]=hpfilter(EG',1600); 
% x3=egc;

y=ecg';
x1=exp(Svol');
x2=EPI';
x3=EG';

x=[ones(HH,1) x1 x2];
[betahat,BINT,residual,rint,STATS]=regress(y,x);
estbeta(mm,:)=betahat(2:3)';
Rsq(mm,1)= STATS(1);

x=[ones(HH,1) x1 x2 x3];
[betahat,BINT,residual,rint,STATS]=regress(y,x);
estbeta1(mm,:)=betahat(2:4)';
Rsq1(mm,1)= STATS(1);
end
ebeta=mean(estbeta,1)';
Rsq=mean(Rsq);
ebeta1=mean(estbeta1,1)';
Rsq1=mean(Rsq1);
end