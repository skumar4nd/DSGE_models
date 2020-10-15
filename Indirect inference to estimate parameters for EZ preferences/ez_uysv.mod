% Real consumption growth on real wage uncertainty and inflation
var V C B Ecg t2 g Eg s PI expPI;
%var V C B Ecg t2 g PI expPI;

varexo eg es epi;
%varexo eg epi;

parameters sigma beta mug rhog bstar Rstar cstar vstar gamma rhos mus ss sg mupi spi rhopi;

load list;
set_param_value('sigma',sigma);
set_param_value('beta',beta);
set_param_value('gamma',gamma);
set_param_value('mug',mug);
set_param_value('rhog',rhog);
set_param_value('bstar',bstar);
set_param_value('Rstar',Rstar);
set_param_value('cstar',cstar);
set_param_value('vstar',vstar);
set_param_value('mus',mus);
set_param_value('ss',ss);
set_param_value('sg',sg);
set_param_value('rhos',rhos);
set_param_value('mupi',mupi);
set_param_value('spi',spi);
set_param_value('rhopi',rhopi);

model;

% (1) Value function
V^(1-(1/sigma))= (1-beta)*C^(1-(1/sigma)) + beta* (1+g)^(1-(1/sigma)) * t2^((1-(1/sigma))/(1-gamma));

% (2) Euler equation capital
C^(-1/sigma) = beta* (1+g)^(-1/sigma) * C(+1)^(-1/sigma) * Rstar * (1+PI(+1))^-1 * (V(+1)/t2^(1/(1-gamma)))^((1-sigma*gamma)/sigma);
%Rstar is gross nominal interest rate

% (3) Budget constraint
C+B*(1+g)= 1+g+ Rstar * B(-1);

% (4) Expected real consumption growth
Ecg=(C(+1)/C)-1;

% (5) EZ term
t2=V(+1)^(1-gamma);

% (6) Expected inflation
expPI=PI(+1);

% (7) Inflation
PI=mupi *(1-rhopi)+rhopi*PI(-1) + epi;

% (8) Growth rate of income
g=mug *(1-rhog)+rhog*g(-1) + exp(s(-1))*eg;
%g=mug *(1-rhog)+rhog*g(-1) + eg;

% (9) Volatility of growth
s=mus *(1-rhos)+rhos*s(-1) + es;

% (10) Expected mean growth of income
Eg=g(+1);
end;

initval;
V=vstar;
C=cstar;
B=bstar;
Ecg=0;
g=mug;
s=mus;
PI=mupi;
expPI=mupi;
t2=vstar^(1-gamma);
Eg=mug;
end;

shocks;
var eg; stderr sg; 
var es; stderr ss;
var epi; stderr spi; 
end;

options_.noprint=1;
steady;
check;

stoch_simul(order=3,irf=20,noprint,nograph,nomoments,nofunctions,nocorr,pruning);
%stoch_simul(order=3,irf=20,noprint,pruning);
%stoch_simul(order=3,pruning,irf=0,periods=5200,simul_replic=1000);
%simul_replic: Number of series to simulate when empirical moments are requested