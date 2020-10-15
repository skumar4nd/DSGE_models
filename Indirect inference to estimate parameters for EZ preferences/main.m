%% Estimate Income and inlfation process for Structural model using SMM
close all; clear all;
global key
% Moments come from Survey fo Consumer Expectations
% A total of 10 parameters are to be jointly estimated (3 equations parameters and 2 structural parameters)
% Inflation and income process parameters + RRA and IES
% Targets for risk aversion and ies from SCE regressions
format short % 4 digits after decimal
sd= xlsread('sce_regcoef.xlsx','sd_real'); % column: beta, sd, varcovar
sd_beta= xlsread('sce_regcoef.xlsx','sdb_real'); % column: beta, sd, varcovar
%iqr_beta= xlsread('sce_regcoef.xlsx','iqr_beta'); % column: beta, sd, varcovar
%save('sce_regcoef.mat');
%load('sce_regcoef.mat');
truebeta_sd=[sd(1,1), sd(2,1)];  % regression coefficient from the data to be matched
weights_sd=inv([sd(1,2)^2,sd(1,3); sd(1,3),sd(2,2)^2]); %inverse of var-covar matrix=optimal weighting matrix
truebeta_sdbeta=[sd_beta(1,1), sd_beta(2,1)];  % regression coefficient from the data to be matched
weights_sdbeta=inv([sd_beta(1,2)^2,sd_beta(1,3); sd_beta(1,3),sd_beta(2,2)^2]); %inverse of var-covar matrix
% truebeta_iqrbeta=[iqr_beta(1,1), iqr_beta(2,1)];  % regression coefficient from the data to be matched
% weights_iqrbeta=inv([iqr_beta(1,2)^2,iqr_beta(1,3); iqr_beta(1,3),iqr_beta(2,2)^2]); %inverse of var-covar matrix
clearvars sd sd_beta
key=2; % exogenous process parametrization: parametric- beta distribution
if key==1
    % sd as a measure of risk
    exo= xlsread('exog process.xlsx','5est','A3:G3');
elseif key==2
    % sd_beta as a measure of risk
    exo= xlsread('exog process.xlsx','5est','A8:G8');
end

rra=13;
ies=0.9;
parameterguess(1)= rra;
parameterguess(2)= ies;

options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
A=[];B=[];Aeq=[];Beq=[];NONLCON=[];
LB=[0.00001, 0.00001];
UB=[Inf, Inf];

% % Indirect inference
% % 1. sd
% targets(1)= truebeta_sd(1);
% targets(2)= truebeta_sd(2);
% % key=1; % exogenous process parametrization: non parametric
% %[est_sd,distance_sd] = fmincon(@bd,parameterguess,A,B,Aeq,Beq,LB,UB,NONLCON,options,targets,weights_sd);
% % save('sce_sd_bnd.mat');
% [est_sd,dis_sd] = fminsearch(@bd,parameterguess,options,targets,weights_sd,exo);
% save('sce_sd_unbnd_allreal.mat');


% 2. sd_beta
targets(1)= truebeta_sdbeta(1);
targets(2)= truebeta_sdbeta(2);
% [est_sdbeta,distance_sdbeta] = fmincon(@bd,parameterguess,A,B,Aeq,Beq,LB,UB,NONLCON,options,targets,weights_sdbeta,exo);
% %save('sce_sdbeta_bnd.mat');
% save('sce_sdbeta_bnd_allreal.mat');
[est_sdbeta,distance_sdbeta] = fminsearch(@bd,parameterguess,options,targets,weights_sdbeta,exo);
save('sce_sdbeta_unbnd_allreal.mat');

% % 3. iqr_beta
% targets(1)= truebeta_iqrbeta(1);
% targets(2)= truebeta_iqrbeta(2);
% [parametersestimation_iqrbeta,distance_iqrbeta] = fminsearch(@paramdistance,parameterguess,options,targets,weights_iqrbeta);
% save('sce_iqrbeta.mat');
