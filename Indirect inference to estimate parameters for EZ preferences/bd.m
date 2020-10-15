function [distance,ebeta] = bd(parameterguess,truebeta,optwt,exo)
global gamma sigma ebeta
parameterguess
truebeta
gamma=parameterguess(1); % risk aversion guess
sigma=parameterguess(2); % ies guess
[ebeta,Rsq,ebeta1,Rsq1]=call_sce_exo(gamma,sigma,exo);
ebeta1
dis=truebeta'-ebeta1(1:2); % distance between the true and simulated beta
distance=dis'*optwt*dis %optimal weighting matrix
%distance=dis*dis'



