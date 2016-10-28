parameters;
theta0 = lb_param_ub(:,2);
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
priormean = (ub+lb)'/2;
priorsdev = sqrt(((ub-lb).^2)/12);
priorsdev = priorsdev';
priorbias = priormean - theta0';
priorrmse = sqrt(priorbias.^2 + priorsdev.^2);


relbs = zeros(9,1);
relrmses = relbs;
rmses = relbs;
incis = relbs;
i = 1;
load ./NN1PDM1_1000/tuned_local.out;
contrib = thetahatsLL;
m = mean(contrib);
s = std(contrib);
e = contrib - repmat(theta0',rows(contrib),1);
b = mean(e);
e = e.^2;
mse = mean(e);
rmse = sqrt(mse);
relb = b ./ theta0';
relrmse = rmse ./ theta0';
relbs(:,i) = relb;
relrmses(:,i) = relrmse;
rmses(:,i) = rmse;
in_ci = (cilower <= theta0') & (ciupper >= theta0');
incis(:,i) = mean(100*in_ci)';
format('bank');
100*relbs
mean(abs(ans))
100*relrmses
mean(ans)
%rmses
incis
mean(incis)

