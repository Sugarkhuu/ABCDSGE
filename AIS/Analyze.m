load tuned_local.out;
parameters;
theta0 = lb_param_ub(:,2);
rep = rows(thetahatsLC);
for j = 1:6
    if j == 1
        contrib = thetahatsLC(1:rep,:);
    elseif j == 2    
        contrib = thetahatsLC50(1:rep,:);
    elseif j == 3    
        contrib = thetahatsLL(1:rep,:);
    elseif j == 4    
        contrib = thetahatsLL50(1:rep,:);
    elseif j == 5    
        contrib = thetahatsLQ(1:rep,:);
    elseif j == 6    
        contrib = thetahatsLQ50(1:rep,:);
    endif
    m = mean(contrib);
    s = std(contrib);
    e = contrib - repmat(theta0',rows(contrib),1);
    b = mean(e);
    e = e.^2;
    mse = mean(e);
    rmse = sqrt(mse);
    lb = lb_param_ub(:,1);
    ub = lb_param_ub(:,3);
    priormean = (ub+lb)'/2;
    priorsdev = sqrt(((ub-lb).^2)/12);
    priorsdev = priorsdev';
    priorbias = priormean - theta0';
    priorrmse = sqrt(priorbias.^2 + priorsdev.^2);
    percb = 100*b./theta0';
    percr = 100*rmse./theta0';
    mae = mean(abs(e));
    clabels = char("true", "mean", "pmean", "sdev.","psdev", "bias", "pbias","percb","rmse", "prmse","percr");
    rlabels = char(
    "alpha",
    "beta",
    "delta",
    "gam",
    "rho1",
    "sigma1",
    "rho2",
    "sigma2",
    "nss"
    );
    if j == 1
        printf("\n\nEstimation results (LC mean): rep %d\n", rep);
    elseif j == 2    
        printf("\n\nEstimation results (LC median): rep %d\n", rep);
    elseif j == 3    
        printf("\n\nEstimation results (LL mean): rep %d\n", rep);
    elseif j == 4    
        printf("\n\nEstimation results (LL median): rep %d\n", rep);
    elseif j == 5    
        printf("\n\nEstimation results (LQ mean): rep %d\n", rep);
    elseif j == 6    
        printf("\n\nEstimation results (LC median): rep %d\n", rep);
    endif
    prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; percb; rmse; priorrmse;percr]', rlabels, clabels);
    printf("\n\n");
endfor
printf("\n\n");
printf("90%% CI coverage: \n");
in_ci = (cilower <= theta0') & (ciupper >= theta0');
in_ci = mean(in_ci,1);
disp(in_ci);
printf("\n");

