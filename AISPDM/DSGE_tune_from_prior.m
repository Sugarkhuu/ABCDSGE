outfile = "tunePRIOR.out";
mc_reps = 100; % number of MC reps
nworkers = 30;  % number of worker MPI ranks

% controls for creating the adaptive importance sampling density
iters = 10;
initialparticles = nworkers*round(300/nworkers); % number to take from sample from prior
nparticles = nworkers*round(300/nworkers); % number per round
particlequantile = 20; % keep the top % of particles
verbose = false;

% controls for drawing the final sample from mixture of AIS and prior
mixture = 0.5; % proportion sampled from original prior 
AISdraws = nworkers*round(5000/nworkers); # number of draws from final AIS density

% controls for the nonparametric fits
%nneighbors = 300;    
nbw = 30;

% design
parameters; % loaded from Common to ensure sync with Gendata

lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
prior_params = [lb ub];
theta0 = lb_param_ub(:,2); % original form
nparams = rows(theta0);

% which statistics to use
load selected; % selected statistics
asbil_selected = selected;

setupmpi; % sets comm world, nodes, node, etc.
asbil_theta = theta0; setupdynare; % sets structures and RNG for simulations
MPI_Barrier(CW);
warning ( "off") ;


% number of particles for each node
particles_per_node = floor(nparticles/(nodes-1));

%  frontend: load data and make containers
if !node
    % here, you need to provide code that defines USERthetaZ,
    % which is a reasonably large number set of
    % [theta  Z] where theta is a draw from prior
    % and Z is the output of aux_stat
    load simdata.paramspace;
    USERthetaZ = clean_data(simdata);
    % containers
    errors = zeros(mc_reps, nparams,nbw);
    in_ci = zeros(9,nbw);
    rmses = zeros(9,nbw);
    cicoverage = rmses;
endif



for rep = 1:mc_reps
    % the 'true' Zn
    if node==1 % simulate on node1 (can't do it on 0, no *.modfile
        % next line for tuning to true theta0
        %asbil_theta = theta0;
        % next two lines for tuning to draws from prior
        ok = false;
        while !ok    
            asbil_theta = sample_from_prior();
            theta0 = asbil_theta;
            USERsimulation;
            Zn = aux_stat(data);
            ok = Zn(1,:) != -1000;
        endwhile	
        Zn = Zn(asbil_selected,:);
        for i = 2:nodes-1
            MPI_Send(Zn, i, mytag, CW);
        endfor	
        MPI_Send(Zn, 0, mytag, CW);
        MPI_Send(theta0, 0, mytag+1, CW);
        MPI_Send(data, 0, mytag+2, CW);
    else % receive it on the other nodes
        Zn = MPI_Recv(1, mytag, CW);
        if  !node
            theta0 = MPI_Recv(1, mytag+1, CW);
            realdata = MPI_Recv(1, mytag+2, CW);
        endif    
    endif
    MPI_Barrier(CW);    
    Zn = Zn';
    
    % call the algoritm that gets AIS particles
    if ! node
            tic;
            printf("Starting CreateAIS\n");
    endif
    CreateAIS; # this gets the particles

    % now draw from the AIS density
    reps_per_node = round(AISdraws/(nodes-1));
    if ! node
            toc;
            printf("Starting SampleFromAIS\n");
            tic;
    endif

    SampleFromAIS; # this samples from AIS density


    % see the results
    if !node
        toc;
        printf("Starting fit and CI\n");
        thetas = contribs(:,1:nparams);
        Zs = contribs(:, nparams+1:end);
        test = sum(Zs,2) != 0;
        thetas = thetas(test,:);
        Zs = Zs(test,:);
        Z = [Zn; Zs];
function pdm = makepdm(thetas, realdata)        
        % add the parameter-dependent moments
        y = realdata(:,1);
        c = realdata(:,2);
        n = realdata(:,3);
        r = realdata(:,4);
        w = realdata(:,5);
        pdm = zeros(rows(thetas),7);
        for i = 1:rows(thetas)
            % break out params    
            alpha = thetas(i,1);
            beta = thetas(i,2);
            delta = thetas(i,3);
            gam = thetas(i,4);
            rho_z = thetas(i,5);
            sig_z = thetas(i,6);
            rho_eta = thetas(i,7);
            sig_eta = thetas(i,8);
            nss = thetas(i,9);
            % recover psi
            c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
            kss = nss/c1;
            iss = delta*kss;
            yss = kss^alpha * nss^(1-alpha);
            css = yss - iss;
            psi =  (css^(-gam)) * (1-alpha) * (kss^alpha) * (nss^(-alpha));
            % get white noise shock to prefs
            e = log(w)-rho_eta*lag(log(w),1) - gam*(log(c) - rho_eta*lag(log(c),1))-log(psi)*(1-rho_eta);
            e = e/sig_eta;
            e = e(3:end,:);
            pdm(i,1) = mean(e);
            pdm(i,2) = mean(e.^2 - 1);
            pdm(i,3) = mean(e(1:end-1,:).*e(2:end,:));
            e1 = e;
            % now the Euler eqn
            %cc = c ./lag(c,1);
            %e = (beta*cc.^(-gam).*(1 + (1-alpha)*lag(r,1).*y./(lag(n,1).*lag(w,1))-delta))-1;
            %e = e(2:end,:);
            %pdm(i,5) = mean(e);
            % production function
            k = alpha*lag(n,1).*lag(w,1)./lag(r,1)/(1-alpha);
            lnz = log(y) - alpha*log(k) - (1-alpha)*log(n);
            e = lnz - rho_z*lag(lnz,1);
            e = e(3:end,:);
            e = e/sig_z;
            pdm(i,4) = mean(e);
            pdm(i,5) = mean(e.^2 - 1);
            pdm(i,6) = mean(e(1:end-1,:).*e(2:end,:));
            pdm(i,7) = mean(e.*e1); % no cross correlation between errors
            % law of motion k: good for delta
            %invest = y - c;
            %e = lag(invest,1) + (1 - delta)*lag(k,1) - k;
            %e = e(2:end,:);
            %pdm(i,5) = mean(e);
        endfor
        pdm = [zeros(1,size(pdm,2)); pdm]; % stack real data stats (zeros) on top
endfunction


        pdm = makepdm(thetas, realdata);
        Z = [Z pdm]; 

        % first pre-whiten using all draws
        %q = quantile(Z,0.99);
        %test = Z < q;
        %Z = test.*Z + (1-test).*q;
        %q = quantile(-Z,0.99);
        %test = -Z < q;
        %Z = test.*Z - (1-test).*q;
	    stdZ = std(Z);
        Z = Z ./stdZ;
        Zs = Z(2:end,:);
		Zn = Z(1,:);
     
        % loop over bandwidths: they go from 0.1 to 10, quadratically
        for i = 1:nbw
                bandwidths(i,:) = 0.1 + (10-0.1)*((i-1)/(nbw-1))^2;
        endfor        
        for bwiter = 1:nbw
            bandwidth = bandwidths(bwiter,:);        
            % now the fit using mean and mediani
            weight = __kernel_normal((Zs-Zn)/bandwidth);
            weight = weight/sum(weight(:));
            % the nonparametric fits, use local linear
            r = LocalPolynomial(thetas, Zs, Zn,  weight, false, 1);
            thetahat = r.mean';
            thetahat = keep_in_support(thetahat); % chop off estimates that go out of support (rare, but happens)
            % now CIs
            weight = __kernel_normal((Zs-Zn)/bandwidth);
            weight = weight/sum(weight(:));
            % the nonparametric fits, use local linear
            %r = LocalPolynomial(thetas, Zs, Zn,  weight, true, 1);
            r = LocalConstant(thetas, weight, true);

            % CI coverage
            in10 = ((theta0 > r.c') & (theta0 < r.d'));
            in_ci(:,bwiter) = in_ci(:,bwiter) + in10;
            errors(rep,:,bwiter) = thetahat'- theta0';
            rmse = zeros(9,1);
            if rep > 1
                printf("bandwidth = %f\n", bandwidth);    
                contrib = errors(1:rep,:,bwiter);
                m = mean(contrib);
                s = std(contrib);
                e = contrib;
                b = mean(e);
                e = e.^2;
                mse = mean(e);
                rmse = sqrt(mse);
                clabels = char("bias","rmse");
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
                printf("\n\nEstimation results (LL mean): rep %d\n", rep);
                prettyprint([b ; rmse]', rlabels, clabels);
                printf("\n\n");
                toc;
            endif
            rmses(:,bwiter) = rmse;
            cicoverage(:,bwiter) = in_ci(:,bwiter)/rep;
        endfor   
        save(outfile, "bandwidths", "cicoverage", "rmses");
    endif
endfor
        
if not(MPI_Finalized) MPI_Finalize; endif

