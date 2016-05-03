#{
How to use this:
Set DO_NN and DO_PRM as desired.
1 run DSGE_tune, setting DO_LOCAL = false
2 run DSGE_montecarlo
3 change DO_LOCAL to true in DSGE_tune and DSGE_montecarlo
4 repeat steps 1 and 2
5 create an appropriate output directory given DO_NN and DO_PDM,
  and move *.out to that dir.
6 repeat 1-5 for all the NN and PDM combos  
#}
DO_NN = true;
DO_PDM = false;
DO_LOCAL = true;

if !DO_LOCAL
    outfile = "tuned_from_prior.out";
else
    outfile = "tuned_local.out";
endif    
mc_reps = 100; % number of MC reps
nworkers = 25;  % number of worker MPI ranks

SetupAIS; # done in one place for uniformity

% design
parameters; % loaded from Common to ensure sync with Gendata

lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
prior_params = [lb ub];
theta0 = lb_param_ub(:,2); % original form
nparams = rows(theta0);

setupmpi; % sets comm world, nodes, node, etc.
asbil_theta = theta0; setupdynare; % sets structures and RNG for simulations
load selected;
asbil_selected = selected;
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
	thetahatsLC = zeros(mc_reps, nparams);
	cilower = zeros(mc_reps, nparams);
	ciupper = zeros(mc_reps, nparams);
endif


for rep = 1:mc_reps
    % the 'true' Zn
    if node==1 % simulate on node1 (can't do it on 0, no *.modfile
        asbil_theta = theta0;
        ok = false;
        while !ok    
            USERsimulation;
            Zn = aux_stat(data);
            ok = Zn(1,:) != -1000;
        endwhile
        realdata = data;
        if DO_NN
            Zn = NNstat(Zn');
        else
            Zn = Zn(asbil_selected,:)';    
        endif
        if DO_PDM
            pdm = makepdm(asbil_theta', realdata);
            Zn = [Zn pdm];
            n_pdm = size(pdm,2);
        else
            n_pdm = 0;    
        endif 
        for i = 2:nodes-1
            MPI_Send(Zn, i, mytag, CW);
            MPI_Send(realdata, i, mytag+1, CW);
        endfor	
        MPI_Send(Zn, 0, mytag, CW);
        MPI_Send(realdata, 0, mytag+1, CW);
    else % receive it on the other nodes
        Zn = MPI_Recv(1, mytag, CW);
        realdata = MPI_Recv(1, mytag+1, CW);
    endif
    if DO_PDM % nodes need to know the size
        pdm = makepdm(asbil_theta', realdata);
        n_pdm = size(pdm,2);
    else
        n_pdm = 0;    
    endif 
    MPI_Barrier(CW);    
	% call the algoritm that gets AIS particles
	if !node
            printf("starting AIS\n");
            tic;
    endif  
    CreateAIS; # this gets the particles

    % now draw from the AIS density
    reps_per_node = round(AISdraws/(nodes-1));
   	if !node
            toc;
            printf("starting AIS2\n");
            tic;
    endif  
    SampleFromAIS; # this samples from AIS density


	% see the results

   	if !node
        toc;    
        % create the bandwidths used in tuning
        makebandwidths;
        load SelectedBandwidths; # these come from most recent run of DSGE_tune.m
        bandwidthsCI = bandwidths(bwselectCI,:);
        bandwidths = bandwidths(bwselect,:);
	    printf("starting fit and CI\n");
		thetas = contribs(:,1:nparams);
        Zs = contribs(:, nparams+1:end);
        test = Zs(:,1) != -1000;
        thetas = thetas(test,:);
        Zs = Zs(test,:);
        Z = [Zn; Zs];
	    stdZ = std(Z);
        Z = Z ./stdZ;
        Zs = Z(2:end,:);
		Zn = Z(1,:);
        %AISweights = prior(thetass) ./ (mixture*prior(thetass) +(1-mixture)*AIS_density(thetass, particles(:,1:nparams)));
        AISweights = 1;
        weights = zeros(rows(Zs),9);
        for i = 1:9
            weights(:,i) = __kernel_normal((Zs-Zn)/bandwidths(i,:));
        endfor    
        weights = AISweights.*weights; # AIS_weights != 1 is for SBIL by AIS
        weights = weights./sum(weights);
        thetahatLC = zeros(1,9);
        for i = 1:9
                r = LocalConstant(thetas(:,i), weights(:,i), false);
                thetahatLC(:,i) = r.mean;
        endfor
        thetahatLC = keep_in_support(thetahatLC);
        % confidence intervals
        % weights
        %AISweights = prior(thetas) ./ (mixture*prior(thetas) + (1-mixture)*AIS_density(thetas, particles(:,1:nparams)));
        AISweights = 1;
        weights = zeros(rows(Zs),9);
        for i = 1:9
            weights(:,i) = __kernel_normal((Zs-Zn)/bandwidthsCI(i,:));
        endfor         
        weights = AISweights.*weights; # AIS_weights != 1 is for SBIL by AIS
        weights = weights./sum(weights);
        % for CIs, use local constant
        lower = zeros(9,1);
        upper = lower;
        for i = 1:9
            r = LocalConstant(thetas(:,i), weights(:,i), true);
            lower(i,:) = r.c;
            upper(i,:) = r.d;
        endfor    
		% results
        if rep == 1
                in_ci = zeros(rows(theta0),1);
        endif
        % CI coverage
        in10 = ((theta0 > lower) & (theta0 < upper));
        in_ci = in_ci + in10;
        thetahatsLC(rep,:) = thetahatLC;
        cilower(rep,:) = lower';
        ciupper(rep,:) = upper';
        save(outfile, "thetahatsLC", "cilower", "ciupper");
		if rep > 1
			contrib = thetahatsLC(1:rep,:);
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
			mae = mean(abs(e));
			clabels = char("true", "mean", "pmean", "sdev.","psdev", "bias", "pbias","rmse", "prmse");
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
			printf("\n\nEstimation results (LC mean): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
			printf("\n\n");
            printf("90%% CI coverage: \n");
            disp(in_ci/rep);
            disp([lower upper]);
            printf("\n");
		endif
    endif
endfor

if not(MPI_Finalized) MPI_Finalize; endif

