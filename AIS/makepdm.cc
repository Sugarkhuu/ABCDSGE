#include <oct.h>
#include <octave/parse.h>

DEFUN_DLD(makepdm, args, ,"makepdm: parameter dependent moments")
{
	int nargin = args.length();
	Matrix thetas (args(0).matrix_value());
	Matrix data (args(1).matrix_value());

	int nn = thetas.rows();
    int T = data.rows();
	Matrix pdm(nn,8);
    pdm.fill(0.0);
	octave_value_list f_return;
    
    int i, t;
    double alpha, beta, delta, gam, rho_z, sig_z, rho_eta, sig_eta, nss;
    double c1, kss, iss, yss, css, psi;
    double y, lagy, c, lagc, n, lagn, r, lagr, w, lagw;
    double lneta, laglneta, e;
    double lag2n, lag2r, lag2w, k, lagk, lnz, laglnz, invest, laginvest;

	for (i = 0; i < nn; i++) {
            // break out params    
            alpha = thetas(i,0);
            beta = thetas(i,1);
            delta = thetas(i,2);
            gam = thetas(i,3);
            rho_z = thetas(i,4);
            sig_z = thetas(i,5);
            rho_eta = thetas(i,6);
            sig_eta = thetas(i,7);
            nss = thetas(i,8);
            // recover psi
            c1 = pow(((1.0/beta + delta - 1.0)/alpha),(1.0/(1.0-alpha)));
            kss = nss/c1;
            iss = delta*kss;
            yss = pow(kss,alpha) * pow(nss,(1.0-alpha));
            css = yss - iss;
            psi =  pow(css,-gam) * (1.0-alpha) * pow(kss,alpha) * pow(nss,-alpha);

            // loop over observations: averages are from t=3 to T,
            // because we loose one when computing k, and we use
            // a lag of computed k
            for (t = 2; t < T; t++) {
                // the needed variables and lags
                y = data(t,0);
                lagy = data(t-1,0);
                c = data(t,1);
                lagc = data(t-1,1);
                n = data(t,2);
                lagn = data(t-1,2);
                lag2n = data(t-2,2);
                r = data(t,3);
                lagr = data(t-1,3);
                lag2r = data(t-2,3);
                w = data(t,4);
                lagw = data(t-1,4);
                lag2w = data(t-2,4);

                // MRS-MUC-MUL
                lneta = log(w) - gam*log(c) - log(psi);
                laglneta = log(lagw) - gam*log(lagc) - log(psi);
                e = lneta - rho_eta*laglneta;
                e = e/sig_eta;
                pdm(i,0) += e;
                pdm(i,1) += e*e - 1.0;
                // now the Euler eqn
                e = beta*pow(c,-gam)*(1.0 + r -delta) - pow(lagc,-gam);
                pdm(i,2) += e;
                // production function
                k = alpha*lagn * lagw / lagr / (1.0-alpha);
                lagk = alpha*lag2n * lag2w / lag2r / (1.0-alpha);
                lnz = log(y) - alpha*log(k) - (1.0-alpha)*log(n);
                laglnz = log(lagy) - alpha*log(lagk) - (1.0-alpha)*log(lagn);
                e = lnz - rho_z*laglnz;
                e = e/sig_z;
                pdm(i,3) += e;
                pdm(i,4) += e*e - 1.0;
                // MPL
                lnz = log(w) + alpha*(log(n)-log(k)) - log(1.0-alpha);
                laglnz = log(lagw) + alpha*(log(lagn)-log(lagk)) - log(1.0-alpha);
                e = lnz - rho_z*laglnz;
                e = e/sig_z;
                pdm(i,5) += e;
                pdm(i,6) += e*e - 1.0;
                // law of motion k: good for delta
                invest = y - c;
                laginvest = lagy - lagc;
                e = laginvest + (1.0 - delta)*lagk - k;
                pdm(i,7) += e;
            }
    }    
    pdm = pdm /(double(T)-2.0);        
	f_return(0) = pdm;
	return f_return;
}

