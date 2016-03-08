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
            lneta = log(w) - gam*log(c) -log(psi);
            e = lneta - rho_eta*lag(lneta,1);
            e = e/sig_eta;
            e = e(3:end,:);
            e1 = e;
            pdm(i,1) = mean(e);
            pdm(i,2) = mean(e.^2 - 1);
            %e1 = e;
            % now the Euler eqn
            cc = c ./lag(c,1);
            e = (beta*cc.^(-gam).*(1 + r -delta))-1;
            e = e(2:end,:);
            pdm(i,3) = mean(e);
            % production function
            k = alpha*lag(n,1).*lag(w,1)./lag(r,1)/(1-alpha);
            lnz = log(y) - alpha*log(k) - (1-alpha)*log(n);
            e = lnz - rho_z*lag(lnz,1);
            e = e(3:end,:);
            e = e/sig_z;
            pdm(i,4) = mean(e);
            pdm(i,5) = mean(e.^2 - 1);
            pdm(i,6) = mean(e1.*e); % no cross error corr.
            % MPL
            %lnz = alpha*log(r/alpha) + (1-alpha)*log((1-alpha)./w);
            %e = lnz - rho_z*lag(lnz,1);
            %e = e(2:end,:);
            %e = e/sig_z;
            %pdm(i,5) = mean(e);
            %pdm(i,6) = mean(e.^2 - 1);
            % law of motion k: good for delta
            invest = y - c;
            e = lag(invest,1) + (1 - delta)*lag(k,1) - k;
            e = e(3:end,:);
            pdm(i,7) = mean(e);
        endfor
        pdm = [zeros(1,size(pdm,2)); pdm]; % stack real data stats (zeros) on top
endfunction
