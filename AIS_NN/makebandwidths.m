% create the bandwidths used in tuning
nbw = 20;
bandwidths = zeros(nbw,1);
for i = 1:nbw
        bandwidths(i,:) = 0.001 + (10-0.001)*((i-1)/(nbw-1))^2;
endfor 
