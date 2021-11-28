function [] = Table_4()

% set the number of variables "n", the number of experiments "runs" to the
% values proposed in the paper and the solve the resulting problems for
% different densities and condition numbers

n = [500,1000,2000,4000];
density = [1];
runs = 10;
cond = [10^6,10^10,10^14];

% problems with one-sided bounds and densities 0.1 and 0.01
fprintf('One-sided case\n\n');
for i = 1:4
    
    n_act = n(i);
    
    for j = 1:1
        
        dens_act = density(j);
        
        for k = 1:3
            
            cond_act = cond(k);
        
            gen_densrand_one_sided(n_act,dens_act,runs,cond_act);
            
        end

    end
end
end