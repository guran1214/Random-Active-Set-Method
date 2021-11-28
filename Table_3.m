function [] = Table_3()

% set the number of variables "n", the number of experiments "runs" to the
% values proposed in the paper and the solve the resulting problems for
% different densities and condition numbers

n = [1000,5000,10000];
density = [0.1,0.01,0.001];
runs = 10;
cond = [10^2,10^6,10^10,10^14];

% problems with one-sided bounds and densities 0.1 and 0.01
fprintf('One-sided case\n\n');
for i = 1:2
    
    n_act = n(i);
    
    for j = 1:2
        
        dens_act = density(j);
        
        for k = 1:4
            
            cond_act = cond(k);
        
            gen_sprand_one_sided(n_act,dens_act,runs,cond_act);
            
        end

    end
end

% problems with one-sided bounds and density 0.001
for i = 2:3
    
    n_act = n(i);
    
    dens_act = density(3);
        
    for k = 1:4

        cond_act = cond(k);

        gen_sprand_one_sided(n_act,dens_act,runs,cond_act);

    end
end



end