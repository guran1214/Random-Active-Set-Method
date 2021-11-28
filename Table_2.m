function [] = Table_2()

% set the number of variables "n", the number of experiments "runs" to the
% values proposed in the paper and the solve the resulting problems for
% different values of epsilon (pot) and for both empty (cA=0) and full
% (cA=1) initial active set


n = 2000; 
runs = 10;
cA = [1];
pot = [0,5,10,14];

for i = 1:1
    
    cA_act = cA(i);
    
    for j = 1:4
        
        pot_act = pot(j);
        
        tridiag_data(runs,n,pot_act,cA_act);

    end
end
