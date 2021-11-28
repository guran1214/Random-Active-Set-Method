function [ x,s,A,iter,avgI,exitflag] = ras( Q,g,I,tol,sil,maxiter)
% Random Active Set Method
% solve: min 0.5x'Qx + g'x, s.t. x>=0
% enter data vector b as column vector
% Q is assumed to be symmetric positive definite
% KKT:    Qx + g - s = 0;  x >= 0, s >=0, x's = 0
% input:    (Q,g) problem data,
%           I ... guess on initial inactive set, which is a logistic(binary) vector, and the indexes of ones stand for the inactive set
%           tol ... tolerance on nonnegativity violation of s
%           sil ... if 2 in each iteration some output is shown, if 1 output is shown at the end, if 0 output is not allowed to be shown
%           maxiter ... maximum value of the number of iterations
% output:   (x,s) optimal solution,
%           I ... inactive set at optimal solution x, which is a logistic(binary) vector
%           iter ... number of iterations
%           avgI ... average size of the linear systems solved
%           exitflag ... if 1 problem is solved succesfully, if 0 it fails to solving in the given maxiter iterations
% coded by Ran Gu (Nankai Universtiy)
% 2021-11-07
%% initialization
n = length(g); % n is the dimension
if nargin <= 5 || isempty(maxiter)
    maxiter = 200; 
end
if nargin <= 4 || isempty(sil)
    sil = 1; % default: only print the final information
end
if nargin <= 3 || isempty(tol)
    tol = 1e-10;
end
if nargin <= 2 || isempty(I)
    I = false(n,1);
end
if any(I>1) && all(I-floor(I)==0)
    I0 = I; I = false(n,1); I(I0) = true;
end
A = ~I; % A stands for the active set, which is also a logistic vector
sumI = 0; % count the size of solved linear equations for output information
exitflag = 0; % if 1 problem is solved succesfully, if 0 it fails to solving in the given maxiter iterations
nzeros = zeros(n,1); nfalse = false(n,1);
%% when iter = 1
xI = -Q(I,I)\g(I); % solve KKT equalities at I
sA = Q(A,I) * xI + g(A); % solve KKT equalities at A
x = nzeros; x(I) = xI; s = nzeros; s(A) = sA; % construct x and s
iter = 1; % the number of sovling linear equations
sumI = sumI + nnz(I); % sum up the sizes of linear equations solved
Im = x < 0; % Im is a logistic vector standing for the negative indexes of x
Ip = I & (~Im); % Ip is a logistic vector standing for the non-negative indexes of x
Am = s < -tol; % Am is a logistic vector standing for the negative indexes of s
Ap = A & (~Am); % Ap is a logistic vector standing for the non-negative indexes of s
if sil==2
    fprintf('iter       I       A      Im     Am  \n')
    fprintf(' %3.0d  %6.0f  %6.0f  %6.0f %6.0f\n',iter,nnz(I),nnz(A),nnz(Im),nnz(Am))
end
if ~(any(Im)|| any(Am)) % stopping criteria: KKT system solved
    exitflag = 1;
    if sil > 0
        fprintf('problem is solved in %i iterations\n',iter);
    end
    avgI = sumI / iter;
    return;
end
NImp0 = nfalse; NAmp0 = nfalse; NImf = Im; NAmf = Am; NImc = nfalse; NAmc=nfalse;
Ip0 = Ip; Ap0 = Ap;

%% when iter > 1

while iter <=maxiter
    for loop = 1:1000 % find a different active set
        Imf = (NImp0&(rand(n,1)<=0.5))| (NImf&(rand(n,1)<=0.02))| (NImc&(rand(n,1)<=0.02));
        Amf = (NAmp0&(rand(n,1)<=0.99))| (NAmf&(rand(n,1)<=0.07))| (NAmc&(rand(n,1)<=0.04));
        Imc = Im&(~Imf); Amc = Am&(~Amf);
        if any(Imc)||any(Amc)
            I = Ip|Imf|Amc;
            A = ~I;
            break;
        else
            NImf = NImf|NImp0|NImc; NAmf = NAmf|NAmp0|NAmc; NImp0 = nfalse; NAmp0 = nfalse; NImc = nfalse; NAmc = nfalse;
        end
    end
    xI = -Q(I,I)\g(I); % solve KKT equalities at I
    sA = Q(A,I) * xI + g(A); % solve KKT equalities at A
    x = nzeros; x(I) = xI; s = nzeros; s(A) = sA; % construct x and s
    iter = iter + 1; % the number of sovling linear equations
    sumI = sumI + nnz(I); % sum up the sizes of linear equations solved
    Im = x < 0; % Im is a logistic vector standing for the negative indexes of x
    Ip = I & (~Im); % Ip is a logistic vector standing for the non-negative indexes of x
    Am = s < -tol; % Am is a logistic vector standing for the negative indexes of s
    Ap = A & (~Am); % Ap is a logistic vector standing for the non-negative indexes of s
    NImp0 = Im&Ip0; NImf = Im&Imf; NImc = Im&Amc; NAmp0 = Am&Ap0; NAmf = Am&Amf; NAmc = Am&Imc;
    
    if sil == 2 % show iterative information
        fprintf(' %3.0d  %6.0f  %6.0f  %6.0f %6.0f\n',iter,nnz(I),nnz(A),nnz(Im),nnz(Am))
    end
    if ~(any(Im)|| any(Am)) % stopping criteria: KKT system solved
        exitflag = 1;
        break;
    end
    Ip0 = Ip; Ap0 = Ap;
end
%% exit
avgI = sumI / iter;
if exitflag
    if sil > 0
        fprintf('problem is solved in %i iterations\n',iter);
    end
else
    if sil > 0
        fprintf('failed because no solution has been found in %i iterations\n',maxiter);
    end
end

end
