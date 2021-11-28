function [x, s, Aopt, time, k, maxfloor, totalsolves, avgI, objnew, I1, Iopt, floor, tol, sil] = mkr( Q, q, b, cA, A, I1, I2, floor, time, maxfloor, totalsolves, avgI, tol, sil)

% call: [x,alpha,Aopt,time,iter,objvalue] = hr_performance(Q,q,b)
% default starting active set is full
% enter data vectors q and b as column vectors
% solves: min J(x) = q'x + (1/2) x'Qx  subject to: x  <= b
% Q is assumed to be symmetric positive definite
% KKT:    Qx + s + q = 0;  x <= b, s >=0, s'(x-b) = 0
% input:    (Q,q,b) problem data,
%           cA=0 ... initial active set empty set
%           cA=1 ... initial acitve set full set
%           cA=2 ... initial active set defined by A
%           A ... guess on initial active set,
%           tol ...  allow small tol to be also able to deal with data lacking strict complementarity
%           sil ... if 0 in each iteration some output is shown
%           I1, I2, floor, time, maxfloor, totalsolves, avgI ... variables needed in the recursive calls
% output:   (x,s) optimal solution,
%           Aopt = active set at opt. sol. 
%           time ... time in seconds
%           k ... # iterations
%           maxfloor ... maximal depth of levels of subproblems (recursive calls)
%           totalsolves ... totat solves of a linear system on the inactive variables (computational bootleneck)
%           avgI ... average size of the linear systems solved
%           objnew ... optimal objective value
%           I1, Iopt, floor, sil, tol ... variables needed in the recursive calls

n = length(q);                  % problem size
done = 0;                        % not yet done
k = 0;
nosub = 0;
if(nargin == 3)
    cA = 1; A = []; I1 = []; I2 = []; floor= 0; time = 0; tol = 10^-10; sil = 1;
end

if (cA == 1 & floor == 0) % as initial active set is full, just compute the dual variables and define new active set
    tic;
    s = -q - Q*b;
    A = find(s>=0);
    I1 = ones(n,1);
    I1(A) = 0;
    I1 = find(I1>0);
    I2 = [];
    k = k+1;
    totalsolves = 0;
    maxfloor = 0;
    avgI = 0;
    else if (cA == 0 & floor == 0) % initial active set is empty
        tic;
        A = [];
        I1 = 1:n;
        I2 = [];
        totalsolves = 0;
        maxfloor = 0;
        avgI = 0;
        else if(cA == 2 & floor == 0) % initial active set is given by the input data
             tic;
             I1 = ones(n,1);
             I1(A) = 0;
             I1 = find(I1>0);
             I2 = [];
             totalsolves = 0;
             maxfloor = 0;
             avgI = 0;
            end
        end
end

I = [I1;I2];
N = [A;I];
N = sort(N);  % sorted N speeds up the computations
objnew = 10^20;

% main loop 
while done < 1;                  % while not optimal
     k = k + 1;                  % start a new iteration
     innerdone = 0;
     Aold = A;
     A1sum = [];
     A2sum = [];
     innerit = 0;
     objold = objnew;
% solve system KKT(A):
    while innerdone < 1
     innerit = innerit + 1;
     if innerit > 1             % update inactive sets
         help = zeros(n,1);
         help(I1) = 1;
         help(I1check) = 0;
         I1 = find(help>0);
         help = zeros(n,1);
         help(I2) = 1;
         help(I2check) = 0;
         I2 = find(help>0);
     end
     I = [I1;I2];
     avgI = avgI + length(I);
     I = sort(I);                % sorted I speeds up the computations
     x = b;                      % x(A) = b(A); 
     if isempty(I) == 0;         % compute primal variables on the inactive set
         rhs = -q(I);
         if isempty(A) == 0;       % update right hand side rhs
           rhs = rhs - Q(I,A)*x( A);   
         end; 	   
         x(I) = Q(I,I) \ rhs;   % solve for inactive variables	
     end;
         innerdone = (max(x-b) <= 0);   % exit inner loop if primal feasible x has been found
         help1 = x>=b;                  % update active set: all x >= b
         help2 = zeros(n,1);
         help2(I1) = 1;
         I1check = help1.*help2;
         I1check = find(I1check>0);
         A1sum = [A1sum;I1check];
         help2 = zeros(n,1);
         help2(I2) = 1;
         I2check = help1.*help2;
         I2check = find(I2check>0);
         A2sum = [A2sum;I2check];
         Anew = [I1check;I2check];
         A = [A;Anew];
    end
    
     I2 = [I2;I1];
     s = zeros( n,1);  
     if isempty(A) == 0;            % backsubsitute for s(A), if |A|>0 
         s(A) = -q(A) - Q(A,N)* x(N);
     end;
     
     
     objnew = 0.5*x(N)'*Q(N,N)*x(N)+q(N)'*x(N);   % compute objective value of the new point
     
     if(nosub == 0) % if number of constraints is reduced then "nosub = 1" as we start with a "new" problem
     if (objnew >= objold)  % solve a smaller subproblem if the objective value did not improve
         help1 = s < 0;     % define active and inactive set for the subproblem
         help2 = zeros(n,1); % take all variables expect the ones in Aold to the subproblem
         help2(A2sum) = 1;
         help2(A1sum) = 1;
         Add = [];
         if(isempty(Aold)) % if Aold (in paper B_s) is empty take one element from A1sum (in paper B_1) 
             % or if A1sum is also empty then from I_1 which is the last part of I2. 
             % Note that |I_1 \cup B_1| > 1: \not = \emptyset -> see Lemma 6 from paper; \not = 1 
             % because then nosub = 1 (one bound removed) -> see below
             if(isempty(A1sum))
                Add = I2(end);
                I2 = I2(1:end-1);
             else
                Add = A1sum(1);
                help2(A1sum(1)) = 0;
             end
         end
         Aleft = [Aold,Add];
         I1inner = help1.*help2;
         I1inner = find(I1inner>0);
         help1 =  s>= 0;
         Ainner = help1.*help2;
         Ainner = find(Ainner>0);
             floor = floor + 1;
             [x, s, A2sum, time, ki, maxfloor, totalsolves, avgI, objnew, I1, I2, floor, tol, sil] = mkr( Q, q+Q(:,Aleft)*b(Aleft),b,cA, Ainner, I1inner, I2, floor, time, maxfloor, totalsolves, avgI, tol, sil);
             floor = floor - 1;
             A = [Aleft;A2sum];
             I = I2; % I1 is empty by construction
             I = sort(I);
             x = b; % recompute x und s because d is different in the subproblem
             if isempty(I) == 0;          
                 rhs = -q(I);
                 if isempty(A) == 0;       % update right hand side rhs
                   rhs = rhs - Q( I,A)*x(A);   
                 end
                 x( I) = Q(I,I) \  rhs;   % solve for inactive variables
             end;
              s = zeros(n,1);
             if isempty(A) == 0            % backsubsitute for s(A), if |A|>0 
                s(A) = -q(A) - Q(A,N)*x(N);  
             end;
     objnew = 0.5*x(N)'*Q(N,N)*x(N)+q(N)'*x(N); % update objective value for new primal feasible point
     end
     end
     
     Apre = A;          %define new active and inactive sets
     I1 = find(s<0);
     I2 = find(x<b);
     A = zeros(n,1);
     A(N) = 1;
     A(I1) = 0;
     A(I2) = 0;
     A = find(A>0);
     nosub = 0;
     
     if (sil == 0)  % show some information on the progress on the screen
         err = log10(max( norm(Q(N,N)*x(N)+q(N)+s(N),1), 1e-99));
         err_p = log10(max( max(x-b),1e-99));
         err_d = log10(max( -min(s), 1e-99));
         fprintf(' %3.0d %3.0d %6.0f  %6.0f  %6.0f  %10.2f %10.2f %10.2f  %10.10f  \n',[k innerit length(I1)+length(I2)   length(A)  floor   err_p   err_d   err  objnew]);
     end
     
     totalsolves = totalsolves + innerit;
     
    if(isempty(A) == 1 && length(Apre)==1) % in this special case, one bound can be removed
        b(Apre) = 10^20;
        nosub = 1;  % "new problem" - in the first iteration no subproblem is allowed
    end
    
    done = ~any(s<-tol); % if I1 is empty the new point is also dual feasible and hence optimal -
                         % we allow for a small tolerance tol to be also
                         % able to deal with not strictly complementary
                         % data - our default value for tol is 10^-10
end;                 % end while
 
if(floor == 0)
    time = toc; %if finished, stop time
    if(cA == 1) % if finished compute average size of linear system
        avgI = avgI/(totalsolves+1);
    else
        avgI = avgI/totalsolves;
    end
end

    Aopt = A;
    Iopt = I2;
    maxfloor = max(maxfloor,floor);
%     A