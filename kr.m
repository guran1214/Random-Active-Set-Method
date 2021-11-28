function [ x, s, k, Aopt, avgI] = kr( Q, q, b, cA, A, tol, sil)

% solves: min q'x + (1/2) x'Qx  subject to: x  <= b
% Q is assumed to be symmetric positive definite
% KKT:    Qx + s + q = 0;  x <= b, s >=0, s'(x-b) = 0
% input:  (Q,q,b) problem data,
%           cA=0 ... initial active set empty set
%           cA=1 ... initial acitve set full set
%           cA=2 ... initial active set defined by A
%         A ... guess on initial active set
%         tol ...  allow small tol to be also able to deal with data lacking strict complementarity
%         sil ... if 0 in each iteration some output is shown
% output: (x,s) optimal solution, k= # iterations
%         Aopt = active set at opt. sol. 
%         avgI ... average size of the linear systems solved
% call:   [ x, s, iter, Aopt] = kr( Q, d, b, A, 10^-10, 0);
  
n = length(q);                  % problem size
k = 0;                           
done = 0;                        % not yet done                     %
avgI = 0;

if(cA == 1) % as initial active set is full, just compute the dual variables and define new active set
    s = -q - Q*b;
    A = find(s>0);
    k = 1;
else if(cA == 0)
      A = [];
    end
end

% main loop 
while done < 1;                  % while not optimal
     k = k + 1;                  % start a new iteration
% solve system KKT(A):
     I = ones(n,1);              % compute I, the complement of A
     I(A) = zeros( length(A),1); 
     I = find(I>0);              % complement of A
     avgI = avgI + length(I);
     x = b;                      % x( A) = b( A); 
     I = sort(I);                % sorted I speeds up the computations
     s = zeros( n,1);            % s( I) = 0
     if ~isempty( I);          
          rhs = -q( I);
         if ~isempty(A);       % update right hand side rhs
           rhs = rhs - Q(I,A)*x( A);   
         end; 	   
     x(I) = Q(I,I) \ rhs;      % solve for inactive variables	   
     end;                               
     if ~isempty( A);            % backsubsitute for s(A), if |A|>0 
         s( A) = -q( A) - Q( A, :)* x;  
     end;
     if(sil == 0)   % show some information on the progress on the screen
         f2alt = 0.5*x'*Q*x+q'*x + 0.5*max(x-b,0)'*Q*max(x-b,0);
         err = log10(max( norm(Q*x+q+s,1), 1e-99));
         err_p = log10(max( max(x-b),1e-99));
         err_d = log10(max( -min(s), 1e-99));
         fprintf(' %3.0d %6.0f %6.0f %10.2f %10.2f %10.2f %10.10f \n',[k  (n-length(A)) length(A) err_p err_d err f2alt]);
     end
     A = [find(x>b); find(s>0)]; % active constraints of new solution
     done = (max(x-b)<= tol) & (min(s)>= -tol); % if I1 is empty the new point is also dual feasible and hence optimal -
                         % we allow for a small tolerance tol to be also
                         % able to deal with not strictly complementary
                         % data - our default value for tol is 10^-10
     if k > 200; done = 1;       % emergency exit to if kr-algorithm cycles
        if(sil == 0) 
         disp(' max number of iterations reached.');
        end
     end;
end;                 % end while
Aopt = A;
avgI = avgI/k;












