function [] = tridiag_data(runs,n,pot,cA)
% n ... number of variables
% pot ... factor for addition of I
% cA=0 ... initial active set empty set
% cA=1 ... initial acitve set full set

atime_hr = 0;
anoi_hr = 0;
ainact_hr = 0;
totalmax = 0;
atime_kr = 0;
anoi_kr = 0;
ainact_kr = 0;
asolves = 0;
fail = 0;
atime_rg = 0;
anoi_rg = 0;
ainact_rg = 0;
acss = 'not defined';
if(cA==0)
    acss = 'empty';
end
if(cA==1)
    acss = 'full';
end

for i = 1:runs
    
    % generate random data
    rand('seed',n*10000+i*1000);
    randn('seed',n*10000+i*1000);
    p = sprandn(n,n,.1)+speye(n);
    p = tril(triu(p,-100));
    Q0 = p*p';
    q=rand(n,1)*20*n-10*n;
    b=ones(n,1);
    Q = Q0 + 10^-pot*eye(n);
    Q = sparse(Q);
    
    % call KR-Algorithm
    outertime = tic;
    [xopt, sopt, oiter, Aopt, avg_inact] = kr(Q, q, b, cA, [], 1e-8, 1);
    time = toc(outertime);
    if (oiter > 200)
        fail = fail + 1;
    else
            anoi_kr = anoi_kr + oiter;
            ainact_kr = ainact_kr + avg_inact;
            atime_kr = atime_kr + time;
    end
    
    % call mKR-Algorithm
    outertime = tic;
    [xopt, sopt, Aopt, dummy, oiter, maxfloor, totalsolves, avginact] = mkr( Q, q, b, cA, [], [], [], 0, 0, 0, 0, 0, 1e-8, 1);
    time = toc(outertime);
    anoi_hr = anoi_hr + oiter;
    atime_hr = atime_hr + time;    
    ainact_hr = ainact_hr+avginact;
    asolves = asolves + totalsolves;
    totalmax = max(totalmax,maxfloor);

    % call random active set method
    outertime = tic;
    [ ~,~,~,iter,avgI,exitflag] = ras( Q,-(Q*b+q),zeros(size(b))==cA,1e-8,0);
    time = toc(outertime);
    anoi_rg = anoi_rg + iter;
    atime_rg = atime_rg + time;    
    ainact_rg = ainact_rg+avgI;
end

anoi_hr = anoi_hr/runs;
atime_hr = atime_hr/runs;
ainact_hr = ainact_hr/runs;
asolves = asolves/runs;

anoi_kr = anoi_kr/(runs-fail);
atime_kr = atime_kr/(runs-fail);
ainact_kr = ainact_kr/(runs-fail);

anoi_rg = anoi_rg/(runs);
atime_rg = atime_rg/(runs);
ainact_rg = ainact_rg/(runs);

fprintf('\n\n Kunisch-Rendl-Algorithm\n'); % display some info on screen 
disp('niter epsilon   atime   anoi  fail  acss   average inact');        % display some info on screen 
fprintf('%4d   %3.2f   %3.2f    %3.2f  %d  %s   %3.2f\n\n',runs,pot,atime_kr,anoi_kr,fail,acss,ainact_kr);


fprintf('\n\n Modified Kunisch-Rendl-Algorithm\n'); % display some info on screen 
disp('niter epsilon  atime anoi  total_solves  mds   acss   average inact');
fprintf('%4d    %3.2f   %3.2f  %3.2f     %3.2f   %d    %s      %3.2f \n\n',runs,pot,atime_hr,anoi_hr,asolves,totalmax,acss,ainact_hr);

fprintf('\n\n Random Active Set method\n'); % display some info on screen 
disp('niter epsilon   atime   anoi  acss   average inact');        % display some info on screen 
fprintf('%4d   %3.2f   %3.2f    %3.2f  %s   %3.2f\n\n',runs,pot,atime_rg,anoi_rg,acss,ainact_rg);
