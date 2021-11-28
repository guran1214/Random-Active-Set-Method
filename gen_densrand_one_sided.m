function [] = gen_sprand_one_sided_rg(n,density,runs,cond)

% initialize variables for output
anoi_kr = 0;
atime_kr = 0;
inact_var_kr = 0;
cyc_kr = 0;
anoi_hr = 0;
atime_hr = 0;
anoi_mat = 0;
atime_mat = 0;
inact_var_hr = 0;
allmaxfloor = 0;
atotal = 0;
arel_err = 0;
atime_rg = 0;
anoi_rg = 0;
ainact_rg = 0;
fprintf('\n Parameters \n');
fprintf(' n    runs  density cond \n');
fprintf('%d  %d    %3.3f  %1.1e\n',n,runs,density,cond);

for i = 1:runs
    
    % generate data
    
    g = rand(n,1)-0.5;
    O=randn(n); 
    [O,~]=qr(O);
    p=1; D=diag(cond.^(((0:n-1)/(n-1)).^p)); Q=O*D*O'; Q=(Q+Q')/2;
    q=-g; b=zeros(n,1);
    
    % call KR-Algorithm
    solvetime = tic;
    %fprintf('\n\n KR-Algorithm\n');
    [x,alpha,iter,Aopt,avg_inact] = kr( Q, q, b, 1, [], 10^-10, 1);
    kr_time = toc(solvetime);
    if (iter > 200)
            cyc_kr = cyc_kr + 1;
    else
            anoi_kr = anoi_kr + iter;
            atime_kr = atime_kr + kr_time;
            inact_var_kr = inact_var_kr + avg_inact;
    end

    % call mKR-Algorithm
    solvetime = tic;
    %fprintf('\n\n Modified KR-Algorithm\n');
    [x, alpha, Aopt, dummy, iter, maxfloor, totalsolves, avg_inact] = mkr(Q, q, b, 1, [], [], [], 0, 0, 0, 0, 0, 10^-10, 1);
    hr_time = toc(solvetime);
    anoi_hr = anoi_hr + iter;
    atotal = atotal + totalsolves;
    atime_hr = atime_hr + hr_time;
    inact_var_hr = inact_var_hr + avg_inact;
    allmaxfloor = max(allmaxfloor,maxfloor);
    obj_opt = 0.5*x'*Q*x+q'*x;
    %fprintf('Error Objective: %10.10f:\n\n',abs(0.5*x'*Q*x+q'*x-obj_opt));

    % call random activeset algorithm
    outertime = tic;
    [ ~,~,~,iter,avgI,exitflag] = ras( Q,-(Q*b+q),false(size(b)),1e-10,0);
    time = toc(outertime);
    anoi_rg = anoi_rg + iter;
    atime_rg = atime_rg + time;    
    ainact_rg = ainact_rg+avgI;
    
    % call MATLAB quadprog interior point method
    %fprintf('\n\n MATLAB quadprog interior point method\n');
    options = optimset('quadprog');
    options = optimset(options,'Display','off'); 
     solvetime = tic;
    [x,dummy1,dummy2,output] = quadprog(Q, q, [], [], [], [], [], b, [], options);
    anoi_mat = anoi_mat + output.iterations;
    matlab_time = toc(solvetime);
    atime_mat = atime_mat + matlab_time;
    arel_err = arel_err + (0.5*x'*Q*x+q'*x-obj_opt)/abs(obj_opt);
end

anoi_hr = anoi_hr/runs;
atime_hr = atime_hr/runs;
inact_var_hr = inact_var_hr/runs;
atotal = atotal/runs;

anoi_kr = anoi_kr/(runs-cyc_kr);
atime_kr = atime_kr/(runs-cyc_kr);
inact_var_kr = inact_var_kr/(runs-cyc_kr);

anoi_rg = anoi_rg/(runs);
atime_rg = atime_rg/(runs);
ainact_rg = ainact_rg/(runs);

anoi_mat = anoi_mat/runs;
atime_mat = atime_mat/runs;
arel_err = arel_err/runs;

fprintf('\n\n Kunisch-Rendl-Algorithm\n'); % display some info on screen 
fprintf(' atime   anoi   kr_cyc  inact_var\n');        
fprintf(' %3.3f   %3.3f   %4d    %3.3f \n\n',atime_kr,anoi_kr,cyc_kr,inact_var_kr);

fprintf('\n\n Modified Kunisch-Rendl-Algorithm\n'); % display some info on screen 
fprintf('  atime   anoi   total_solves   mds   inact_var\n');        
fprintf(' %3.3f   %3.3f   %3.3f    %4d   %3.3f \n\n',atime_hr,anoi_hr,atotal,allmaxfloor,inact_var_hr);

fprintf('\n\n Random Active Set Method\n'); % display some info on screen 
fprintf(' atime   anoi  inact_var\n');        % display some info on screen 
fprintf(' %3.3f   %3.3f    %3.3f \n\n',atime_rg,anoi_rg,ainact_rg);

fprintf('\n\n MATLAB quadprog interior point method\n'); % display some info on screen 
fprintf('  atime   anoi   relative error\n');        
fprintf(' %3.3f   %3.3f   %10.10e \n\n',atime_mat,anoi_mat,arel_err);
