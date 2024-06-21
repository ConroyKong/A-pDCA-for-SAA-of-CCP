function [x, time, iter] = SCA(S, U, alpha, t, opts)

fprintf('***************** SCA for DC approximation ***************** \n');

%% parameter setting
[N, n, m] = size(S);
maxiter = 2*1e2;

if isfield(opts,'solver')
    solver = opts.solver;
else
    solver = 'gurobi';
end
if isfield(opts,'maxitime')
    maxitime = opts.maxitime;
else
    maxitime = 1800;
end
if isfield(opts,'x0')
    x0 = opts.x0;
else
    x0 = randn(n,1);
end
if isfield(opts, 'tol')
    tol = opts.tol;
else
    tol = 1e-6;
end

xk = x0;
fval = -sum(x0);
time = 0;

for iter = 1 : maxiter
    fval_old = fval;
    %% compute subgradient of \sum_{i=1}^N \max{C(x,\hat{\xi}^i),0} at xk
    for j = 1:m
        C(:,j) = (S(:,:,j).^2)*(xk.^2) - U;
    end
    [CN, inx_set]= max(C,[],2);
    %% compute subgradient of H(x) at xk
    gk = zeros(n,1);
    gval = 0;
    for i = 1:N
        if CN(i) >= 0
            inx = inx_set(i);
            sk = 2*(S(i,:,inx).^2)'.*xk;
            gk = gk + sk; gval = gval + CN(i);
        end
    end

    if strcmp(solver, 'gurobi') == 1
        %% Model parameters in GUROBI
        model.A = sparse([-gk' ones(1,N)]);
        model.rhs = N*alpha*t+gval-gk'*xk;
        model.lb = [zeros(1,n) zeros(1,N)];
        model.ub = inf(1,n+N);
        % quad terms
        counter = 1;
        for j = 1:m
            S_j_squared = S(:,:,j).^2;
            S_j_squared = round(S_j_squared, 6);
            for i = 1:N
                model.quadcon(counter).Qc = sparse(1:(n + N), 1:(n + N), [S_j_squared(i, :), zeros(1, N)]);
                model.quadcon(counter).q  = sparse(n + N, 1); model.quadcon(counter).q(n + i) = -1;
                model.quadcon(counter).rhs = U - t;
                counter = counter + 1;
            end
        end
        model.obj = [-ones(1, n) zeros(1,N)];
        model.modelsense = 'Min';
        model.sense = '<';

        %% Solve the problem using GUROBI
        params.outputflag = 0;
        params.timelimit = maxitime;
        params.optimalitytol = tol;
        result = gurobi(model, params);
        xk = result.x(1:n);
        time = time + result.runtime;
    else
        %% define variables
        x = sdpvar(n,1); z = sdpvar(N,1);

        %% intial point
        assign(x, x0);

        %% Define constraints
        Constraints = [sum(z) - (gval + gk'*(x-xk)) <= N*alpha*t, x >= zeros(n,1),...
            z >= zeros(N,1)];
        for j = 1:m
            Constraints = [Constraints, (S(:,:,j).^2)*(x.^2) - U*ones(N,1) + t*ones(N,1)  <= z];
        end

        %% Define an objective
        Objective = -sum(x);

        %% solve by GUROBI solver
        ops = sdpsettings('solver', solver, 'verbose', 0, 'usex0', 1, 'gurobi.timelimit', maxitime);
        sol = optimize(Constraints, Objective, ops);
        time = time + sol.solvertime;

        %% report iterate information
        x = value(x);
    end

    fval = -sum(xk);
    fprintf('Iternum: %d, fval: %.4f\n', iter, fval);

    if abs(fval - fval_old)/max(1,abs(fval_old)) <= tol || time > 1800
        break;
    end

end
x = xk;
end