function [x, time] = CVaR(S, U, alpha, opts)

%     fprintf('***************** CVaR ***************** \n');

%% parameter setting
[N, n, m] = size(S);
% M = ceil((1-alpha)*N);

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

if strcmp(solver, 'gurobi') == 1
    %% Model parameters in GUROBI
    model.A = sparse([sparse(1,n) 1 ones(1,N)/(N*alpha)]);
    model.lb = [zeros(1,n) -inf zeros(1,N)];
    model.ub = inf(1,n+1+N);
    model.rhs = 0;
    % quad terms
    counter = 1;
    for j = 1:m
        S_j_squared = S(:,:,j).^2;
        S_j_squared = round(S_j_squared, 6);
        for i = 1:N
            model.quadcon(counter).Qc = sparse(1:(n + N + 1), 1:(n + N + 1), [S_j_squared(i, :), zeros(1, N + 1)]);
            model.quadcon(counter).q  = sparse(n + N + 1, 1);
            model.quadcon(counter).q(n + 1) = -1; model.quadcon(counter).q(n + 1 + i) = -1;
            model.quadcon(counter).rhs = U;
            counter = counter + 1;
        end
    end
    model.obj = [-ones(1, n) zeros(1,N+1)];
    model.modelsense = 'Min';
    model.sense = '<';

    %% Solve the problem using GUROBI
    params.outputflag = 0;
    params.timelimit = maxitime;
    params.optimalitytol = 1e-6;
    result = gurobi(model, params);
    if ~strcmp(result.status(1:7),'OPTIMAL')&&~strcmp(result.status(1:7),'SUBOPTI')
        fprintf('GUROBI_status: %s \n', result.status);
    else
        x = result.x(1:n);
        time = result.runtime;
    end

else
    %% define variables
    x = sdpvar(n,1); z = sdpvar(N,1); t = sdpvar(1,1);

    %% Define constraints
    Constraints = [t + sum(z)/(N*alpha) <= 0, z >= zeros(N,1), x >= zeros(n,1)];
    for j = 1:m
        Constraints = [Constraints, (S(:,:,j).^2)*(x.^2) - U*ones(N,1) - t*ones(N,1) <= z];
    end

    %% Define an objective
    Objective = -sum(x);

    %% solve by GUROBI solver
    ops = sdpsettings('solver', solver, 'verbose', 0, 'usex0', 0, 'gurobi.timelimit', maxitime);
    sol = optimize(Constraints, Objective, ops);
    x = value(x);
    time = sol.solvertime;
end
end
