function [x, time] = MIP(S, U, upper, alpha, opts)

fprintf('***************** MIP ***************** \n')

%% parameter setting
[N, n, m] = size(S);

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
    model.A = sparse([sparse(1,n) ones(1,N)]);
    model.lb = [zeros(1,n) zeros(1,N)];
    model.ub = [inf(1,n) ones(1,N)];
    model.rhs = alpha*N;
    model.sense = '<';
    % quad terms
    counter = 1;
    for j = 1:m
        S_j_squared = S(:,:,j).^2;
        S_j_squared = round(S_j_squared, 6);
        for i = 1:N
            model.quadcon(counter).Qc = sparse(1:(n + N), 1:(n + N), [S_j_squared(i,:), zeros(1, N)]);
            model.quadcon(counter).q  = sparse(n + N, 1); model.quadcon(counter).q(n+i) = -upper;
            model.quadcon(counter).rhs = U;
            counter = counter + 1;
        end
    end

    model.obj = [-ones(1, n) zeros(1,N)];
    model.modelsense = 'Min';
    model.vtype = [repmat('C',n,1); repmat('B',N,1)];

    %% Solve the problem using GUROBI
    params.outputflag = 0;
    params.optimalitytol = 1e-6;
    params.timelimit = maxitime;
    result = gurobi(model, params);
    x = result.x(1:n);
    time = result.runtime;
else
    %% define variables
    x = sdpvar(n,1); y = binvar(N,1);

    %% Define constraints
    Constraints = [sum(y) <= alpha*N, x >= zeros(n,1)];
    for j = 1:m
        Constraints = [Constraints, (S(:,:,j).^2)*(x.^2) - U <= upper*y];
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
