function [x, time, iter] = DCA(S, U, alpha, opts)

fprintf('***************** DCA ***************** \n');

%% parameter setting
[N, n, m] = size(S);
M = ceil((1-alpha)*N); T = N - M; maxiter = 2*1e2;

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
    x0 = randn(n,m);
end
if isfield(opts, 'tol')
    tol = opts.tol;
else
    tol = 1e-6;
end

xk = x0; fval = -sum(x0); time = 0;

for iter = 1 : maxiter
    fval_old = fval;
    %% compute subgradient of C(x,\xi^l) for l = 1,...,N at xk
    for j = 1:m
        C(:,j) = (S(:,:,j).^2)*(xk.^2) - U;
    end
    [CN, inx_set]= max(C,[],2);
    %% compute subgradient of H(x) at xk
    [mz, I] = maxk(CN, N-M);
    mz_fval = sum(mz);
    gk = zeros(n,1);
    for i = 1:size(I,1)
        inx = inx_set(I(i));
        sk = 2*(S(I(i),:,inx).^2)'.*xk;
        gk = gk + sk;
    end

    if strcmp(solver, 'gurobi') == 1
        %% Model parameters in GUROBI
        % x z lam mu
        model.A = sparse([-gk' sparse(1,N) ones(1,N) T+1;...
            sparse(N,n) speye(N) -speye(N) -ones(N,1)]);
        model.rhs = [mz_fval-gk'*xk zeros(1,N)];
        model.lb = [zeros(1,n) -inf(1,N) zeros(1,N) -inf];
        model.ub = inf(1,n+2*N+1);
        % quad terms
        counter = 1;
        for j = 1:m
            S_j_squared = S(:,:,j).^2;
            S_j_squared = round(S_j_squared, 6);
            for i = 1:N
                model.quadcon(counter).Qc = sparse(1:(n+2*N+1), 1:(n+2*N+1), [S_j_squared(i,:), zeros(1, 2*N+1)]);
                model.quadcon(counter).q  = sparse(n+2*N+1, 1); model.quadcon(counter).q(n+i) = -1;
                model.quadcon(counter).rhs = U;
                counter = counter + 1;
            end
        end
        model.obj = [-ones(1, n) zeros(1,2*N+1)];
        model.modelsense = 'Min';
        model.sense = '<';

        %% Solve the problem using GUROBI
        params.outputflag = 0;
        params.optimalitytol = tol;
        result = gurobi(model, params);
        time = time + result.runtime;
        xk = result.x(1:n);
    else
        %% define variables
        x = sdpvar(n,1); z = sdpvar(N,1); lambda = sdpvar(N,1); muv = sdpvar(1);

        %% intial point
        assign(x, x0);

        %% Define constraints
        Constraints = [
            sum(lambda) + (T+1)*muv - mz_fval - trace(gk'*(x-xk)) <= 0,...
            z - lambda - muv*ones(N,1) <= zeros(N,1), lambda >= zeros(N,1), x >= zeros(n,1)];
        for j = 1:m
            Constraints = [Constraints, (S(:,:,j).^2)*(x.^2) - U <= z];
        end

        %% Define an objective
        Objective = -sum(x);

        %% solve by GUROBI solver
        ops = sdpsettings('solver', solver, 'verbose', 0, 'usex0', 1, 'gurobi.timelimit', maxitime);
        sol = optimize(Constraints, Objective, ops);
        time = time + sol.solvertime;

        %% report iterate information
        xk = value(x);
    end

    fval = -sum(xk);
    fprintf('Iternum: %d, fval: %.4f\n', iter, fval);

    if abs(fval - fval_old)/max(1,fval_old) <= tol || time > 1800
        break;
    end

end
x = xk;
end