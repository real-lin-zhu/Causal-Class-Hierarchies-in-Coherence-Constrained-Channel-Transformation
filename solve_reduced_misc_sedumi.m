function solution = solve_reduced_misc_sedumi( ...
    Jresource, Hclass, use_two_sided_diamond, sedumi_pars, mode)
%SOLVE_REDUCED_MISC_SEDUMI Reduced MISC primal or explicit dual via SeDuMi.
%
%   solution = solve_reduced_misc_sedumi(Jresource, Hclass, two_sided,
%       sedumi_pars, mode)
%
% Jresource is the 16-by-16 Choi matrix of the two resource channels in
% [Ai Ao Bi Bo] order. Hclass is 'Seq', 'ICO', or 'Par'. mode is 'primal'
% or 'dual'. The target is the qubit identity channel. Complex random Choi
% matrices are handled with SeDuMi complex-Hermitian PSD cones. mode may
% also be 'both', which returns solution.primal and solution.dual.

    if nargin < 5 || isempty(mode)
        mode = 'primal';
    end
    Hclass = validatestring(Hclass, {'Par', 'Seq', 'ICO'});
    mode = validatestring(mode, {'primal', 'dual', 'both'});
    if exist('sedumi', 'file') ~= 2
        error('SeDuMi was not found on the MATLAB path.');
    end

    Jresource = herm(full(Jresource));
    if ~isequal(size(Jresource), [16, 16])
        error('Jresource must be 16-by-16 for two qubit resource channels.');
    end

    model = build_reduced_misc_model( ...
        Jresource, Hclass, logical(use_two_sided_diamond));
    if strcmp(mode, 'primal')
        solution = solve_primal_model(model, sedumi_pars);
    elseif strcmp(mode, 'dual')
        solution = solve_explicit_dual_model(model, sedumi_pars);
    else
        solution = struct();
        solution.primal = solve_primal_model(model, sedumi_pars);
        solution.dual = solve_explicit_dual_model(model, sedumi_pars);
    end
    solution.Hclass = Hclass;
    solution.use_two_sided_diamond = logical(use_two_sided_diamond);
end


function model = build_reduced_misc_model(Jresource, Hclass, two_sided)

    DR = 16;
    nineq = 4 + 3 * double(two_sided);
    nlin = 4 + nineq;

    K = struct();
    K.l = nlin;
    K.s = [DR, DR, DR];
    K.scomplex = 1:3;
    nvar = nlin + 3 * DR^2;

    M1 = sparse(DR^2, nvar);
    Mp = sparse(DR^2, nvar);
    Mm = sparse(DR^2, nvar);
    offset = nlin;
    M1(:, offset + (1:DR^2)) = speye(DR^2);
    offset = offset + DR^2;
    Mp(:, offset + (1:DR^2)) = speye(DR^2);
    offset = offset + DR^2;
    Mm(:, offset + (1:DR^2)) = speye(DR^2);

    M0 = 0.5 * (Mp + Mm);
    Mc = 0.5 * (Mp - Mm);
    MS = M0 + M1;

    % trace(Jresource.'*T) = vec(Jresource).'*vec(T).
    jrow = reshape(Jresource, 1, []);
    ax0 = jrow * M0;
    ax1 = jrow * M1;
    axc = jrow * Mc;

    A = sparse(0, nvar);
    b = zeros(0, 1);
    complex_rows = false(0, 1);

    [Lcausal, bcausal, causal_complex] = ...
        reduced_causal_linear_system(Hclass);
    A = [A; Lcausal * MS]; %#ok<AGROW>
    b = [b; bcausal]; %#ok<AGROW>
    complex_rows = [complex_rows; causal_complex]; %#ok<AGROW>

    diag_selector = sparse(1:DR, 1:DR+1:DR^2, 1, DR, DR^2);
    A = [A; diag_selector * Mc]; %#ok<AGROW>
    b = [b; zeros(DR, 1)]; %#ok<AGROW>
    complex_rows = [complex_rows; false(DR, 1)]; %#ok<AGROW>

    iq1 = 1; iqp = 2; iqm = 3; iell = 4;
    slack = 5;
    [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
        A, b, complex_rows, sparse(1, iell, 1, 1, nvar) ...
        - 0.5*sparse(1, iqp, 1, 1, nvar) ...
        - 0.5*sparse(1, iqm, 1, 1, nvar) ...
        - sparse(1, iq1, 1, 1, nvar), 0, slack);
    [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
        A, b, complex_rows, sparse(1, iq1, 1, 1, nvar) - ax1, 0, slack);
    [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
        A, b, complex_rows, sparse(1, iqm, 1, 1, nvar) - ax0 + axc, 0, slack);
    [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
        A, b, complex_rows, sparse(1, iqp, 1, 1, nvar) - ax0 - axc, -2, slack);

    if two_sided
        [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
            A, b, complex_rows, sparse(1, iq1, 1, 1, nvar) + ax1, 0, slack);
        [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
            A, b, complex_rows, sparse(1, iqm, 1, 1, nvar) + ax0 - axc, 0, slack);
        [A, b, complex_rows, slack] = add_nonnegative_inequality( ...
            A, b, complex_rows, sparse(1, iqp, 1, 1, nvar) + ax0 + axc, 2, slack);
    end
    assert(slack == nlin + 1, 'Internal SeDuMi slack count mismatch.');

    [A, b] = expand_complex_equalities(A, b, complex_rows);
    A = hermitianize_cone_coefficients(A, K);
    [A, b] = independent_real_equalities(A, b);

    c = sparse(nvar, 1);
    c(iell) = 1;
    if two_sided
        c(iell) = 0.5;
    end

    model = struct();
    model.A = A;
    model.b = b;
    model.c = c;
    model.K = K;
    model.M0 = M0;
    model.M1 = M1;
    model.Mc = Mc;
    model.Jresource = Jresource;
    model.iq1 = iq1;
    model.iqp = iqp;
    model.iqm = iqm;
    model.iell = iell;
    model.two_sided = two_sided;
    model.Hclass = Hclass;
end


function solution = solve_primal_model(model, pars)

    [x, y, info] = sedumi(model.A, model.b, model.c, model.K, pars);
    status = sedumi_status(info);
    optval = real(full(model.c.' * x));

    DR = 16;
    T0 = herm(reshape(full(model.M0 * x), DR, DR));
    T1 = herm(reshape(full(model.M1 * x), DR, DR));
    Tc = herm(reshape(full(model.Mc * x), DR, DR));
    S = T0 + T1;
    Jt = model.Jresource.';
    x0 = real(trace(Jt * T0));
    x1 = real(trace(Jt * T1));
    xc = real(trace(Jt * Tc));
    z0 = real(0.5 * (x(model.iqp) + x(model.iqm)));
    z1 = real(x(model.iq1));
    zc = real(0.5 * (x(model.iqp) - x(model.iqm)));
    ell = real(x(model.iell));

    O0 = omega0(); O1 = omega1(); Oc = omegac();
    Jout = x0 * O0 + x1 * O1 + xc * Oc;
    Z = z0 * O0 + z1 * O1 + zc * Oc;
    Jidentity = O0 + Oc;
    [lhs, rhs] = reduced_causal_constraint_cells(S, solution_class(model));

    solution = struct();
    solution.status = status;
    solution.distance = optval;
    solution.sedumi_optval = optval;
    solution.sedumi_info = info;
    solution.sedumi_dual = y;
    solution.T0 = T0;
    solution.T1 = T1;
    solution.Tc = Tc;
    solution.S = S;
    solution.z0 = z0;
    solution.z1 = z1;
    solution.zc = zc;
    solution.ell = ell;
    solution.x0 = x0;
    solution.x1 = x1;
    solution.xc = xc;
    solution.Jout = Jout;
    solution.Z = Z;
    solution.J_identity = Jidentity;
    solution.min_eig_T1 = min_eig_herm(T1);
    solution.min_eig_T0_plus_Tc = min_eig_herm(T0 + Tc);
    solution.min_eig_T0_minus_Tc = min_eig_herm(T0 - Tc);
    solution.min_eig_Z = min_eig_herm(Z);
    solution.min_eig_Z_minus_diff = min_eig_herm(Z - (Jout - Jidentity));
    solution.min_eig_Z_plus_diff = min_eig_herm(Z + (Jout - Jidentity));
    solution.ell_trace_margin = ell - (z0 + z1);
    solution.causal_residual_fro = cell_residual_fro(lhs, rhs);
    solution.misc_diag_residual = norm(diag(Tc), 2);
end


function solution = solve_explicit_dual_model(model, pars)

    % Dual of min c'*x, A*x=b, x in K:
    % max b'*y subject to c-A'*y in K. Introduce the cone slack explicitly
    % and solve it as a second, independent SeDuMi problem.
    A = model.A; b = model.b; c = model.c; K = model.K;
    m = size(A, 1);
    n = size(A, 2);

    Kd = struct();
    Kd.f = m;
    Kd.l = K.l;
    Kd.s = K.s;
    Kd.scomplex = K.scomplex;

    Bfull = [A', speye(n)];
    [selector, complex_equalities] = cone_upper_selector(K);
    Bd = selector * Bfull;
    bd = selector * c;
    Kd.ycomplex = find(complex_equalities).';
    cd = [-b; zeros(n, 1)];

    [xd, yd, info] = sedumi(Bd, bd, cd, Kd, pars);
    status = sedumi_status(info);
    y = real(full(xd(1:m)));
    cone_slack = xd(m + (1:n));
    lower_bound = real(full(b.' * y));

    solution = struct();
    solution.status = status;
    solution.lower_bound = lower_bound;
    solution.sedumi_optval = lower_bound;
    solution.sedumi_info = info;
    solution.equality_multiplier = y;
    solution.cone_slack = cone_slack;
    solution.dual_of_explicit_model = yd;
    solution.dual_equality_residual = norm( ...
        selector * (A' * y + cone_slack - c));
    solution.min_cone_eigenvalue = min_cone_eigenvalue(cone_slack, K);
end


function Hclass = solution_class(model)
    Hclass = model.Hclass;
end


function [A, b, complex_rows, next_slack] = add_nonnegative_inequality( ...
    A, b, complex_rows, expression_row, rhs, slack_index)

    expression_row(slack_index) = expression_row(slack_index) - 1;
    A = [A; expression_row]; %#ok<AGROW>
    b = [b; rhs]; %#ok<AGROW>
    complex_rows = [complex_rows; false]; %#ok<AGROW>
    next_slack = slack_index + 1;
end


function [Areal, breal] = expand_complex_equalities(A, b, complex_rows)

    nout = numel(b) + sum(complex_rows);
    Areal = sparse(nout, size(A, 2));
    breal = zeros(nout, 1);
    cursor = 0;
    for row = 1:numel(b)
        cursor = cursor + 1;
        Areal(cursor, :) = A(row, :);
        breal(cursor) = real(b(row));
        if complex_rows(row)
            cursor = cursor + 1;
            Areal(cursor, :) = -1i * A(row, :);
            breal(cursor) = imag(b(row));
        end
    end
end


function A = hermitianize_cone_coefficients(A, K)

    offset = K.l;
    for block = 1:numel(K.s)
        d = K.s(block);
        cols = offset + (1:d^2);
        for row = 1:size(A, 1)
            R = reshape(full(A(row, cols)), d, d);
            C = R.';
            C = 0.5 * (C + C');
            A(row, cols) = reshape(C.', 1, []);
        end
        offset = offset + d^2;
    end
    A = sparse(A);
end


function [Aind, bind] = independent_real_equalities(A, b)

    effective = [real(A), imag(A)];
    [~, R, pivot] = qr(full(effective.'), 'vector');
    diagonal = abs(diag(R));
    if isempty(diagonal)
        numerical_rank = 0;
    else
        tolerance = max(size(effective)) * eps(max(diagonal));
        numerical_rank = sum(diagonal > tolerance);
    end
    keep = sort(pivot(1:numerical_rank));
    Aind = sparse(A(keep, :));
    bind = b(keep);
end


function [selector, complex_rows] = cone_upper_selector(K)

    n = K.l + sum(K.s.^2);
    row_indices = (1:K.l).';
    complex_rows = false(K.l, 1);
    offset = K.l;
    for block = 1:numel(K.s)
        d = K.s(block);
        [r, c] = find(triu(ones(d)));
        local = sub2ind([d, d], r, c);
        row_indices = [row_indices; offset + local]; %#ok<AGROW>
        complex_rows = [complex_rows; r ~= c]; %#ok<AGROW>
        offset = offset + d^2;
    end
    selector = sparse(1:numel(row_indices), row_indices, 1, ...
        numel(row_indices), n);
end


function value = min_cone_eigenvalue(x, K)

    value = min(real(x(1:K.l)));
    offset = K.l;
    for block = 1:numel(K.s)
        d = K.s(block);
        X = herm(reshape(x(offset + (1:d^2)), d, d));
        value = min(value, min(real(eig(X))));
        offset = offset + d^2;
    end
end


function status = sedumi_status(info)

    if isfield(info, 'pinf') && info.pinf ~= 0
        status = 'Primal infeasible';
    elseif isfield(info, 'dinf') && info.dinf ~= 0
        status = 'Dual infeasible';
    elseif ~isfield(info, 'numerr') || info.numerr == 0
        status = 'Solved';
    elseif info.numerr == 1
        status = 'Inaccurate/Solved';
    else
        status = 'Failed';
    end
end


function [L, b, complex_rows] = reduced_causal_linear_system(Hclass)

    persistent cache
    key = validatestring(Hclass, {'Par', 'Seq', 'ICO'});
    if ~isempty(cache) && isfield(cache, key)
        entry = cache.(key);
        L = entry.L; b = entry.b; complex_rows = entry.complex_rows;
        return;
    end

    DR = 16;
    [base, complex_rows] = causal_residual_vector(zeros(DR), key);
    L = zeros(numel(base), DR^2);
    for j = 1:DR^2
        Ej = zeros(DR);
        Ej(j) = 1;
        L(:, j) = causal_residual_vector(Ej, key) - base;
    end
    b = -base;
    entry = struct('L', sparse(L), 'b', b, 'complex_rows', complex_rows);
    if isempty(cache), cache = struct(); end
    cache.(key) = entry;
    L = entry.L;
end


function [residual, complex_rows] = causal_residual_vector(S, Hclass)

    [lhs, rhs] = reduced_causal_constraint_cells(S, Hclass);
    residual = zeros(0, 1);
    complex_rows = false(0, 1);
    for k = 1:numel(lhs)
        difference = full(lhs{k} - rhs{k});
        if isscalar(difference)
            residual = [residual; difference]; %#ok<AGROW>
            complex_rows = [complex_rows; false]; %#ok<AGROW>
        else
            mask = triu(true(size(difference, 1)));
            [row, col] = find(mask);
            residual = [residual; difference(mask)]; %#ok<AGROW>
            complex_rows = [complex_rows; row ~= col]; %#ok<AGROW>
        end
    end
end


function [lhs, rhs] = reduced_causal_constraint_cells(S, Hclass)

    dimsR = [2, 2, 2, 2];
    lhs = {}; rhs = {};
    if strcmp(Hclass, 'Par')
        lhs{end+1} = S;
        core = ptrace_kron(S, dimsR, [2, 4]);
        rhs{end+1} = insert_identity_kron(core, dimsR, [1, 3], [2, 4], 1/4);
        lhs{end+1} = ptrace_kron(S, dimsR, [1, 3]);
        rhs{end+1} = eye(4);
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;
    elseif strcmp(Hclass, 'Seq')
        lhs{end+1} = S;
        core = ptrace_kron(S, dimsR, 4);
        rhs{end+1} = insert_identity_kron(core, dimsR, [1, 2, 3], 4, 1/2);
        lhs{end+1} = ptrace_kron(S, dimsR, [3, 4]);
        core = ptrace_kron(S, dimsR, [2, 3, 4]);
        rhs{end+1} = insert_identity_kron(core, [2, 2], 1, 2, 1/2);
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;
    else
        lhs{end+1} = S;
        coreA = ptrace_kron(S, dimsR, 2);
        termA = insert_identity_kron(coreA, dimsR, [1, 3, 4], 2, 1/2);
        coreB = ptrace_kron(S, dimsR, 4);
        termB = insert_identity_kron(coreB, dimsR, [1, 2, 3], 4, 1/2);
        coreAB = ptrace_kron(S, dimsR, [2, 4]);
        termAB = insert_identity_kron(coreAB, dimsR, [1, 3], [2, 4], 1/4);
        rhs{end+1} = termA + termB - termAB;
        lhs{end+1} = ptrace_kron(S, dimsR, [1, 2]);
        core = ptrace_kron(S, dimsR, [1, 2, 4]);
        rhs{end+1} = insert_identity_kron(core, [2, 2], 1, 2, 1/2);
        lhs{end+1} = ptrace_kron(S, dimsR, [3, 4]);
        core = ptrace_kron(S, dimsR, [2, 3, 4]);
        rhs{end+1} = insert_identity_kron(core, [2, 2], 1, 2, 1/2);
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;
    end
end


function O = omega0()
    O = diag([1, 0, 0, 1]);
end

function O = omega1()
    O = diag([0, 1, 1, 0]);
end

function O = omegac()
    O = zeros(4); O(1, 4) = 1; O(4, 1) = 1;
end


function Y = ptrace_kron(X, dims, trace_sys)
    selectors = partial_trace_selectors_kron(dims, trace_sys);
    Y = X(1, 1) * zeros(size(selectors{1}, 1));
    for k = 1:numel(selectors)
        E = selectors{k}; Y = Y + E * X * E';
    end
end


function selectors = partial_trace_selectors_kron(dims, trace_sys)
    dims = double(dims(:).'); n = numel(dims);
    trace_sys = sort(double(trace_sys(:).'));
    keep_sys = setdiff(1:n, trace_sys, 'stable');
    Dfull = prod(dims); Dkeep = prod(dims(keep_sys));
    if isempty(trace_sys), selectors = {speye(Dfull)}; return; end
    trace_subs = all_subscripts_kron(dims(trace_sys));
    keep_subs = all_subscripts_kron(dims(keep_sys));
    selectors = cell(size(trace_subs, 1), 1);
    for t = 1:size(trace_subs, 1)
        cols = zeros(Dkeep, 1);
        for r = 1:Dkeep
            sub = zeros(1, n); sub(keep_sys) = keep_subs(r, :);
            sub(trace_sys) = trace_subs(t, :);
            cols(r) = subscript_kron_to_ind(dims, sub);
        end
        selectors{t} = sparse(1:Dkeep, cols, 1, Dkeep, Dfull);
    end
end


function Y = insert_identity_kron(T, full_dims, keep_sys, id_sys, scale)
    full_dims = double(full_dims(:).');
    keep_sys = double(keep_sys(:).'); id_sys = double(id_sys(:).');
    temp_order = [keep_sys, id_sys]; dims_temp = full_dims(temp_order);
    Ttemp = kron(T, scale * speye(prod(full_dims(id_sys))));
    perm = zeros(1, numel(full_dims));
    for s = 1:numel(full_dims), perm(s) = find(temp_order == s); end
    P = permutation_unitary_kron(dims_temp, perm);
    Y = P * Ttemp * P';
end


function P = permutation_unitary_kron(dims, permutation)
    dims = double(dims(:).'); permutation = double(permutation(:).');
    D = prod(dims); row_idx = zeros(D, 1);
    for col = 1:D
        sub = ind_to_subscript_kron(dims, col);
        row_idx(col) = subscript_kron_to_ind(dims(permutation), sub(permutation));
    end
    P = sparse(row_idx, (1:D).', 1, D, D);
end


function subs = all_subscripts_kron(dims)
    dims = double(dims(:).'); D = prod(dims);
    if isempty(dims), subs = zeros(1, 0); return; end
    subs = zeros(D, numel(dims));
    for idx = 1:D, subs(idx, :) = ind_to_subscript_kron(dims, idx); end
end


function sub = ind_to_subscript_kron(dims, idx)
    dims = double(dims(:).'); sub = zeros(1, numel(dims)); idx = idx - 1;
    for k = numel(dims):-1:1
        sub(k) = mod(idx, dims(k)) + 1; idx = floor(idx / dims(k));
    end
end


function idx = subscript_kron_to_ind(dims, sub)
    idx = 0;
    for k = 1:numel(dims), idx = idx * dims(k) + sub(k) - 1; end
    idx = idx + 1;
end


function X = herm(X)
    X = 0.5 * (X + X');
end


function value = min_eig_herm(X)
    value = min(real(eig(herm(full(X)))));
end


function residuals = cell_residual_fro(lhs, rhs)
    residuals = zeros(1, numel(lhs));
    for k = 1:numel(lhs), residuals(k) = norm(lhs{k} - rhs{k}, 'fro'); end
end
