%% reduced_misc_dual_sdp.m
% Reduced dual SDP corresponding to reduced_misc_primal_sdp.m.
%
% This implements the symmetry-reduced dual described in Appendix C of
% Causal-Class Hierarchies in Coherence-Constrained Channel Transformation.
% The paper writes the common dual after the substitutions
%
%   U = T0 + Tc,  V = T0 - Tc,  W = T1,
%   S = (U + V)/2 + W.
%
% The dual variables are:
%   r1, rplus, rminus: scalar Watrous multipliers,
%   YH: reduced comb-constraint multipliers,
%   Hc: reduced MISC multiplier.
%
% The causal adjoints below use the same [Ai Ao Bi Bo] tensor convention as
% reduced_misc_primal_sdp.m and experiment_AD_to_ID.m.

clear; clc;

%% Settings

error_rates = (0:19) / 20;
classes_to_solve = {'Par', 'Seq', 'ICO'};

% For Par and ICO, check the three dual PSD slacks only on the
% slot-exchange-invariant 10+6 block subspace, matching the slot reduction.
use_slot_blocks = true;

% SeDuMi accuracy and display settings. Set fid = 1 for solver output.
sedumi_pars = struct('eps', 1e-9, 'bigeps', 1e-7, 'fid', 0);

save_outputs = true;
output_file = fullfile('results', 'data', 'AD_to_identity_reduced_misc_dual.mat');

if exist('sedumi', 'file') ~= 2
    error(['SeDuMi was not found on the MATLAB path. Install SeDuMi and ', ...
        'run its setup script before running this file.']);
end

if save_outputs
    output_dir = fileparts(output_file);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
end

fprintf('\n====================================================\n');
fprintf('Symmetry-reduced MISC dual SDP\n');
fprintf('AD noise epsilon grid = [%s]\n', sprintf(' %.2f', error_rates));
fprintf('Slot block reduction for Par/ICO = %d\n', use_slot_blocks);
fprintf('====================================================\n');

sweep = repmat(struct('error_rate', [], 'solutions', []), numel(error_rates), 1);

for e = 1:numel(error_rates)
    error_rate = error_rates(e);
    JN = choi_amplitude_damping(error_rate);
    GN = herm(kron(JN, JN).');  % scalar link functional ell_N(T)=<GN,T>

    fprintf('\n----------------------------------------------------\n');
    fprintf('AD noise epsilon = %.2f  (%d/%d)\n', ...
        error_rate, e, numel(error_rates));
    fprintf('----------------------------------------------------\n');

    solutions = struct();
    for k = 1:numel(classes_to_solve)
        Hclass = classes_to_solve{k};
        fprintf('\nSolving reduced %s dual SDP...\n', Hclass);

        solutions.(Hclass) = solve_reduced_misc_dual( ...
            GN, Hclass, use_slot_blocks, sedumi_pars);

        sol = solutions.(Hclass);
        fprintf('  status       : %s\n', sol.status);
        fprintf('  lower bound  : %.12g\n', sol.lower_bound);
        fprintf('  r1,r+,r-     : %.12g  %.12g  %.12g\n', ...
            sol.r1, sol.rplus, sol.rminus);
        fprintf('  min eig slack: U %.3e, V %.3e, W %.3e\n', ...
            sol.min_eig_slack_U, sol.min_eig_slack_V, sol.min_eig_slack_W);
    end

    sweep(e).error_rate = error_rate;
    sweep(e).solutions = solutions;
end

if save_outputs
    solutions_by_error = sweep; %#ok<NASGU>
    save(output_file, ...
        'error_rates', 'use_slot_blocks', 'classes_to_solve', ...
        'sedumi_pars', 'sweep', 'solutions_by_error');
    fprintf('\nSaved reduced-dual results to %s\n', output_file);
end


%% ========================================================================
%%                              DUAL SOLVER
%% ========================================================================

function solution = solve_reduced_misc_dual( ...
    GN, Hclass, use_slot_blocks, sedumi_pars)

    Hclass = validatestring(Hclass, {'Par', 'Seq', 'ICO'});

    dimsR = [2, 2, 2, 2];  % [Ai Ao Bi Bo]
    DR = prod(dimsR);

    slot_P = slot_exchange_unitary();
    use_blocks_here = use_slot_blocks && any(strcmp(Hclass, {'Par', 'ICO'}));

    if use_blocks_here
        [Qplus, Qminus] = slot_exchange_bases(slot_P);
    else
        Qplus = [];
        Qminus = [];
    end

    % Free variables contain the diagonal MISC multiplier and symmetric
    % causal multipliers. Nonnegative variables are r1,r+,r- and their
    % three upper-bound slacks. The three PSD slacks are explicit SeDuMi
    % cone blocks tied to the affine dual expressions by equalities.
    [MA, objective_row, meta] = dual_affine_maps(Hclass);
    nfree = size(MA, 2);

    if use_blocks_here
        projection = {kron(Qplus.', Qplus'), kron(Qminus.', Qminus')};
        one_family_dims = [size(Qplus, 2), size(Qminus, 2)];
    else
        projection = {speye(DR^2)};
        one_family_dims = DR;
    end

    K = struct();
    K.f = nfree;
    K.l = 6;
    K.s = repmat(one_family_dims, 1, 3);
    nvar = nfree + K.l + sum(K.s.^2);

    ir1 = nfree + 1;
    irplus = nfree + 2;
    irminus = nfree + 3;
    is1 = nfree + 4;
    isplus = nfree + 5;
    isminus = nfree + 6;

    MAfull = sparse(DR^2, nvar);
    MAfull(:, 1:nfree) = MA;
    MD = sparse(DR^2, nvar);
    MD(:, meta.hc) = diagonal_matrix_map(DR);
    g = sparse(GN(:));

    MU = -0.5 * MAfull - 0.5 * MD;
    MU(:, irplus) = g;
    MV = -0.5 * MAfull + 0.5 * MD;
    MV(:, irminus) = g;
    MW = -MAfull;
    MW(:, ir1) = g;
    slack_maps = {MU, MV, MW};

    A = sparse(3, nvar);
    A(1, [ir1, is1]) = 1;
    A(2, [irplus, isplus]) = 1;
    A(3, [irminus, isminus]) = 1;
    b = [1; 0.5; 0.5];

    offset = nfree + K.l;
    cone_block = 0;
    for family = 1:3
        for parity = 1:numel(projection)
            cone_block = cone_block + 1;
            dblock = K.s(cone_block);
            cols = offset + (1:dblock^2);
            equality_map = sparse(dblock^2, nvar);
            equality_map(:, cols) = speye(dblock^2);
            equality_map = equality_map - projection{parity} * slack_maps{family};
            upper = upper_triangle_selector(dblock);
            A = [A; upper * equality_map]; %#ok<AGROW>
            b = [b; zeros(size(upper, 1), 1)]; %#ok<AGROW>
            offset = offset + dblock^2;
        end
    end

    A = symmetrize_semidefinite_coefficients(A, K);
    [A, b] = independent_equalities(A, b);
    c = sparse(nvar, 1);
    c(1:nfree) = -objective_row.';
    c(irplus) = 2;

    [x, y, info] = sedumi(A, b, c, K, sedumi_pars);
    status = sedumi_status(info);
    primal_minimum = full(c.' * x);
    lower_bound = -primal_minimum;

    Hc_value = diag(full(x(meta.hc)));
    Y0_value = herm(reshape(meta.Y0_map * x(meta.y0), DR, DR));
    Y1_value = herm(reshape(meta.Y1_map * x(meta.y1), 4, 4));
    if isempty(meta.y2)
        Y2_value = [];
    else
        Y2_value = herm(reshape(meta.Y2_map * x(meta.y2), 4, 4));
    end
    if isempty(meta.eta)
        eta_value = [];
    else
        eta_value = full(x(meta.eta));
    end

    Astar_value = herm(reshape(MA * x(1:nfree), DR, DR));
    DAB_Hc_value = Hc_value;
    r1_value = full(x(ir1));
    rplus_value = full(x(irplus));
    rminus_value = full(x(irminus));
    slack_U_value = rplus_value * GN - 0.5 * Astar_value - 0.5 * DAB_Hc_value;
    slack_V_value = rminus_value * GN - 0.5 * Astar_value + 0.5 * DAB_Hc_value;
    slack_W_value = r1_value * GN - Astar_value;

    solution = struct();
    solution.Hclass = Hclass;
    solution.status = status;
    solution.lower_bound = lower_bound;
    solution.cvx_optval = lower_bound; % Retained for old result consumers.
    solution.sedumi_optval = lower_bound;
    solution.sedumi_info = info;
    solution.sedumi_dual = y;
    solution.use_slot_blocks = use_blocks_here;

    solution.r1 = r1_value;
    solution.rplus = rplus_value;
    solution.rminus = rminus_value;
    solution.Hc = Hc_value;
    solution.Astar = Astar_value;
    if strcmp(Hclass, 'Par')
        solution.dual_vars = struct('Y0', Y0_value, 'Y1', Y1_value);
    elseif strcmp(Hclass, 'Seq')
        solution.dual_vars = struct( ...
            'Y0', Y0_value, 'Y1', Y1_value, 'eta', eta_value);
    elseif strcmp(Hclass, 'ICO')
        solution.dual_vars = struct( ...
            'Y0', Y0_value, 'Y1', Y1_value, 'Y2', Y2_value, ...
            'eta', eta_value);
    end

    solution.slack_U = slack_U_value;
    solution.slack_V = slack_V_value;
    solution.slack_W = slack_W_value;

    [solution.min_eig_slack_U, solution.min_eig_slack_U_blocks] = ...
        min_eig_for_slack(slack_U_value, Qplus, Qminus, use_blocks_here);
    [solution.min_eig_slack_V, solution.min_eig_slack_V_blocks] = ...
        min_eig_for_slack(slack_V_value, Qplus, Qminus, use_blocks_here);
    [solution.min_eig_slack_W, solution.min_eig_slack_W_blocks] = ...
        min_eig_for_slack(slack_W_value, Qplus, Qminus, use_blocks_here);

    solution.r_bounds_residual = [
        solution.r1, 1 - solution.r1, ...
        solution.rplus, 0.5 - solution.rplus, ...
        solution.rminus, 0.5 - solution.rminus];
end


%% ========================================================================
%%                         REDUCED CAUSAL ADJOINTS
%% ========================================================================

function [MA, objective_row, meta] = dual_affine_maps(Hclass)

    persistent cache
    key = validatestring(Hclass, {'Par', 'Seq', 'ICO'});
    if ~isempty(cache) && isfield(cache, key)
        entry = cache.(key);
        MA = entry.MA;
        objective_row = entry.objective_row;
        meta = entry.meta;
        return;
    end

    DR = 16;
    Y0_map = symmetric_matrix_map(DR);
    Ysmall_map = symmetric_matrix_map(4);

    cursor = 0;
    meta = struct();
    meta.hc = cursor + (1:DR);
    cursor = cursor + DR;
    meta.y0 = cursor + (1:size(Y0_map, 2));
    cursor = cursor + size(Y0_map, 2);
    meta.y1 = cursor + (1:size(Ysmall_map, 2));
    cursor = cursor + size(Ysmall_map, 2);
    if strcmp(key, 'ICO')
        meta.y2 = cursor + (1:size(Ysmall_map, 2));
        cursor = cursor + size(Ysmall_map, 2);
    else
        meta.y2 = [];
    end
    if any(strcmp(key, {'Seq', 'ICO'}))
        meta.eta = cursor + 1;
        cursor = cursor + 1;
    else
        meta.eta = [];
    end
    nfree = cursor;

    MA = zeros(DR^2, nfree);
    zero4 = zeros(4);
    for j = 1:numel(meta.y0)
        Y0 = reshape(Y0_map(:, j), DR, DR);
        MA(:, meta.y0(j)) = dual_astar_value(Y0, zero4, zero4, 0, key);
    end
    for j = 1:numel(meta.y1)
        Y1 = reshape(Ysmall_map(:, j), 4, 4);
        MA(:, meta.y1(j)) = dual_astar_value(zeros(DR), Y1, zero4, 0, key);
    end
    for j = 1:numel(meta.y2)
        Y2 = reshape(Ysmall_map(:, j), 4, 4);
        MA(:, meta.y2(j)) = dual_astar_value(zeros(DR), zero4, Y2, 0, key);
    end
    if ~isempty(meta.eta)
        MA(:, meta.eta) = dual_astar_value(zeros(DR), zero4, zero4, 1, key);
    end

    objective_row = zeros(1, nfree);
    if strcmp(key, 'Par')
        for j = 1:numel(meta.y1)
            objective_row(meta.y1(j)) = trace( ...
                reshape(Ysmall_map(:, j), 4, 4));
        end
    else
        objective_row(meta.eta) = 4;
    end

    meta.Y0_map = Y0_map;
    meta.Y1_map = Ysmall_map;
    meta.Y2_map = Ysmall_map;
    entry = struct('MA', sparse(MA), ...
        'objective_row', objective_row, 'meta', meta);
    if isempty(cache)
        cache = struct();
    end
    cache.(key) = entry;
    MA = entry.MA;
end


function vector = dual_astar_value(Y0, Y1, Y2, eta, Hclass)

    dimsR = [2, 2, 2, 2];
    DR = 16;
    if strcmp(Hclass, 'Par')
        Astar = Y0 - replace_subsystems_kron(Y0, dimsR, [2, 4]) ...
            + insert_identity_kron(Y1, dimsR, [2, 4], [1, 3], 1);
    elseif strcmp(Hclass, 'Seq')
        K1 = Y1 - replace_subsystems_kron(Y1, [2, 2], 2);
        Astar = Y0 - replace_subsystems_kron(Y0, dimsR, 4) ...
            + insert_identity_kron(K1, dimsR, [1, 2], [3, 4], 1) ...
            + eta * eye(DR);
    else
        K1 = Y1 - replace_subsystems_kron(Y1, [2, 2], 2);
        K2 = Y2 - replace_subsystems_kron(Y2, [2, 2], 2);
        Astar = Y0 ...
            - replace_subsystems_kron(Y0, dimsR, 2) ...
            - replace_subsystems_kron(Y0, dimsR, 4) ...
            + replace_subsystems_kron(Y0, dimsR, [2, 4]) ...
            + insert_identity_kron(K1, dimsR, [3, 4], [1, 2], 1) ...
            + insert_identity_kron(K2, dimsR, [1, 2], [3, 4], 1) ...
            + eta * eye(DR);
    end
    vector = Astar(:);
end


function map = symmetric_matrix_map(n)

    [row, col] = find(triu(ones(n)));
    map = sparse(n^2, numel(row));
    for j = 1:numel(row)
        map(sub2ind([n, n], row(j), col(j)), j) = 1;
        if row(j) ~= col(j)
            map(sub2ind([n, n], col(j), row(j)), j) = 1;
        end
    end
end


function map = diagonal_matrix_map(n)

    map = sparse(1:n+1:n^2, 1:n, 1, n^2, n);
end


function selector = upper_triangle_selector(n)

    [row, col] = find(triu(ones(n)));
    linear = sub2ind([n, n], row, col);
    selector = sparse(1:numel(linear), linear, 1, numel(linear), n^2);
end


function A = symmetrize_semidefinite_coefficients(A, K)

    offset = K.f + K.l;
    for block = 1:numel(K.s)
        d = K.s(block);
        cols = offset + (1:d^2);
        for row = 1:size(A, 1)
            coefficient = reshape(full(A(row, cols)), d, d);
            coefficient = 0.5 * (coefficient + coefficient.');
            A(row, cols) = coefficient(:).';
        end
        offset = offset + d^2;
    end
    A = sparse(A);
end


function [Aind, bind] = independent_equalities(A, b)

    [~, R, pivot] = qr(full(A.'), 'vector');
    diagonal = abs(diag(R));
    if isempty(diagonal)
        numerical_rank = 0;
    else
        tolerance = max(size(A)) * eps(max(diagonal));
        numerical_rank = sum(diagonal > tolerance);
    end
    keep = sort(pivot(1:numerical_rank));
    Aind = sparse(A(keep, :));
    bind = b(keep);
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


%% ========================================================================
%%                             CHANNELS
%% ========================================================================

function J = choi_amplitude_damping(eps_ad)

    A0 = [1, 0; 0, sqrt(1 - eps_ad)];
    A1 = [0, sqrt(eps_ad); 0, 0];
    J = choi_from_kraus({A0, A1});
end


function J = choi_from_kraus(Ks)

    % Unnormalized Choi matrix:
    %   J_N = sum_ij |i><j| tensor N(|i><j|).
    %
    % Tensor order is [input output], in kron-order basis.

    d_in = size(Ks{1}, 2);
    d_out = size(Ks{1}, 1);
    J = zeros(d_in * d_out, d_in * d_out);

    for i = 1:d_in
        for j = 1:d_in
            Eij = zeros(d_in, d_in);
            Eij(i, j) = 1;

            Nij = zeros(d_out, d_out);
            for a = 1:numel(Ks)
                K = Ks{a};
                Nij = Nij + K * Eij * K';
            end

            for o1 = 1:d_out
                for o2 = 1:d_out
                    row = subscript_kron_to_ind([d_in, d_out], [i, o1]);
                    col = subscript_kron_to_ind([d_in, d_out], [j, o2]);
                    J(row, col) = Nij(o1, o2);
                end
            end
        end
    end

    J = herm(J);
end


%% ========================================================================
%%                         SLOT-EXCHANGE REDUCTION
%% ========================================================================

function P = slot_exchange_unitary()

    % [Ai Ao Bi Bo] -> [Bi Bo Ai Ao]
    P = permutation_unitary_kron([2, 2, 2, 2], [3, 4, 1, 2]);
end


function [Qplus, Qminus] = slot_exchange_bases(P)

    P = full(P);
    Pplus = 0.5 * (eye(size(P)) + P);
    Pminus = 0.5 * (eye(size(P)) - P);

    Qplus = orth(Pplus);
    Qminus = orth(Pminus);

    if size(Qplus, 2) ~= 10 || size(Qminus, 2) ~= 6
        error('Unexpected slot-symmetry block dimensions: %d and %d.', ...
            size(Qplus, 2), size(Qminus, 2));
    end
end


%% ========================================================================
%%                      KRON-ORDER LINEAR ALGEBRA HELPERS
%% ========================================================================

function Y = replace_subsystems_kron(X, dims, repl)

    repl = double(repl(:).');
    keep = setdiff(1:numel(dims), repl, 'stable');
    drepl = prod(dims(repl));

    Xkeep = ptrace_kron(X, dims, repl);
    Y = insert_identity_kron(Xkeep, dims, keep, repl, 1 / drepl);
end


function Y = ptrace_kron(X, dims, trace_sys)

    selectors = partial_trace_selectors_kron(dims, trace_sys);
    Dkeep = size(selectors{1}, 1);
    Y = X(1, 1) * zeros(Dkeep, Dkeep);

    for k = 1:numel(selectors)
        E = selectors{k};
        Y = Y + E * X * E';
    end
end


function selectors = partial_trace_selectors_kron(dims, trace_sys)

    dims = double(dims(:).');
    n = numel(dims);
    trace_sys = sort(double(trace_sys(:).'));
    keep_sys = setdiff(1:n, trace_sys, 'stable');

    Dfull = prod(dims);
    Dkeep = prod(dims(keep_sys));

    if isempty(trace_sys)
        selectors = {speye(Dfull)};
        return;
    end

    trace_subs = all_subscripts_kron(dims(trace_sys));
    keep_subs = all_subscripts_kron(dims(keep_sys));

    selectors = cell(size(trace_subs, 1), 1);
    for t = 1:size(trace_subs, 1)
        rows = zeros(Dkeep, 1);
        cols = zeros(Dkeep, 1);

        for r = 1:Dkeep
            full_sub = zeros(1, n);
            full_sub(keep_sys) = keep_subs(r, :);
            full_sub(trace_sys) = trace_subs(t, :);

            rows(r) = r;
            cols(r) = subscript_kron_to_ind(dims, full_sub);
        end

        selectors{t} = sparse(rows, cols, ones(Dkeep, 1), Dkeep, Dfull);
    end
end


function Y = insert_identity_kron(T, full_dims, keep_sys, id_sys, scale)

    full_dims = double(full_dims(:).');
    keep_sys = double(keep_sys(:).');
    id_sys = double(id_sys(:).');

    temp_order = [keep_sys, id_sys];
    dims_temp = full_dims(temp_order);

    Did = prod(full_dims(id_sys));
    Ttemp = kron(T, scale * speye(Did));

    perm = zeros(1, numel(full_dims));
    for s = 1:numel(full_dims)
        perm(s) = find(temp_order == s);
    end

    P = permutation_unitary_kron(dims_temp, perm);
    Y = P * Ttemp * P';
end


function P = permutation_unitary_kron(dims, permutation)

    dims = double(dims(:).');
    permutation = double(permutation(:).');

    if numel(permutation) ~= numel(dims) || ...
            ~isequal(sort(permutation), 1:numel(dims))
        error('permutation must be a permutation of 1:length(dims).');
    end

    D = prod(dims);
    dims_out = dims(permutation);
    row_idx = zeros(D, 1);
    col_idx = (1:D).';

    for col = 1:D
        sub_in = ind_to_subscript_kron(dims, col);
        sub_out = sub_in(permutation);
        row_idx(col) = subscript_kron_to_ind(dims_out, sub_out);
    end

    P = sparse(row_idx, col_idx, ones(D, 1), D, D);
end


function subs = all_subscripts_kron(dims)

    dims = double(dims(:).');
    D = prod(dims);

    if isempty(dims)
        subs = zeros(1, 0);
        return;
    end

    subs = zeros(D, numel(dims));
    for idx = 1:D
        subs(idx, :) = ind_to_subscript_kron(dims, idx);
    end
end


function sub = ind_to_subscript_kron(dims, idx)

    dims = double(dims(:).');
    sub = zeros(1, numel(dims));
    idx0 = idx - 1;

    for k = numel(dims):-1:1
        sub(k) = mod(idx0, dims(k)) + 1;
        idx0 = floor(idx0 / dims(k));
    end
end


function idx = subscript_kron_to_ind(dims, sub)

    dims = double(dims(:).');
    sub = double(sub(:).');

    idx0 = 0;
    for k = 1:numel(dims)
        idx0 = idx0 * dims(k) + (sub(k) - 1);
    end

    idx = idx0 + 1;
end


function X = herm(X)

    X = 0.5 * (X + X');
end


function [lambda_min, block_values] = ...
    min_eig_for_slack(X, Qplus, Qminus, use_blocks)

    X = herm(full(X));

    if use_blocks
        block_values = [
            min(real(eig(herm(Qplus' * X * Qplus)))), ...
            min(real(eig(herm(Qminus' * X * Qminus))))];
        lambda_min = min(block_values);
    else
        block_values = min(real(eig(X)));
        lambda_min = block_values;
    end
end
