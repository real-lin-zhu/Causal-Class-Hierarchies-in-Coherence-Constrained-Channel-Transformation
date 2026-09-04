%% reduced_misc_primal_sdp.m
% Output-symmetry-reduced primal SDP from Appendix C of
% Causal-Class Hierarchies in Coherence-Constrained Channel Transformation.
%
% Reduced variables:
%   Theta = T0 kron Omega0 + T1 kron Omega1 + Tc kron Omegac,
%   Z     = z0 Omega0 + z1 Omega1 + zc Omegac.
%
% Tensor order for T0, T1, Tc:
%   [Ai Ao Bi Bo], all qubits.
%
% This script is self-contained apart from SeDuMi. It uses kron-order indexing
% (the last tensor factor is fastest), matching MATLAB kron.
%
% The Appendix subsection writes the output-symmetry reduction. The causal
% equations below are the same reduction in the tensor convention used by
% experiment_AD_to_ID.m / ChannelSimulation_2_to_1.m, where the comb
% identities are inserted on the output registers Ao and Bo.

clear; clc;

%% Settings

% Use an integer index so the intended 20-point grid is unambiguous:
% 0, 0.05, 0.10, ..., 0.95

error_rates = (0:19) / 20;

classes_to_solve = {'Par', 'Seq', 'ICO'};

% Exact slot-exchange block reduction for Par and ICO:
%   T = Qplus*Tplus*Qplus' + Qminus*Tminus*Qminus'
% with block sizes 10 and 6. This is not applied to the fixed-order Seq SDP.
use_slot_blocks = true;

% Exact diamond-distance SDP for a Hermitian channel difference:
%   Z >= Diff, Z >= -Diff, Z >= 0, minimize ||Tr_output Z||_inf / 2.
% Set this false only to reproduce the one-sided relaxation.
use_two_sided_diamond = false;

% SeDuMi accuracy and display settings.  Set fid = 1 for solver output.
sedumi_pars = struct('eps', 1e-9, 'bigeps', 1e-7, 'fid', 0);

save_outputs = true;
output_file = fullfile('results', 'data', 'AD_to_identity_reduced_misc_primal.mat');
save_figure = true;
figure_output_file = fullfile('results', 'figures', 'AD_to_identity_reduced_misc_primal_sweep.png');

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
fprintf('Symmetry-reduced MISC primal SDP\n');
fprintf('AD noise epsilon grid = [%s]\n', sprintf(' %.2f', error_rates));
fprintf('Slot block reduction for Par/ICO = %d\n', use_slot_blocks);
fprintf('Two-sided diamond constraints = %d\n', use_two_sided_diamond);
fprintf('====================================================\n');

sweep = repmat(struct('error_rate', [], 'solutions', []), numel(error_rates), 1);

for e = 1:numel(error_rates)
    error_rate = error_rates(e);
    JN = choi_amplitude_damping(error_rate);
    JN2 = kron(JN, JN);  % [Ai Ao Bi Bo]

    fprintf('\n----------------------------------------------------\n');
    fprintf('AD noise epsilon = %.2f  (%d/%d)\n', ...
        error_rate, e, numel(error_rates));
    fprintf('----------------------------------------------------\n');

    solutions = struct();
    for k = 1:numel(classes_to_solve)
        Hclass = classes_to_solve{k};
        fprintf('\nSolving reduced %s SDP...\n', Hclass);

        solutions.(Hclass) = solve_reduced_misc_primal( ...
            JN2, Hclass, use_slot_blocks, use_two_sided_diamond, ...
            sedumi_pars);

        sol = solutions.(Hclass);
        fprintf('  status       : %s\n', sol.status);
        fprintf('  distance     : %.12g\n', sol.distance);
        fprintf('  ell          : %.12g\n', sol.ell);
        fprintf('  x0, x1, xc   : %.12g  %.12g  %.12g\n', sol.x0, sol.x1, sol.xc);
        fprintf('  min eig tests: T1 %.3e, T0+Tc %.3e, T0-Tc %.3e\n', ...
            sol.min_eig_T1, sol.min_eig_T0_plus_Tc, sol.min_eig_T0_minus_Tc);
        fprintf('  max residual : causal %.3e, diag(Tc) %.3e, slot %.3e\n', ...
            max(sol.causal_residual_fro), sol.misc_diag_residual, sol.slot_residual_fro);
    end

    sweep(e).error_rate = error_rate;
    sweep(e).solutions = solutions;
end

if save_outputs
    solutions_by_error = sweep; %#ok<NASGU>
    save(output_file, ...
        'error_rates', 'use_slot_blocks', 'use_two_sided_diamond', ...
        'classes_to_solve', 'sedumi_pars', 'sweep', 'solutions_by_error');
    fprintf('\nSaved reduced-MISC SDP sweep results to %s\n', output_file);
end

distances = sweep_distances(sweep, classes_to_solve);
analytical_grid = linspace(min(error_rates), max(error_rates), 501);
analytical_ico = analytical_ico_bound(analytical_grid);
if save_figure
    figure_dir = fileparts(figure_output_file);
    if ~isempty(figure_dir) && ~exist(figure_dir, 'dir')
        mkdir(figure_dir);
    end

    fig = figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    plot(error_rates, distances(:, 1), '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
    hold on;
    plot(error_rates, distances(:, 2), '-s', 'LineWidth', 1.8, 'MarkerSize', 5);
    plot(error_rates, distances(:, 3), '-^', 'LineWidth', 1.8, 'MarkerSize', 5);
    plot(analytical_grid, analytical_ico, 'k--', 'LineWidth', 2.0);
    hold off;
    grid on;
    xlabel('AD error rate \epsilon');
    ylabel('Diamond distance');
    title('Reduced MISC simulation of identity from two AD channels');
    legend([classes_to_solve, {'Analytical general-process bound'}], ...
        'Location', 'northwest');
    xlim([min(error_rates), max(error_rates)]);
    xticks(error_rates);
    xticklabels(compose('%.2f', error_rates));
    xtickangle(45);
    ylim([0, 0.5]);

    exportgraphics(fig, figure_output_file, 'Resolution', 200);
    fprintf('Saved reduced-MISC sweep figure to %s\n', figure_output_file);
end


%% ========================================================================
%%                              SDP SOLVER
%% ========================================================================

function solution = solve_reduced_misc_primal( ...
    JN2, Hclass, use_slot_blocks, use_two_sided_diamond, sedumi_pars)

    Hclass = validatestring(Hclass, {'Par', 'Seq', 'ICO'});

    dimsR = [2, 2, 2, 2];  % [Ai Ao Bi Bo]
    DR = prod(dimsR);
    JN2_t = JN2.';

    slot_P = slot_exchange_unitary();
    use_blocks_here = use_slot_blocks && any(strcmp(Hclass, {'Par', 'ICO'}));

    if use_blocks_here
        [Qplus, Qminus] = slot_exchange_bases(slot_P);
        nplus = size(Qplus, 2);
        nminus = size(Qminus, 2);
    else
        Qplus = [];
        Qminus = [];
        nplus = 0;
        nminus = 0;
    end

    % Use PSD variables P1=T1, Pp=T0+Tc, and Pm=T0-Tc.  This makes
    % every matrix inequality a native SeDuMi semidefinite cone.
    if use_blocks_here
        block_dims = repmat([nplus, nminus], 1, 3);
        embed = {kron(Qplus, Qplus), kron(Qminus, Qminus)};
    else
        block_dims = [DR, DR, DR];
        embed = {speye(DR^2)};
    end

    % Linear cone entries: q1=z1, qp=z0+zc, qm=z0-zc, ell, then one
    % nonnegative slack for each remaining scalar inequality.
    nineq = 4 + 3 * double(use_two_sided_diamond);
    nlin = 4 + nineq;
    K = struct();
    K.l = nlin;
    K.s = block_dims;
    nvar = nlin + sum(block_dims.^2);

    M1 = sparse(DR^2, nvar);
    Mp = sparse(DR^2, nvar);
    Mm = sparse(DR^2, nvar);
    offset = nlin;
    for family = 1:3
        for parity = 1:numel(embed)
            dblock = block_dims((family-1)*numel(embed) + parity);
            cols = offset + (1:dblock^2);
            if family == 1
                M1(:, cols) = embed{parity};
            elseif family == 2
                Mp(:, cols) = embed{parity};
            else
                Mm(:, cols) = embed{parity};
            end
            offset = offset + dblock^2;
        end
    end

    M0 = 0.5 * (Mp + Mm);
    Mc = 0.5 * (Mp - Mm);
    MS = M0 + M1;
    jrow = reshape(JN2_t.', 1, []);
    ax0 = jrow * M0;
    ax1 = jrow * M1;
    axc = jrow * Mc;

    A = sparse(0, nvar);
    b = zeros(0, 1);

    [Lcausal, bcausal] = reduced_causal_linear_system(Hclass);
    A = [A; Lcausal * MS]; %#ok<AGROW>
    b = [b; bcausal]; %#ok<AGROW>

    diag_selector = sparse(1:DR, 1:DR+1:DR^2, 1, DR, DR^2);
    A = [A; diag_selector * Mc]; %#ok<AGROW>
    b = [b; zeros(DR, 1)]; %#ok<AGROW>

    if ~use_blocks_here && any(strcmp(Hclass, {'Par', 'ICO'}))
        symop = speye(DR^2) - kron(slot_P, slot_P);
        upper = upper_triangle_selector(DR);
        A = [A; upper * symop * M0; upper * symop * M1; ...
            upper * symop * Mc]; %#ok<AGROW>
        b = [b; zeros(3 * size(upper, 1), 1)]; %#ok<AGROW>
    end

    iq1 = 1; iqp = 2; iqm = 3; iell = 4;
    slack = 5;
    [A, b, slack] = add_nonnegative_inequality( ...
        A, b, sparse(1, iell, 1, 1, nvar) ...
        - 0.5*sparse(1, iqp, 1, 1, nvar) ...
        - 0.5*sparse(1, iqm, 1, 1, nvar) ...
        - sparse(1, iq1, 1, 1, nvar), 0, slack);
    [A, b, slack] = add_nonnegative_inequality( ...
        A, b, sparse(1, iq1, 1, 1, nvar) - ax1, 0, slack);
    [A, b, slack] = add_nonnegative_inequality( ...
        A, b, sparse(1, iqm, 1, 1, nvar) - ax0 + axc, 0, slack);
    [A, b, slack] = add_nonnegative_inequality( ...
        A, b, sparse(1, iqp, 1, 1, nvar) - ax0 - axc, -2, slack);

    if use_two_sided_diamond
        [A, b, slack] = add_nonnegative_inequality( ...
            A, b, sparse(1, iq1, 1, 1, nvar) + ax1, 0, slack);
        [A, b, slack] = add_nonnegative_inequality( ...
            A, b, sparse(1, iqm, 1, 1, nvar) + ax0 - axc, 0, slack);
        [A, b, slack] = add_nonnegative_inequality( ...
            A, b, sparse(1, iqp, 1, 1, nvar) + ax0 + axc, 2, slack);
    end
    assert(slack == nlin + 1, 'Internal SeDuMi slack count mismatch.');

    A = symmetrize_semidefinite_coefficients(A, K);
    [A, b] = independent_equalities(A, b);
    c = sparse(nvar, 1);
    c(iell) = 1;
    if use_two_sided_diamond
        c(iell) = 0.5;
    end

    [x, y, info] = sedumi(A, b, c, K, sedumi_pars);
    status = sedumi_status(info);
    optval = full(c.' * x);

    T0_value = herm(reshape(full(M0 * x), DR, DR));
    T1_value = herm(reshape(full(M1 * x), DR, DR));
    Tc_value = herm(reshape(full(Mc * x), DR, DR));
    S_value = T0_value + T1_value;

    x0_value = real(trace(JN2_t * T0_value));
    x1_value = real(trace(JN2_t * T1_value));
    xc_value = real(trace(JN2_t * Tc_value));

    [causal_lhs_value, causal_rhs_value] = ...
        reduced_causal_constraint_cells(S_value, Hclass);

    solution = struct();
    solution.Hclass = Hclass;
    solution.status = status;
    solution.cvx_optval = optval; % Retained for compatibility with old files.
    solution.sedumi_optval = optval;
    solution.sedumi_info = info;
    solution.sedumi_dual = y;
    solution.ell = full(x(iell));
    solution.distance = optval;
    solution.use_two_sided_diamond = use_two_sided_diamond;
    solution.use_slot_blocks = use_blocks_here;
    solution.slot_block_dims = [nplus, nminus];

    solution.T0 = T0_value;
    solution.T1 = T1_value;
    solution.Tc = Tc_value;
    solution.S = S_value;
    solution.z0 = full(0.5 * (x(iqp) + x(iqm)));
    solution.z1 = full(x(iq1));
    solution.zc = full(0.5 * (x(iqp) - x(iqm)));
    solution.x0 = x0_value;
    solution.x1 = x1_value;
    solution.xc = xc_value;

    solution.Omega0 = omega0();
    solution.Omega1 = omega1();
    solution.Omegac = omegac();
    solution.Jout = x0_value * solution.Omega0 ...
        + x1_value * solution.Omega1 ...
        + xc_value * solution.Omegac;
    solution.Z = solution.z0 * solution.Omega0 ...
        + solution.z1 * solution.Omega1 ...
        + solution.zc * solution.Omegac;
    solution.J_identity = solution.Omega0 + solution.Omegac;

    solution.min_eig_T1 = min_eig_herm(T1_value);
    solution.min_eig_T0_plus_Tc = min_eig_herm(T0_value + Tc_value);
    solution.min_eig_T0_minus_Tc = min_eig_herm(T0_value - Tc_value);
    solution.min_eig_Z = min_eig_herm(solution.Z);
    solution.min_eig_Z_minus_diff = min_eig_herm( ...
        solution.Z - (solution.Jout - solution.J_identity));
    solution.min_eig_Z_plus_diff = min_eig_herm( ...
        solution.Z + (solution.Jout - solution.J_identity));
    solution.ell_trace_margin = solution.ell - (solution.z0 + solution.z1);

    solution.causal_residual_fro = cell_residual_fro( ...
        causal_lhs_value, causal_rhs_value);
    solution.misc_diag_residual = norm(diag(Tc_value), 2);

    if any(strcmp(Hclass, {'Par', 'ICO'}))
        solution.slot_residual_fro = max([ ...
            norm(T0_value - slot_P*T0_value*slot_P', 'fro'), ...
            norm(T1_value - slot_P*T1_value*slot_P', 'fro'), ...
            norm(Tc_value - slot_P*Tc_value*slot_P', 'fro')]);
    else
        solution.slot_residual_fro = 0;
    end
end


function [A, b, next_slack] = add_nonnegative_inequality( ...
    A, b, expression_row, rhs, slack_index)

    % expression_row*x >= rhs is represented by
    % expression_row*x - nonnegative_slack = rhs.
    expression_row(slack_index) = expression_row(slack_index) - 1;
    A = [A; expression_row]; %#ok<AGROW>
    b = [b; rhs]; %#ok<AGROW>
    next_slack = slack_index + 1;
end


function [Aind, bind] = independent_equalities(A, b)

    % SeDuMi expects a full-row-rank equality system.  The causal matrix
    % equations deliberately contain some redundant scalar equations.
    if isempty(A)
        Aind = A;
        bind = b;
        return;
    end

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


function A = symmetrize_semidefinite_coefficients(A, K)

    % Only the symmetric part of an SDP-block coefficient contributes to
    % its inner product with a symmetric cone variable.  Make that explicit
    % before the equality-rank test so hidden antisymmetric dependencies are
    % removed as well.
    offset = K.l;
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


function selector = upper_triangle_selector(n)

    [row, col] = find(triu(ones(n)));
    linear = sub2ind([n, n], row, col);
    selector = sparse(1:numel(linear), linear, 1, numel(linear), n^2);
end


function [L, b] = reduced_causal_linear_system(Hclass)

    % Return L and b such that L*vec(S)=b.  Building this operator once per
    % class keeps the tensor-indexing implementation in one tested place.
    persistent cache
    key = validatestring(Hclass, {'Par', 'Seq', 'ICO'});
    if ~isempty(cache) && isfield(cache, key)
        L = cache.(key).L;
        b = cache.(key).b;
        return;
    end

    DR = 16;
    zeroS = zeros(DR);
    base = causal_residual_vector(zeroS, key);
    L = zeros(numel(base), DR^2);
    for j = 1:DR^2
        Ej = zeros(DR);
        Ej(j) = 1;
        L(:, j) = causal_residual_vector(Ej, key) - base;
    end
    b = -base;
    entry = struct('L', sparse(L), 'b', b);
    if isempty(cache)
        cache = struct();
    end
    cache.(key) = entry;
    L = entry.L;
end


function residual = causal_residual_vector(S, Hclass)

    [lhs, rhs] = reduced_causal_constraint_cells(S, Hclass);
    residual = zeros(0, 1);
    for k = 1:numel(lhs)
        difference = full(lhs{k} - rhs{k});
        if isscalar(difference)
            residual = [residual; difference]; %#ok<AGROW>
        else
            mask = triu(true(size(difference, 1)));
            residual = [residual; difference(mask)]; %#ok<AGROW>
        end
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


function distances = sweep_distances(sweep, classes_to_solve)

    distances = NaN(numel(sweep), numel(classes_to_solve));

    for e = 1:numel(sweep)
        for k = 1:numel(classes_to_solve)
            Hclass = classes_to_solve{k};
            distances(e, k) = sweep(e).solutions.(Hclass).distance;
        end
    end
end


function U = analytical_ico_bound(eps_ad)

    eps_ad = min(max(eps_ad, 0), 1);
    s = sqrt(1 - eps_ad);
    eps_c = (-9 + sqrt(657)) / 32;

    U = zeros(size(eps_ad));
    first_branch = eps_ad <= eps_c;
    second_branch = ~first_branch;

    U(first_branch) = 0.5 * (1 - s(first_branch) .* ...
        (eps_ad(first_branch) + ...
        sqrt(eps_ad(first_branch).^2 + 1 - eps_ad(first_branch))));

    U(second_branch) = 0.5 * (1 - 2 * s(second_branch) .* ...
        (4 * eps_ad(second_branch) / 5 + 2 * s(second_branch) / 5));
end


%% ========================================================================
%%                         REDUCED CAUSAL CONSTRAINTS
%% ========================================================================

function [lhs, rhs] = reduced_causal_constraint_cells(S, Hclass)

    dimsR = [2, 2, 2, 2];  % [Ai Ao Bi Bo]
    lhs = {};
    rhs = {};

    if strcmp(Hclass, 'Par')
        % Experiment convention:
        % S = Tr_AoBo S tensor I_AoBo/4.
        lhs{end+1} = S;
        core = ptrace_kron(S, dimsR, [2, 4]);  % keep [Ai Bi]
        rhs{end+1} = insert_identity_kron(core, dimsR, [1, 3], [2, 4], 1/4);

        % Tr_AiBi S = I_AoBo.
        lhs{end+1} = ptrace_kron(S, dimsR, [1, 3]);  % keep [Ao Bo]
        rhs{end+1} = eye(4);

        % Tr S = 4.
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;

    elseif strcmp(Hclass, 'Seq')
        % Fixed order A before B.
        % S = Tr_Bo S tensor I_Bo/2.
        lhs{end+1} = S;
        core = ptrace_kron(S, dimsR, 4);  % keep [Ai Ao Bi]
        rhs{end+1} = insert_identity_kron(core, dimsR, [1, 2, 3], 4, 1/2);

        % Tr_BiBo S = Tr_AoBiBo S tensor I_Ao/2.
        dims_AiAo = [2, 2];
        lhs{end+1} = ptrace_kron(S, dimsR, [3, 4]);      % keep [Ai Ao]
        core = ptrace_kron(S, dimsR, [2, 3, 4]);         % keep [Ai]
        rhs{end+1} = insert_identity_kron(core, dims_AiAo, 1, 2, 1/2);

        % Tr S = 4.
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;

    elseif strcmp(Hclass, 'ICO')
        % Experiment convention:
        % S = Tr_Ao S tensor I_Ao/2
        %   + Tr_Bo S tensor I_Bo/2
        %   - Tr_AoBo S tensor I_AoBo/4.
        lhs{end+1} = S;
        coreA = ptrace_kron(S, dimsR, 2);       % keep [Ai Bi Bo]
        termA = insert_identity_kron(coreA, dimsR, [1, 3, 4], 2, 1/2);
        coreB = ptrace_kron(S, dimsR, 4);       % keep [Ai Ao Bi]
        termB = insert_identity_kron(coreB, dimsR, [1, 2, 3], 4, 1/2);
        coreAB = ptrace_kron(S, dimsR, [2, 4]); % keep [Ai Bi]
        termAB = insert_identity_kron(coreAB, dimsR, [1, 3], [2, 4], 1/4);
        rhs{end+1} = termA + termB - termAB;

        % Tr_AiAo S = Tr_AiAoBo S tensor I_Bo/2.
        dims_BiBo = [2, 2];
        lhs{end+1} = ptrace_kron(S, dimsR, [1, 2]);  % keep [Bi Bo]
        core = ptrace_kron(S, dimsR, [1, 2, 4]);     % keep [Bi]
        rhs{end+1} = insert_identity_kron(core, dims_BiBo, 1, 2, 1/2);

        % Tr_BiBo S = Tr_AoBiBo S tensor I_Ao/2.
        dims_AiAo = [2, 2];
        lhs{end+1} = ptrace_kron(S, dimsR, [3, 4]);  % keep [Ai Ao]
        core = ptrace_kron(S, dimsR, [2, 3, 4]);     % keep [Ai]
        rhs{end+1} = insert_identity_kron(core, dims_AiAo, 1, 2, 1/2);

        % Tr S = 4.
        lhs{end+1} = trace(S);
        rhs{end+1} = 4;

    else
        error('Unknown Hclass. Use Par, Seq, or ICO.');
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
%%                       OUTPUT-SYMMETRY MATRICES
%% ========================================================================

function O = omega0()

    O = zeros(4, 4);
    O(1, 1) = 1;  % |00><00|
    O(4, 4) = 1;  % |11><11|
end


function O = omega1()

    O = zeros(4, 4);
    O(2, 2) = 1;  % |01><01|
    O(3, 3) = 1;  % |10><10|
end


function O = omegac()

    O = zeros(4, 4);
    O(1, 4) = 1;  % |00><11|
    O(4, 1) = 1;  % |11><00|
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

    % Make the dimensions explicit; this catches basis-order mistakes.
    if size(Qplus, 2) ~= 10 || size(Qminus, 2) ~= 6
        error('Unexpected slot-symmetry block dimensions: %d and %d.', ...
            size(Qplus, 2), size(Qminus, 2));
    end
end


%% ========================================================================
%%                      KRON-ORDER LINEAR ALGEBRA HELPERS
%% ========================================================================

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


function lambda_min = min_eig_herm(X)

    X = herm(full(X));
    lambda_min = min(real(eig(X)));
end


function residuals = cell_residual_fro(lhs, rhs)

    residuals = zeros(1, numel(lhs));

    for k = 1:numel(lhs)
        residuals(k) = norm(full(lhs{k} - rhs{k}), 'fro');
    end
end
