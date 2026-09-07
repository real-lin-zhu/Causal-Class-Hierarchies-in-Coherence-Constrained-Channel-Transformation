function results = general_nonreduced_channel_distance_sdp(JN1, varargin)
%% general_nonreduced_channel_distance_sdp.m
% Full non-reduced SDP distance calculator for two resource channels
% transformed into one target channel.
%
% Usage:
%   results = general_nonreduced_channel_distance_sdp(JN, JT)
%       Uses two copies of the same input channel: JN tensor JN -> JT.
%
%   results = general_nonreduced_channel_distance_sdp(JN1, JN2, JT)
%       Uses two different input channels: JN1 tensor JN2 -> JT.
%
% Optional name-value pairs:
%   'Strategies'                    {'Par','Seq','ICO'}
%   'FreeClasses'                   {'MISC','DISC'}
%   'Dims'                          [dAi dAo dBi dBo dCi dCo]
%   'UseTwoSidedDiamond'            true
%   'EnforceSlotSymmetry'           []
%   'EnforceSequentialSlotSymmetry' false
%   'CvxSolver'                     ''
%   'CvxPrecision'                  'best'
%   'Verbose'                       true
%
% Tensor order for the superchannel Choi variable is
%   [Ai Ao Bi Bo Ci Co].
%
% Choi convention is unnormalized:
%   J_N = sum_ij |i><j| tensor N(|i><j|).
%
% If 'Dims' is omitted, the function assumes all three channels are square
% channels and infers d from sqrt(size(J,1)).
%
% EnforceSlotSymmetry defaults to [] = auto. In auto mode, slot-swap
% symmetry is imposed only when JN1 and JN2 are numerically equal. For
% different input channels, forcing slot symmetry is an extra restriction
% and can incorrectly make ICO worse than Seq.

    if nargin < 2
        error('Provide at least JN1 and JT, or JN1, JN2, and JT.');
    end

    [JN2, JT, opt_args] = split_channel_arguments(JN1, varargin);
    opts = parse_options(opt_args{:});

    if isempty(opts.Dims)
        dims = infer_square_channel_dims(JN1, JN2, JT);
    else
        dims = double(opts.Dims(:).');
        validate_dims(dims, JN1, JN2, JT);
    end

    JN1 = herm_numeric(full(JN1));
    JN2 = herm_numeric(full(JN2));
    JT = herm_numeric(full(JT));

    [opts.EnforceSlotSymmetry, opts.SlotSymmetryAuto] = ...
        resolve_slot_symmetry(opts.EnforceSlotSymmetry, JN1, JN2);

    Jresource = kron(JN1, JN2);

    strategies = normalize_string_list(opts.Strategies);
    free_classes = normalize_string_list(opts.FreeClasses);

    distances = NaN(numel(strategies), numel(free_classes));
    statuses = cell(numel(strategies), numel(free_classes));
    solutions = struct();

    if opts.Verbose
        fprintf('\n====================================================\n');
        fprintf('General full non-reduced MISC/DISC channel-distance SDP\n');
        fprintf('dims [Ai Ao Bi Bo Ci Co] = [%s]\n', sprintf(' %d', dims));
        fprintf('strategies  = %s\n', strjoin(strategies, ', '));
        fprintf('free classes = %s\n', strjoin(free_classes, ', '));
        fprintf('slot symmetry Par/ICO = %d', opts.EnforceSlotSymmetry);
        if opts.SlotSymmetryAuto
            fprintf(' (auto)');
        end
        fprintf(', Seq = %d\n', opts.EnforceSequentialSlotSymmetry);
        fprintf('two-sided diamond constraints = %d\n', opts.UseTwoSidedDiamond);
        fprintf('====================================================\n');
    end

    for r = 1:numel(free_classes)
        free_class = free_classes{r};
        solutions.(free_class) = struct();

        for s = 1:numel(strategies)
            strategy = strategies{s};
            enforce_here = slot_symmetry_setting( ...
                strategy, opts.EnforceSlotSymmetry, ...
                opts.EnforceSequentialSlotSymmetry);

            if opts.Verbose
                fprintf('\nSolving full %s %s SDP...\n', free_class, strategy);
            end

            sol = solve_nonreduced_superchannel_primal_with_retries( ...
                Jresource, JT, dims, strategy, free_class, enforce_here, ...
                opts.UseTwoSidedDiamond, opts.CvxSolver, opts.CvxPrecision, ...
                opts.Verbose);

            solutions.(free_class).(strategy) = sol;
            distances(s, r) = sol.distance;
            statuses{s, r} = sol.status;

            if opts.Verbose
                fprintf('  status       : %s\n', sol.status);
                fprintf('  distance     : %.12g\n', sol.distance);
                fprintf('  min eig tests: Theta %.3e, Z %.3e, Z-Diff %.3e, Z+Diff %.3e\n', ...
                    sol.min_eig_Theta, sol.min_eig_Z, ...
                    sol.min_eig_Z_minus_Diff, sol.min_eig_Z_plus_Diff);
                fprintf('  max residual : causal %.3e, %s %.3e, slot %.3e\n', ...
                    sol.max_causal_residual, free_class, ...
                    sol.free_residual, sol.max_slot_residual);
                fprintf('  free residuals: MISC %.3e, DISC %.3e\n', ...
                    sol.misc_residual, sol.disc_residual);
            end
        end
    end

    raw_distances = distances;
    [distances, hierarchy_adjustments] = apply_hierarchy_upper_bounds( ...
        distances, solutions, strategies, free_classes, opts.Verbose);

    results = struct();
    results.JN1 = JN1;
    results.JN2 = JN2;
    results.JT = JT;
    results.Jresource = Jresource;
    results.dims = dims;
    results.strategies = strategies;
    results.free_classes = free_classes;
    results.distances = distances;
    results.raw_distances = raw_distances;
    results.statuses = statuses;
    results.solutions = solutions;
    results.hierarchy_adjustments = hierarchy_adjustments;
    results.options = opts;

    if opts.Verbose
        print_distance_table(results);
    end
end


%% ========================================================================
%%                              ARGUMENTS
%% ========================================================================

function [JN2, JT, opt_args] = split_channel_arguments(JN1, args)

    if isempty(args)
        error('Missing target Choi matrix JT.');
    end

    first = args{1};
    if ~isnumeric(first)
        error('Second argument must be JT or JN2.');
    end

    if numel(args) >= 2 && isnumeric(args{2})
        JN2 = first;
        JT = args{2};
        opt_args = args(3:end);
    else
        JN2 = JN1;
        JT = first;
        opt_args = args(2:end);
    end
end


function opts = parse_options(varargin)

    opts = struct();
    opts.Strategies = {'Par', 'Seq', 'ICO'};
    opts.FreeClasses = {'MISC', 'DISC'};
    opts.Dims = [];
    opts.UseTwoSidedDiamond = true;
    opts.EnforceSlotSymmetry = [];
    opts.EnforceSequentialSlotSymmetry = false;
    opts.SlotSymmetryAuto = true;
    opts.CvxSolver = '';
    opts.CvxPrecision = 'best';
    opts.Verbose = true;

    if mod(numel(varargin), 2) ~= 0
        error('Optional arguments must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = char(varargin{k});
        value = varargin{k + 1};
        key = lower(strrep(name, '_', ''));

        switch key
            case 'strategies'
                opts.Strategies = value;
            case 'freeclasses'
                opts.FreeClasses = value;
            case 'dims'
                opts.Dims = value;
            case 'usetwosideddiamond'
                opts.UseTwoSidedDiamond = logical(value);
            case 'enforceslotsymmetry'
                if isempty(value)
                    opts.EnforceSlotSymmetry = [];
                else
                    opts.EnforceSlotSymmetry = logical(value);
                end
            case 'enforcesequentialslotsymmetry'
                opts.EnforceSequentialSlotSymmetry = logical(value);
            case 'cvxsolver'
                opts.CvxSolver = char(value);
            case 'cvxprecision'
                opts.CvxPrecision = char(value);
            case 'verbose'
                opts.Verbose = logical(value);
            otherwise
                error('Unknown option: %s.', name);
        end
    end
end


function values = normalize_string_list(values)

    if ischar(values) || isstring(values)
        values = cellstr(values);
    end

    if ~iscell(values)
        error('Expected a string or cell array of strings.');
    end

    values = cellfun(@char, values, 'UniformOutput', false);
end


function [enforce_slot_symmetry, was_auto] = resolve_slot_symmetry( ...
    requested, JN1, JN2)

    was_auto = isempty(requested);
    if ~was_auto
        enforce_slot_symmetry = logical(requested);
        if enforce_slot_symmetry && ~channels_equal(JN1, JN2)
            warning(['Slot symmetry is being forced even though JN1 and JN2 ', ...
                'are different. This adds a restriction and can make ICO ', ...
                'appear worse than Seq.']);
        end
        return;
    end

    enforce_slot_symmetry = channels_equal(JN1, JN2);
end


function tf = channels_equal(JN1, JN2)

    if ~isequal(size(JN1), size(JN2))
        tf = false;
        return;
    end

    scale = max([1, norm(JN1, 'fro'), norm(JN2, 'fro')]);
    tf = norm(JN1 - JN2, 'fro') <= 1e-10 * scale;
end


function dims = infer_square_channel_dims(JN1, JN2, JT)

    dA = infer_square_channel_dim(JN1, 'JN1');
    dB = infer_square_channel_dim(JN2, 'JN2');
    dC = infer_square_channel_dim(JT, 'JT');
    dims = [dA dA dB dB dC dC];
end


function d = infer_square_channel_dim(J, name)

    n = size(J, 1);
    if size(J, 2) ~= n
        error('%s must be square.', name);
    end

    d = round(sqrt(n));
    if d * d ~= n
        error(['Cannot infer square-channel dimension from %s of size %d. ', ...
            'Pass ''Dims'', [dAi dAo dBi dBo dCi dCo].'], name, n);
    end
end


function validate_dims(dims, JN1, JN2, JT)

    if numel(dims) ~= 6
        error('Dims must be [dAi dAo dBi dBo dCi dCo].');
    end

    expected = [prod(dims(1:2)), prod(dims(3:4)), prod(dims(5:6))];
    actual = [size(JN1, 1), size(JN2, 1), size(JT, 1)];

    if any(actual ~= expected) || size(JN1, 2) ~= expected(1) || ...
            size(JN2, 2) ~= expected(2) || size(JT, 2) ~= expected(3)
        error('Choi matrix sizes do not match Dims.');
    end
end


function print_distance_table(results)

    fprintf('\nDistance table:\n');
    fprintf('%10s ', 'class');
    for r = 1:numel(results.free_classes)
        fprintf('%18s', results.free_classes{r});
    end
    fprintf('\n');

    for s = 1:numel(results.strategies)
        fprintf('%10s ', results.strategies{s});
        for r = 1:numel(results.free_classes)
            fprintf('%18.12g', results.distances(s, r));
        end
        fprintf('\n');
    end
    fprintf('\n');
end


function [distances, adjustments] = apply_hierarchy_upper_bounds( ...
    distances, solutions, strategies, free_classes, verbose)

    adjustments = struct('free_class', {}, 'from_strategy', {}, ...
        'to_strategy', {}, 'raw_value', {}, 'adjusted_value', {}, ...
        'feasibility_residual', {});

    hierarchy_tol = 1e-8;
    feasibility_tol = 1e-7;

    for r = 1:numel(free_classes)
        free_class = free_classes{r};

        [distances, adjustments] = adjust_one_hierarchy_edge( ...
            distances, adjustments, solutions, strategies, ...
            free_class, 'Par', 'Seq', r, hierarchy_tol, ...
            feasibility_tol, verbose);

        [distances, adjustments] = adjust_one_hierarchy_edge( ...
            distances, adjustments, solutions, strategies, ...
            free_class, 'Seq', 'ICO', r, hierarchy_tol, ...
            feasibility_tol, verbose);
    end
end


function [distances, adjustments] = adjust_one_hierarchy_edge( ...
    distances, adjustments, solutions, strategies, ...
    free_class, source_strategy, target_strategy, free_idx, ...
    hierarchy_tol, feasibility_tol, verbose)

    source_idx = find(strcmpi(strategies, source_strategy), 1);
    target_idx = find(strcmpi(strategies, target_strategy), 1);
    if isempty(source_idx) || isempty(target_idx) || ...
            ~isfield(solutions, free_class) || ...
            ~isfield(solutions.(free_class), source_strategy)
        return;
    end

    source_value = distances(source_idx, free_idx);
    target_value = distances(target_idx, free_idx);
    if isnan(source_value) || isnan(target_value) || ...
            target_value <= source_value + hierarchy_tol
        return;
    end

    source_sol = solutions.(free_class).(source_strategy);
    residual = Inf;
    if isfield(source_sol, 'causal_residual_by_strategy') && ...
            isfield(source_sol.causal_residual_by_strategy, target_strategy)
        residual = source_sol.causal_residual_by_strategy.(target_strategy);
    end

    if residual > feasibility_tol
        if verbose
            warning(['%s %s raw value is larger than %s, but the %s ', ...
                'solution was not verified feasible for %s; leaving raw ', ...
                'value unchanged. Residual %.3e.'], ...
                free_class, target_strategy, source_strategy, ...
                source_strategy, target_strategy, residual);
        end
        return;
    end

    if verbose
        fprintf(['\nHierarchy correction: %s %s raw %.12g > %s %.12g. ', ...
            'Using the verified %s feasible point as an upper bound ', ...
            'for %s (residual %.3e).\n'], ...
            free_class, target_strategy, target_value, source_strategy, ...
            source_value, source_strategy, target_strategy, residual);
    end

    adj = struct();
    adj.free_class = free_class;
    adj.from_strategy = source_strategy;
    adj.to_strategy = target_strategy;
    adj.raw_value = target_value;
    adj.adjusted_value = source_value;
    adj.feasibility_residual = residual;
    adjustments(end+1) = adj; %#ok<AGROW>

    distances(target_idx, free_idx) = source_value;
end


%% ========================================================================
%%                              SDP SOLVER
%% ========================================================================

function sol = solve_nonreduced_superchannel_primal_with_retries( ...
    Jresource, JT, dims, strategy, free_class, enforce_slot_symmetry, ...
    use_two_sided_diamond, solver_name, precision_setting, verbose)

    precision_attempts = {precision_setting, '', 'medium'};
    sol = [];
    attempted_precision = {};

    for a = 1:numel(precision_attempts)
        attempt_precision = precision_attempts{a};
        if any(strcmp(attempt_precision, attempted_precision))
            continue;
        end
        attempted_precision{end+1} = attempt_precision; %#ok<AGROW>

        if a > 1 && verbose
            if isempty(attempt_precision)
                fprintf('  retrying with CVX default precision...\n');
            else
                fprintf('  retrying with CVX precision %s...\n', attempt_precision);
            end
        end

        sol = solve_nonreduced_superchannel_primal( ...
            Jresource, JT, dims, strategy, free_class, enforce_slot_symmetry, ...
            use_two_sided_diamond, solver_name, attempt_precision);

        if status_is_solved(sol.status)
            return;
        end
    end
end


function sol = solve_nonreduced_superchannel_primal( ...
    Jresource, JT, dims, strategy, free_class, enforce_slot_symmetry, ...
    use_two_sided_diamond, solver_name, precision_setting)

    strategy = validatestring(strategy, {'Par', 'Seq', 'ICO'});
    free_class = validatestring(free_class, {'MISC', 'DISC'});

    cvx_clear;
    if exist('cvx_solver', 'file') && ~isempty(solver_name)
        cvx_solver(solver_name);
    end

    dims = double(dims(:).');  % [Ai Ao Bi Bo Ci Co]
    dAB = prod(dims(1:4));
    dC = prod(dims(5:6));
    dCi = dims(5);
    Dtot = prod(dims);

    Jresource = herm_numeric(Jresource);
    JT = herm_numeric(JT);
    slot_swaps = slot_symmetry_swaps(dims, strategy, enforce_slot_symmetry);

    cvx_begin sdp quiet
        if ~isempty(precision_setting)
            cvx_precision(precision_setting)
        end

        variable Theta(Dtot, Dtot) hermitian
        variable Z(dC, dC) hermitian
        variable ell nonnegative

        expression Jraw(dC, dC)
        expression Jout(dC, dC)
        expression Diff(dC, dC)

        Jraw = ptrace_cvx(kron(Jresource.', eye(dC)) * Theta, [dAB dC], 1);
        Jout = herm_expr(Jraw);
        Diff = herm_expr(Jout - JT);

        [causal_lhs, causal_rhs] = causal_constraint_cells(Theta, dims, strategy);
        [free_lhs, free_rhs] = free_constraint_cells(Theta, dims, free_class);

        if use_two_sided_diamond
            minimize(0.5 * ell)
        else
            minimize(ell)
        end

        subject to
            Theta >= 0;
            Z >= 0;

            for c = 1:numel(causal_lhs)
                causal_lhs{c} == causal_rhs{c};
            end

            for c = 1:numel(free_lhs)
                free_lhs{c} == free_rhs{c};
            end

            for c = 1:numel(slot_swaps)
                Theta == slot_swaps{c} * Theta * slot_swaps{c}';
            end

            Z - Diff >= 0;
            if use_two_sided_diamond
                Z + Diff >= 0;
            end
            ell * eye(dCi) - ptrace_cvx(Z, dims(5:6), 2) >= 0;
    cvx_end

    sol = collect_solution( ...
        Theta, Z, ell, Jresource, JT, dims, strategy, free_class, ...
        slot_swaps, cvx_status, cvx_optval);

    cvx_clear;
end


%% ========================================================================
%%                              CONSTRAINTS
%% ========================================================================

function [lhs, rhs] = causal_constraint_cells(Theta, dims, strategy)

    lhs = {};
    rhs = {};
    dAi = dims(1);
    dAo = dims(2);
    dBi = dims(3);
    dBo = dims(4);
    dCi = dims(5);

    if strcmpi(strategy, 'Par')
        W = ptrace_cvx(Theta, dims, 6);
        dims_W = dims([1 2 3 4 5]);
        lhs{end+1} = W;
        core = ptrace_cvx(W, dims_W, [2 4]);
        rhs{end+1} = insert_identity_cvx( ...
            core, dims_W, [1 3 5], [2 4], 1 / (dAo * dBo));

        lhs{end+1} = 1 / (dAo * dBo) * ...
            ptrace_cvx(Theta, dims, [1 2 3 4 6]);
        rhs{end+1} = eye(dCi);

        lhs{end+1} = ptrace_cvx(Theta, dims, [1 3 6]);
        rhs{end+1} = eye(dAo * dBo * dCi);

    elseif strcmpi(strategy, 'Seq')
        W = ptrace_cvx(Theta, dims, 6);
        dims_W = dims([1 2 3 4 5]);
        lhs{end+1} = W;
        core = ptrace_cvx(W, dims_W, 4);
        rhs{end+1} = insert_identity_cvx( ...
            core, dims_W, [1 2 3 5], 4, 1 / dBo);

        dims_AiAoCi = dims([1 2 5]);
        lhs{end+1} = 1 / dBi * ptrace_cvx(Theta, dims, [3 4 6]);
        core = ptrace_cvx(Theta, dims, [2 3 4 6]);
        rhs{end+1} = insert_identity_cvx( ...
            core, dims_AiAoCi, [1 3], 2, 1 / (dBi * dAo));

        lhs{end+1} = 1 / (dAo * dBo) * ...
            ptrace_cvx(Theta, dims, [1 2 3 4 6]);
        rhs{end+1} = eye(dCi);

    elseif strcmpi(strategy, 'ICO')
        W = ptrace_cvx(Theta, dims, 6);
        dims_W = dims([1 2 3 4 5]);
        lhs{end+1} = W;
        coreA = ptrace_cvx(W, dims_W, 2);
        termA = insert_identity_cvx( ...
            coreA, dims_W, [1 3 4 5], 2, 1 / dAo);
        coreB = ptrace_cvx(W, dims_W, 4);
        termB = insert_identity_cvx( ...
            coreB, dims_W, [1 2 3 5], 4, 1 / dBo);
        coreAB = ptrace_cvx(W, dims_W, [2 4]);
        termAB = insert_identity_cvx( ...
            coreAB, dims_W, [1 3 5], [2 4], 1 / (dAo * dBo));
        rhs{end+1} = termA + termB - termAB;

        dims_BiBoCi = dims([3 4 5]);
        lhs{end+1} = ptrace_cvx(Theta, dims, [1 2 6]);
        core = ptrace_cvx(Theta, dims, [1 2 4 6]);
        rhs{end+1} = insert_identity_cvx( ...
            core, dims_BiBoCi, [1 3], 2, 1 / dBo);

        dims_AiAoCi = dims([1 2 5]);
        lhs{end+1} = ptrace_cvx(Theta, dims, [3 4 6]);
        core = ptrace_cvx(Theta, dims, [2 3 4 6]);
        rhs{end+1} = insert_identity_cvx( ...
            core, dims_AiAoCi, [1 3], 2, 1 / dAo);

        lhs{end+1} = ptrace_cvx(Theta, dims, [1 2 3 4 6]);
        rhs{end+1} = (dAi * dBi) * eye(dCi);

    else
        error('Unknown strategy.');
    end
end


function [lhs, rhs] = free_constraint_cells(Theta, dims, free_class)

    lhs = {};
    rhs = {};

    if strcmpi(free_class, 'MISC')
        lhs{end+1} = dephase_cvx(Theta, dims, 1:6);
        rhs{end+1} = dephase_cvx(Theta, dims, 1:4);
    elseif strcmpi(free_class, 'DISC')
        lhs{end+1} = dephase_cvx(Theta, dims, [5 6]);
        rhs{end+1} = dephase_cvx(Theta, dims, 1:4);
    else
        error('Unknown free class.');
    end
end


%% ========================================================================
%%                              COLLECTION
%% ========================================================================

function sol = collect_solution( ...
    Theta, Z, ell, Jresource, JT, dims, strategy, free_class, ...
    slot_swaps, status, optval)

    sol = struct();
    sol.strategy = strategy;
    sol.free_class = free_class;
    sol.status = status;
    sol.cvx_optval = optval;
    dC = prod(dims(5:6));

    if status_is_solved(status)
        Theta_value = full(Theta);
        Z_value = full(Z);
        ell_value = full(ell);

        Jout_value = link_output_numeric(Theta_value, Jresource, dims);
        Diff_value = herm_numeric(Jout_value - JT);

        [causal_lhs, causal_rhs] = causal_constraint_cells( ...
            Theta_value, dims, strategy);
        causal_residual_by_strategy = struct();
        for h = {'Par', 'Seq', 'ICO'}
            [test_lhs, test_rhs] = causal_constraint_cells( ...
                Theta_value, dims, h{1});
            causal_residual_by_strategy.(h{1}) = ...
                max(cell_residual_fro(test_lhs, test_rhs));
        end

        sol.distance = optval;
        sol.ell = ell_value;
        sol.Theta = Theta_value;
        sol.Z = Z_value;
        sol.Jout = Jout_value;
        sol.Diff = Diff_value;
        sol.min_eig_Theta = min_eig_herm(Theta_value);
        sol.min_eig_Z = min_eig_herm(Z_value);
        sol.min_eig_Z_minus_Diff = min_eig_herm(Z_value - Diff_value);
        sol.min_eig_Z_plus_Diff = min_eig_herm(Z_value + Diff_value);
        sol.min_eig_ell = min_eig_herm( ...
            ell_value * eye(dims(5)) - ptrace_cvx(Z_value, dims(5:6), 2));
        sol.max_causal_residual = max(cell_residual_fro(causal_lhs, causal_rhs));
        sol.causal_residual_by_strategy = causal_residual_by_strategy;
        sol.free_residual = free_constraint_residual( ...
            Theta_value, dims, free_class);
        sol.misc_residual = free_constraint_residual(Theta_value, dims, 'MISC');
        sol.disc_residual = free_constraint_residual(Theta_value, dims, 'DISC');
        sol.max_slot_residual = slot_symmetry_residual(Theta_value, slot_swaps);
    else
        sol.distance = NaN;
        sol.ell = NaN;
        sol.Jout = NaN(dC, dC);
        sol.Diff = NaN(dC, dC);
        sol.min_eig_Theta = NaN;
        sol.min_eig_Z = NaN;
        sol.min_eig_Z_minus_Diff = NaN;
        sol.min_eig_Z_plus_Diff = NaN;
        sol.min_eig_ell = NaN;
        sol.max_causal_residual = NaN;
        sol.free_residual = NaN;
        sol.misc_residual = NaN;
        sol.disc_residual = NaN;
        sol.max_slot_residual = NaN;
    end
end


function Jout = link_output_numeric(Theta, Jresource, dims)

    dAB = prod(dims(1:4));
    dC = prod(dims(5:6));
    Jout = ptrace_cvx(kron(Jresource.', eye(dC)) * Theta, [dAB dC], 1);
    Jout = herm_numeric(Jout);
end


%% ========================================================================
%%                              HELPERS
%% ========================================================================

function enforce_here = slot_symmetry_setting( ...
    strategy, enforce_slot_symmetry, enforce_sequential_slot_symmetry)

    if strcmpi(strategy, 'Seq')
        enforce_here = enforce_sequential_slot_symmetry;
    else
        enforce_here = enforce_slot_symmetry;
    end
end


function slot_swaps = slot_symmetry_swaps(dims, strategy, enforce_slot_symmetry)

    slot_swaps = {};

    if ~enforce_slot_symmetry
        return;
    end

    if ~any(strcmpi(strategy, {'Par', 'Seq', 'ICO'}))
        error('Unknown strategy.');
    end

    P = complete_slot_swap_unitary( ...
        dims, ...
        {'Ai', 'Ao', 'Bi', 'Bo', 'Ci', 'Co'}, ...
        {'Ai', 'Ao'}, {'Bi', 'Bo'});
    validate_slot_swap_operator(P, prod(dims));

    slot_swaps{end+1} = P;
end


function ok = status_is_solved(status)

    ok = contains(status, 'Solved') || contains(status, 'Inaccurate/Solved');
end


function Y = ptrace_cvx(X, dims, traced)

    dims = double(dims(:).');
    n = numel(dims);
    traced = sort(traced(:).');
    keep = setdiff(1:n, traced, 'stable');

    D = prod(dims);
    d_keep = prod(dims(keep));
    subs = all_subs_kron(dims);

    Y = X(1, 1) * zeros(d_keep, d_keep);
    for r = 1:D
        sr = subs(r, :);
        rk = sub_to_ind_kron(dims(keep), sr(keep));

        for c = 1:D
            sc = subs(c, :);
            if all(sr(traced) == sc(traced))
                ck = sub_to_ind_kron(dims(keep), sc(keep));
                Y(rk, ck) = Y(rk, ck) + X(r, c);
            end
        end
    end
end


function Y = insert_identity_cvx(T, full_dims, keep_sys, id_sys, scale)

    full_dims = double(full_dims(:).');
    keep_sys = double(keep_sys(:).');
    id_sys = double(id_sys(:).');

    temp_order = [keep_sys, id_sys];
    dims_temp = full_dims(temp_order);
    Did = prod(full_dims(id_sys));

    Ttemp = kron(T, scale * speye(Did));

    new_order = zeros(1, numel(full_dims));
    for s = 1:numel(full_dims)
        new_order(s) = find(temp_order == s);
    end

    Y = reorder_operator_cvx(Ttemp, dims_temp, new_order);
end


function Y = reorder_operator_cvx(X, dims_old, new_order)

    P = permutation_unitary_kron(dims_old, new_order);
    Y = P * X * P';
end


function Y = dephase_cvx(X, dims, subs)

    mask = dephase_mask(dims, subs);
    Y = mask .* X;
end


function mask = dephase_mask(dims, subs)

    dims = double(dims(:).');
    N = prod(dims);
    labels = all_subs_kron(dims);
    mask = zeros(N, N);

    for r = 1:N
        for c = 1:N
            if all(labels(r, subs) == labels(c, subs))
                mask(r, c) = 1;
            end
        end
    end
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
        sub_in = ind_to_sub_kron(dims, col);
        sub_out = sub_in(permutation);
        row_idx(col) = sub_to_ind_kron(dims_out, sub_out);
    end

    P = sparse(row_idx, col_idx, ones(D, 1), D, D);
end


function subs = all_subs_kron(dims)

    dims = double(dims(:).');
    D = prod(dims);

    if isempty(dims)
        subs = zeros(1, 0);
        return;
    end

    subs = zeros(D, numel(dims));
    for idx = 1:D
        subs(idx, :) = ind_to_sub_kron(dims, idx);
    end
end


function sub = ind_to_sub_kron(dims, idx)

    dims = double(dims(:).');
    sub = zeros(1, numel(dims));
    idx0 = idx - 1;

    for k = numel(dims):-1:1
        sub(k) = mod(idx0, dims(k)) + 1;
        idx0 = floor(idx0 / dims(k));
    end
end


function idx = sub_to_ind_kron(dims, sub)

    dims = double(dims(:).');
    sub = double(sub(:).');

    idx0 = 0;
    for k = 1:numel(dims)
        idx0 = idx0 * dims(k) + (sub(k) - 1);
    end

    idx = idx0 + 1;
end


function H = herm_expr(X)

    H = 0.5 * (X + X');
end


function H = herm_numeric(X)

    H = 0.5 * (X + X');
end


function residuals = cell_residual_fro(lhs, rhs)

    residuals = zeros(1, numel(lhs));
    for c = 1:numel(lhs)
        residuals(c) = norm(full(lhs{c} - rhs{c}), 'fro');
    end
end


function residual = free_constraint_residual(Theta, dims, free_class)

    if strcmpi(free_class, 'MISC')
        residual = norm( ...
            dephase_cvx(Theta, dims, 1:6) ...
            - dephase_cvx(Theta, dims, 1:4), 'fro');
    elseif strcmpi(free_class, 'DISC')
        residual = norm( ...
            dephase_cvx(Theta, dims, [5 6]) ...
            - dephase_cvx(Theta, dims, 1:4), 'fro');
    else
        error('Unknown free class.');
    end
end


function residual = slot_symmetry_residual(Theta, slot_swaps)

    residual = 0;
    for c = 1:numel(slot_swaps)
        residual = max(residual, ...
            norm(Theta - slot_swaps{c} * Theta * slot_swaps{c}', 'fro'));
    end
end


function value = min_eig_herm(X)

    X = herm_numeric(full(X));
    value = min(real(eig(X)));
end
