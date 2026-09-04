%% verify_saved_results.m
% Verify the numerical statements reported with the supplied SeDuMi data.

clear; clc;

data_dir = fullfile('results', 'data');
primal_file = fullfile(data_dir, ...
    'AD_to_identity_reduced_misc_primal.mat');
dual_file = fullfile(data_dir, ...
    'AD_to_identity_reduced_misc_dual.mat');
random_file = fullfile(data_dir, ...
    ['random_resource_pairs_identity_target_MISC_Seq_ICO_', ...
     'gap_distribution_sedumi_1000_primal_dual.mat']);

required_files = {primal_file, dual_file, random_file};
for k = 1:numel(required_files)
    if exist(required_files{k}, 'file') ~= 2
        error('Missing supplied result file: %s', required_files{k});
    end
end

primal = load(primal_file);
dual = load(dual_file);

assert(isequal(primal.error_rates, (0:19) / 20), ...
    'Unexpected amplitude-damping grid in the primal data.');
assert(isequal(primal.error_rates, dual.error_rates), ...
    'Primal and dual amplitude-damping grids differ.');
assert(isequal(primal.classes_to_solve, dual.classes_to_solve), ...
    'Primal and dual causal-class lists differ.');

ad_mismatches = zeros(numel(primal.error_rates), ...
    numel(primal.classes_to_solve));
for e = 1:numel(primal.error_rates)
    for h = 1:numel(primal.classes_to_solve)
        causal_class = primal.classes_to_solve{h};
        primal_value = primal.sweep(e).solutions.(causal_class).distance;
        dual_value = dual.sweep(e).solutions.(causal_class).lower_bound;
        ad_mismatches(e, h) = abs(primal_value - dual_value);
    end
end

[ad_max_mismatch, ad_linear_index] = max(ad_mismatches, [], 'all');
[ad_error_index, ad_class_index] = ind2sub(size(ad_mismatches), ...
    ad_linear_index);

random = load(random_file);
valid = ~random.failed_sets & isfinite(random.seq_minus_ico);
random_mismatches = [abs(random.seq_primal_dual_signed_gap(valid)); ...
    abs(random.ico_primal_dual_signed_gap(valid))];
random_max_mismatch = max(random_mismatches);
above_resolution = sum(random.seq_minus_ico(valid) >= random.gap_min);
below_resolution = sum(random.seq_minus_ico(valid) < random.gap_min);

assert(random.num_sets == 1000, 'The supplied random run is not 1000 samples.');
assert(sum(random.failed_sets) == 0, 'The supplied random run has failed samples.');
assert(sum(valid) == 1000, 'The supplied random run has non-finite gaps.');
assert(above_resolution == 989, 'Expected 989 gaps at least 1e-5.');
assert(below_resolution == 11, 'Expected 11 gaps below 1e-5.');
assert(ad_max_mismatch < 1e-8, ...
    'Amplitude-damping primal/dual mismatch is not below 1e-8.');
assert(random_max_mismatch < 1e-7, ...
    'Random-study primal/dual mismatch is not below 1e-7.');

fprintf('Amplitude-damping grid points       : %d\n', ...
    numel(primal.error_rates));
fprintf('Amplitude-damping max |primal-dual|: %.15g\n', ...
    ad_max_mismatch);
fprintf('  location: epsilon %.2f, class %s\n', ...
    primal.error_rates(ad_error_index), ...
    primal.classes_to_solve{ad_class_index});
fprintf('Random pairs completed              : %d\n', sum(valid));
fprintf('Gaps at least %.1e                  : %d\n', ...
    random.gap_min, above_resolution);
fprintf('Gaps below %.1e                     : %d\n', ...
    random.gap_min, below_resolution);
fprintf('Random-study max |primal-dual|      : %.15g\n', ...
    random_max_mismatch);
fprintf('All supplied numerical checks passed.\n');
