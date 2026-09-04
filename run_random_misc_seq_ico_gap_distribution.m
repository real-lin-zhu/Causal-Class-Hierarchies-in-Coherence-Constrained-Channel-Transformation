%% run_random_misc_seq_ico_gap_distribution.m
% Generate many random resource-channel pairs with the identity target
% channel, then count how often the MISC Seq-ICO distance gap falls into
% bins of width 0.001.
%
% The histogram uses
%   x-axis: increase = Seq distance - ICO distance under MISC
%   y-axis: number of random sets in that increase bin

clear; clc;

addpath(genpath(fullfile(pwd, 'quantinf')));

%% Experiment settings

num_sets = env_positive_integer('RANDOM_MISC_NUM_SETS', 1000);
bin_width = 0.001;
dims = [2 2 2 2 2 2];
rng_seed = 20260623;
num_dual_samples = env_positive_integer('RANDOM_MISC_NUM_DUAL_SAMPLES', 5);
dual_certificate_tolerance = 1e-8;
solve_all_duals = env_logical('RANDOM_MISC_SOLVE_ALL_DUALS', false);
progress_every = env_positive_integer('RANDOM_MISC_PROGRESS_EVERY', 1);
output_tag = strtrim(getenv('RANDOM_MISC_OUTPUT_TAG'));
if ~isempty(output_tag) && isempty(regexp(output_tag, '^[A-Za-z0-9_-]+$', 'once'))
    error('RANDOM_MISC_OUTPUT_TAG may contain only letters, digits, _ and -.');
end

strategies = {'Seq', 'ICO'};
free_classes = {'MISC'};

use_two_sided_diamond = false;
enforce_slot_symmetry = [];
enforce_sequential_slot_symmetry = false;
sedumi_pars = struct('eps', 1e-9, 'bigeps', 1e-7, 'fid', 0);
verbose_sdp = false;

if exist('sedumi', 'file') ~= 2
    error(['SeDuMi was not found on the MATLAB path. Install SeDuMi and ', ...
        'run its setup script before running this file.']);
end

% Full SDP solution structs are large for 1000 samples. Keep this false
% unless you need each Theta/Z solution for later inspection.
save_full_results = strcmpi(getenv('SAVE_FULL_RESULTS'), 'true');

output_data_dir = fullfile(pwd, 'results', 'data');
output_fig_dir = fullfile(pwd, 'results', 'figures');
if ~exist(output_data_dir, 'dir')
    mkdir(output_data_dir);
end
if ~exist(output_fig_dir, 'dir')
    mkdir(output_fig_dir);
end

file_stem = 'random_resource_pairs_identity_target_MISC_Seq_ICO_gap_distribution';
if ~isempty(output_tag)
    file_stem = sprintf('%s_%s', file_stem, output_tag);
end
data_file = fullfile(output_data_dir, [file_stem, '.mat']);
figure_file = fullfile(output_fig_dir, [file_stem, '.png']);

rng(rng_seed, 'twister');

JN1_all = cell(num_sets, 1);
JN2_all = cell(num_sets, 1);
JT = qubit_identity_channel_choi();
results_all = cell(num_sets, 1);

seq_distances = NaN(num_sets, 1);
ico_distances = NaN(num_sets, 1);
seq_minus_ico = NaN(num_sets, 1);
statuses = cell(num_sets, 2);
failed_sets = false(num_sets, 1);
seq_dual_bounds = NaN(num_sets, 1);
ico_dual_bounds = NaN(num_sets, 1);
seq_primal_dual_signed_gap = NaN(num_sets, 1);
ico_primal_dual_signed_gap = NaN(num_sets, 1);
dual_statuses = cell(num_sets, 2);
primal_dual_mismatch_summary = struct();
least_gap_sample_indices = zeros(0, 1);
dual_certificates = struct([]);

fprintf('Running %d random MISC Seq/ICO channel-distance SDPs.\n', num_sets);
fprintf('Each set uses two different random resource channels and the qubit identity target channel.\n');
fprintf('Histogram bin width: %.6g\n\n', bin_width);
fprintf('Explicit dual for every sample: %d\n\n', solve_all_duals);

for set_idx = 1:num_sets
    show_progress = set_idx == 1 || set_idx == num_sets || ...
        mod(set_idx, progress_every) == 0;
    if show_progress
        fprintf('Set %d / %d ... ', set_idx, num_sets);
    end

    [JN1, JN2] = random_distinct_qubit_channel_pair();
    JN1_all{set_idx} = JN1;
    JN2_all{set_idx} = JN2;

    try
        Jresource = kron(JN1, JN2);
        if solve_all_duals
            seq_pair = solve_reduced_misc_sedumi( ...
                Jresource, 'Seq', use_two_sided_diamond, sedumi_pars, 'both');
            ico_pair = solve_reduced_misc_sedumi( ...
                Jresource, 'ICO', use_two_sided_diamond, sedumi_pars, 'both');
            seq_solution = seq_pair.primal;
            ico_solution = ico_pair.primal;
            seq_dual = seq_pair.dual;
            ico_dual = ico_pair.dual;
        else
            seq_solution = solve_reduced_misc_sedumi( ...
                Jresource, 'Seq', use_two_sided_diamond, sedumi_pars, 'primal');
            ico_solution = solve_reduced_misc_sedumi( ...
                Jresource, 'ICO', use_two_sided_diamond, sedumi_pars, 'primal');
        end
        if ~status_is_solved(seq_solution.status) || ...
                ~status_is_solved(ico_solution.status)
            error('SeDuMi primal status: Seq=%s, ICO=%s.', ...
                seq_solution.status, ico_solution.status);
        end

        seq_distances(set_idx) = seq_solution.distance;
        ico_distances(set_idx) = ico_solution.distance;
        seq_minus_ico(set_idx) = seq_distances(set_idx) - ico_distances(set_idx);
        statuses(set_idx, :) = {seq_solution.status, ico_solution.status};

        if solve_all_duals
            seq_dual_bounds(set_idx) = seq_dual.lower_bound;
            ico_dual_bounds(set_idx) = ico_dual.lower_bound;
            seq_primal_dual_signed_gap(set_idx) = ...
                seq_solution.distance - seq_dual.lower_bound;
            ico_primal_dual_signed_gap(set_idx) = ...
                ico_solution.distance - ico_dual.lower_bound;
            dual_statuses(set_idx, :) = {seq_dual.status, ico_dual.status};
        end

        if save_full_results
            results_all{set_idx} = struct( ...
                'Seq', seq_solution, 'ICO', ico_solution);
        end

        if show_progress
            if solve_all_duals
                fprintf('gap %.12g, p-d mismatch Seq %.3e, ICO %.3e\n', ...
                    seq_minus_ico(set_idx), ...
                    abs(seq_primal_dual_signed_gap(set_idx)), ...
                    abs(ico_primal_dual_signed_gap(set_idx)));
            else
                fprintf('gap %.12g\n', seq_minus_ico(set_idx));
            end
        end
    catch ME
        failed_sets(set_idx) = true;
        results_all{set_idx} = struct( ...
            'failed', true, ...
            'identifier', ME.identifier, ...
            'message', ME.message, ...
            'stack', ME.stack);
        fprintf('failed: %s\n', ME.message);
    end

    save(data_file, ...
        'num_sets', 'bin_width', 'dims', 'rng_seed', ...
        'strategies', 'free_classes', ...
        'use_two_sided_diamond', ...
        'enforce_slot_symmetry', ...
        'enforce_sequential_slot_symmetry', ...
        'sedumi_pars', 'num_dual_samples', 'dual_certificate_tolerance', ...
        'solve_all_duals', 'progress_every', ...
        'output_tag', ...
        'verbose_sdp', 'save_full_results', ...
        'JN1_all', 'JN2_all', 'JT', 'results_all', ...
        'seq_distances', 'ico_distances', 'seq_minus_ico', ...
        'statuses', 'failed_sets', ...
        'seq_dual_bounds', 'ico_dual_bounds', ...
        'seq_primal_dual_signed_gap', 'ico_primal_dual_signed_gap', ...
        'dual_statuses', 'primal_dual_mismatch_summary', ...
        'least_gap_sample_indices', 'dual_certificates');
end

if solve_all_duals
    [primal_dual_mismatch_summary, mismatch_text] = ...
        summarize_primal_dual_mismatches( ...
        seq_primal_dual_signed_gap, ico_primal_dual_signed_gap, failed_sets);
    fprintf('\n%s\n', mismatch_text);
end

%% Explicit dual SDPs for the smallest nonnegative primal gaps

dual_candidates = find(~failed_sets & isfinite(seq_minus_ico) & ...
    seq_minus_ico >= 0);
[~, dual_order] = sort(seq_minus_ico(dual_candidates), 'ascend');
dual_count = min(num_dual_samples, numel(dual_candidates));
least_gap_sample_indices = dual_candidates(dual_order(1:dual_count));
dual_certificates = repmat(struct( ...
    'sample_index', [], 'primal_gap', [], ...
    'seq_primal', [], 'ico_primal', [], ...
    'seq_dual', [], 'ico_dual', [], ...
    'seq_primal_dual_gap', [], 'ico_primal_dual_gap', [], ...
    'certified_seq_minus_ico_lower', [], ...
    'certifies_positive_gap', false, 'failed', false, ...
    'error_identifier', '', 'error_message', ''), dual_count, 1);

fprintf('\nSolving explicit dual SDPs for %d least nonnegative-gap samples.\n', ...
    dual_count);
for k = 1:dual_count
    set_idx = least_gap_sample_indices(k);
    cert = dual_certificates(k);
    cert.sample_index = set_idx;
    cert.primal_gap = seq_minus_ico(set_idx);
    cert.seq_primal = seq_distances(set_idx);
    cert.ico_primal = ico_distances(set_idx);
    try
        Jresource = kron(JN1_all{set_idx}, JN2_all{set_idx});
        if solve_all_duals
            cert.seq_dual = struct('status', dual_statuses{set_idx, 1}, ...
                'lower_bound', seq_dual_bounds(set_idx));
            cert.ico_dual = struct('status', dual_statuses{set_idx, 2}, ...
                'lower_bound', ico_dual_bounds(set_idx));
        else
            cert.seq_dual = solve_reduced_misc_sedumi( ...
                Jresource, 'Seq', use_two_sided_diamond, sedumi_pars, 'dual');
            cert.ico_dual = solve_reduced_misc_sedumi( ...
                Jresource, 'ICO', use_two_sided_diamond, sedumi_pars, 'dual');
        end
        cert.seq_primal_dual_gap = ...
            cert.seq_primal - cert.seq_dual.lower_bound;
        cert.ico_primal_dual_gap = ...
            cert.ico_primal - cert.ico_dual.lower_bound;
        cert.certified_seq_minus_ico_lower = ...
            cert.seq_dual.lower_bound - cert.ico_primal;
        cert.certifies_positive_gap = ...
            cert.certified_seq_minus_ico_lower > dual_certificate_tolerance;
        fprintf(['  set %d: primal gap %.6e, Seq p-d %.3e, ', ...
            'ICO p-d %.3e, certified lower gap %.6e\n'], ...
            set_idx, cert.primal_gap, cert.seq_primal_dual_gap, ...
            cert.ico_primal_dual_gap, cert.certified_seq_minus_ico_lower);
    catch ME
        cert.failed = true;
        cert.error_identifier = ME.identifier;
        cert.error_message = ME.message;
        fprintf('  set %d dual failed: %s\n', set_idx, ME.message);
    end
    dual_certificates(k) = cert;
end

%% Count gap intervals and plot

valid_gaps = seq_minus_ico(~isnan(seq_minus_ico));
if isempty(valid_gaps)
    error('No successful sets were solved; no histogram was generated.');
end

gap_min = 1e-5;
negative_count = sum(valid_gaps < 0);
excluded_below_gap_min_count = sum(valid_gaps < gap_min);
histogram_gaps = valid_gaps(valid_gaps >= gap_min);
if isempty(histogram_gaps)
    error('No gaps are at least %.6g; no histogram was generated.', gap_min);
end

max_gap = max(histogram_gaps);
max_edge = max(bin_width, ceil(max_gap / bin_width) * bin_width);
bin_edges = 0:bin_width:max_edge;
if bin_edges(end) <= max_gap
    bin_edges(end + 1) = bin_edges(end) + bin_width;
end
bin_edges(1) = gap_min;

bin_counts = histcounts(histogram_gaps, bin_edges);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

figure('Color', 'w');
bar(bin_centers, bin_counts, 1.0);
grid on;
xlabel('Increase: Seq distance - ICO distance under MISC');
ylabel('Number of sets');
title(sprintf('Distribution of MISC Seq-ICO Increase Across %d Random Sets', num_sets));
xlim([gap_min, bin_edges(end)]);
automatic_xticks = xticks;
automatic_xticks = automatic_xticks(automatic_xticks > gap_min & ...
    automatic_xticks <= bin_edges(end));
xticks([gap_min, automatic_xticks]);
automatic_xticklabels = arrayfun(@(tick) sprintf('%g', tick), ...
    automatic_xticks, 'UniformOutput', false);
xticklabels([{'10^{-5}'}, automatic_xticklabels]);

exportgraphics(gcf, figure_file, 'Resolution', 200);

save(data_file, ...
    'num_sets', 'bin_width', 'dims', 'rng_seed', ...
    'strategies', 'free_classes', ...
    'use_two_sided_diamond', ...
    'enforce_slot_symmetry', ...
    'enforce_sequential_slot_symmetry', ...
    'sedumi_pars', 'num_dual_samples', 'dual_certificate_tolerance', ...
    'solve_all_duals', 'progress_every', ...
    'output_tag', ...
    'verbose_sdp', 'save_full_results', ...
    'JN1_all', 'JN2_all', 'JT', 'results_all', ...
    'seq_distances', 'ico_distances', 'seq_minus_ico', ...
    'statuses', 'failed_sets', ...
    'seq_dual_bounds', 'ico_dual_bounds', ...
    'seq_primal_dual_signed_gap', 'ico_primal_dual_signed_gap', ...
    'dual_statuses', 'primal_dual_mismatch_summary', ...
    'least_gap_sample_indices', 'dual_certificates', ...
    'valid_gaps', 'gap_min', 'negative_count', ...
    'excluded_below_gap_min_count', 'bin_edges', ...
    'bin_centers', 'bin_counts');

fprintf('\nHistogram counts for gaps at least %.6g:\n', gap_min);
for k = 1:numel(bin_counts)
    fprintf('[%.6g, %.6g): %d\n', bin_edges(k), bin_edges(k + 1), bin_counts(k));
end
fprintf('Negative gap count: %d\n', negative_count);
fprintf('Gap count below %.6g (excluded): %d\n', ...
    gap_min, excluded_below_gap_min_count);
fprintf('\nSaved data to:\n  %s\n', data_file);
fprintf('Saved plot to:\n  %s\n', figure_file);

%% Local helpers

function [JN1, JN2] = random_distinct_qubit_channel_pair()

    distinct_tol = 1e-10;
    max_attempts = 100;

    JN1 = random_qubit_channel_choi();
    for attempt = 1:max_attempts
        JN2 = random_qubit_channel_choi();
        scale = max([1, norm(JN1, 'fro'), norm(JN2, 'fro')]);
        if norm(JN1 - JN2, 'fro') > distinct_tol * scale
            return;
        end
    end

    error('Failed to generate two distinct random resource channels.');
end


function J = random_qubit_channel_choi()

    J = randChan(2, 'choi');
    J = 0.5 * (full(J) + full(J)');
end


function J = qubit_identity_channel_choi()

    bell_unnormalized = [1; 0; 0; 1];
    J = bell_unnormalized * bell_unnormalized';
end


function tf = status_is_solved(status)

    tf = strcmp(status, 'Solved') || strcmp(status, 'Inaccurate/Solved');
end


function value = env_positive_integer(name, default_value)

    text = strtrim(getenv(name));
    if isempty(text)
        value = default_value;
        return;
    end
    value = str2double(text);
    if ~isfinite(value) || value < 1 || value ~= floor(value)
        error('%s must be a positive integer.', name);
    end
end


function value = env_logical(name, default_value)

    text = lower(strtrim(getenv(name)));
    if isempty(text)
        value = logical(default_value);
    elseif any(strcmp(text, {'1', 'true', 'yes', 'on'}))
        value = true;
    elseif any(strcmp(text, {'0', 'false', 'no', 'off'}))
        value = false;
    else
        error('%s must be true or false.', name);
    end
end


function [summary, text] = summarize_primal_dual_mismatches( ...
    seq_signed, ico_signed, failed_sets)

    seq_samples = find(~failed_sets & isfinite(seq_signed));
    ico_samples = find(~failed_sets & isfinite(ico_signed));
    signed = [seq_signed(seq_samples); ico_signed(ico_samples)];
    samples = [seq_samples; ico_samples];
    classes = [repmat({'Seq'}, numel(seq_samples), 1); ...
        repmat({'ICO'}, numel(ico_samples), 1)];
    if isempty(signed)
        error('No finite primal-dual comparisons were produced.');
    end

    mismatch = abs(signed);
    [minimum, imin] = min(mismatch);
    [maximum, imax] = max(mismatch);
    summary = struct();
    summary.num_comparisons = numel(mismatch);
    summary.minimum_absolute_mismatch = minimum;
    summary.minimum_sample_index = samples(imin);
    summary.minimum_class = classes{imin};
    summary.minimum_signed_gap = signed(imin);
    summary.maximum_absolute_mismatch = maximum;
    summary.maximum_sample_index = samples(imax);
    summary.maximum_class = classes{imax};
    summary.maximum_signed_gap = signed(imax);
    summary.mean_absolute_mismatch = mean(mismatch);
    summary.median_absolute_mismatch = median(mismatch);

    text = sprintf([ ...
        'Primal-dual comparisons: %d\n', ...
        'Smallest absolute mismatch: %.12g (%s sample %d, signed %.12g)\n', ...
        'Largest absolute mismatch : %.12g (%s sample %d, signed %.12g)\n', ...
        'Mean / median mismatch    : %.12g / %.12g'], ...
        summary.num_comparisons, ...
        summary.minimum_absolute_mismatch, summary.minimum_class, ...
        summary.minimum_sample_index, summary.minimum_signed_gap, ...
        summary.maximum_absolute_mismatch, summary.maximum_class, ...
        summary.maximum_sample_index, summary.maximum_signed_gap, ...
        summary.mean_absolute_mismatch, summary.median_absolute_mismatch);
end
