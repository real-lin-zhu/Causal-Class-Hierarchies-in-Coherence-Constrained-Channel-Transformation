%% run_general_channel_distance_sdp.m
% Reproduce the numerical MISC-DISC separation example from Appendix H.
%
% The matrices below are the four-decimal values printed in the paper. The
% paper's reported table uses the corresponding full-precision data.

clear; clc;

%% Numerical MISC-DISC separation example

JN1 = [
    0.5717,              -0.0947 + 0.1220i,  -0.0474 + 0.0567i,   0.0592 + 0.0537i;
   -0.0947 - 0.1220i,    0.4283,              -0.0045 + 0.2150i,   0.0474 - 0.0567i;
   -0.0474 - 0.0567i,   -0.0045 - 0.2150i,    0.6068,             -0.0755 - 0.2755i;
    0.0592 - 0.0537i,    0.0474 + 0.0567i,   -0.0755 + 0.2755i,    0.3932
];

JN2 = [
    0.4828,              -0.1180 + 0.0428i,  -0.0511 + 0.0630i,   0.1102 + 0.0068i;
   -0.1180 - 0.0428i,    0.5172,               0.0115 - 0.1401i,   0.0511 - 0.0630i;
   -0.0511 - 0.0630i,    0.0115 + 0.1401i,    0.2717,             -0.0660 - 0.0005i;
    0.1102 - 0.0068i,    0.0511 + 0.0630i,   -0.0660 + 0.0005i,    0.7283
];

% Identity target Choi matrix, J_I.
JT = [
    1, 0, 0, 1;
    0, 0, 0, 0;
    0, 0, 0, 0;
    1, 0, 0, 1
];

dims = [2 2 2 2 2 2];

%% Solver and model choices

strategies = {'Par', 'Seq', 'ICO'};
free_classes = {'MISC', 'DISC'};

use_two_sided_diamond = false;

% [] means auto. Since JN1 and JN2 are different, no slot symmetry is
% imposed on Par or ICO for this example.
enforce_slot_symmetry = [];
enforce_sequential_slot_symmetry = false;
cvx_solver_name = 'sedumi';
cvx_precision_setting = 'best';

%% Run the full non-reduced SDP

args = {
    'Strategies', strategies, ...
    'FreeClasses', free_classes, ...
    'UseTwoSidedDiamond', use_two_sided_diamond, ...
    'EnforceSlotSymmetry', enforce_slot_symmetry, ...
    'EnforceSequentialSlotSymmetry', enforce_sequential_slot_symmetry, ...
    'CvxSolver', cvx_solver_name, ...
    'CvxPrecision', cvx_precision_setting, ...
    'Verbose', true
};

results = general_nonreduced_channel_distance_sdp( ...
    JN1, JN2, JT, args{:});

disp(results.distances);
