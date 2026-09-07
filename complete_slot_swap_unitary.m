function [P, permutation] = complete_slot_swap_unitary( ...
    register_dims, register_order, slot_a, slot_b, basis_order)
%COMPLETE_SLOT_SWAP_UNITARY Swap two complete channel slots.
%
%   [P, permutation] = complete_slot_swap_unitary(register_dims,
%   register_order, slot_a, slot_b) swaps every subsystem label in slot_a
%   with the corresponding subsystem label in slot_b. A channel slot should
%   include both its input and output registers, e.g. {'Ai','Ao'}.
%
%   The optional basis_order is passed to subsystem_permutation_unitary.
%   Default 'kron' matches MATLAB kron/tensor and QETLAB. Use
%   'column-major' for code that indexes tensor factors with reshape/sub2ind.

    register_dims = double(register_dims(:).');
    register_order = normalize_labels(register_order, 'register_order');
    slot_a = normalize_labels(slot_a, 'slot_a');
    slot_b = normalize_labels(slot_b, 'slot_b');

    if nargin < 5 || isempty(basis_order)
        basis_order = 'kron';
    end

    n = length(register_dims);

    if length(register_order) ~= n
        error('register_order must have the same length as register_dims.');
    end

    if length(unique(register_order)) ~= length(register_order)
        error('register_order labels must be unique.');
    end

    if length(slot_a) ~= length(slot_b)
        error('The two slots must contain the same number of subsystem labels.');
    end

    if ~isempty(intersect(slot_a, slot_b))
        error('The two slots to be swapped must be disjoint.');
    end

    pos_a = label_positions(register_order, slot_a);
    pos_b = label_positions(register_order, slot_b);

    if any(register_dims(pos_a) ~= register_dims(pos_b))
        error('Corresponding subsystem dimensions in the two slots must match.');
    end

    permutation = 1:n;
    permutation(pos_a) = pos_b;
    permutation(pos_b) = pos_a;

    P = subsystem_permutation_unitary(register_dims, permutation, basis_order);
    validate_slot_swap_operator(P, prod(register_dims));
end


function labels = normalize_labels(labels, name)

    if ischar(labels)
        labels = {labels};
    elseif iscell(labels)
        labels = labels(:).';
    else
        try
            labels = cellstr(labels);
            labels = labels(:).';
        catch
            error('%s must be a char, string array, or cell array of labels.', name);
        end
    end

    for k = 1:length(labels)
        if ~ischar(labels{k})
            try
                labels{k} = char(labels{k});
            catch
                error('%s contains a non-text label.', name);
            end
        end
        labels{k} = strtrim(labels{k});
    end
end


function positions = label_positions(register_order, labels)

    positions = zeros(1, length(labels));

    for k = 1:length(labels)
        match = find(strcmp(register_order, labels{k}));
        if isempty(match)
            error('Unknown subsystem label: %s.', labels{k});
        end
        positions(k) = match;
    end
end
