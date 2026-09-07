function P = subsystem_permutation_unitary(dims, permutation, basis_order)
%SUBSYSTEM_PERMUTATION_UNITARY Permutation unitary for tensor subsystems.
%
%   P = subsystem_permutation_unitary(dims, permutation) returns the sparse
%   matrix implementing
%
%       |i_1, ..., i_n> -> |i_{permutation(1)}, ..., i_{permutation(n)}>.
%
%   The input basis is ordered by dims, and the output basis is ordered by
%   dims(permutation). The total Hilbert-space dimension is unchanged.
%
%   By default, the unitary uses Kronecker-product basis ordering, i.e. the
%   last tensor factor is fastest. This matches MATLAB kron/tensor and
%   QETLAB. Pass 'column-major' to use MATLAB reshape/sub2ind ordering, i.e.
%   the first tensor factor is fastest.

    dims = double(dims(:).');
    permutation = double(permutation(:).');

    if nargin < 3 || isempty(basis_order)
        basis_order = 'kron';
    end

    if isempty(dims)
        error('dims must contain at least one tensor factor.');
    end

    n = length(dims);

    if length(permutation) ~= n || ~isequal(sort(permutation), 1:n)
        error('permutation must be a permutation of 1:length(dims).');
    end

    if any(dims < 1) || any(abs(dims - round(dims)) > 0)
        error('dims must be positive integer dimensions.');
    end

    D = prod(dims);
    row_idx = zeros(D, 1);
    col_idx = (1:D).';

    dims_out = dims(permutation);

    for col = 1:D
        sub_in = ind_to_subscript_row(dims, col, basis_order);
        sub_out = sub_in(permutation);
        row_idx(col) = subscript_row_to_ind(dims_out, sub_out, basis_order);
    end

    P = sparse(row_idx, col_idx, ones(D, 1), D, D);
end


function sub = ind_to_subscript_row(dims, idx, basis_order)

    if strcmpi(basis_order, 'column-major')
        args = cell(1, length(dims));
        [args{:}] = ind2sub(dims, idx);
        sub = zeros(1, length(dims));

        for k = 1:length(dims)
            sub(k) = args{k};
        end
        return;
    end

    if ~strcmpi(basis_order, 'kron')
        error('Unknown basis_order. Use ''kron'' or ''column-major''.');
    end

    sub = zeros(1, length(dims));
    idx0 = idx - 1;

    for k = length(dims):-1:1
        sub(k) = mod(idx0, dims(k)) + 1;
        idx0 = floor(idx0 / dims(k));
    end
end


function idx = subscript_row_to_ind(dims, sub, basis_order)

    if strcmpi(basis_order, 'column-major')
        args = num2cell(sub);
        idx = sub2ind(dims, args{:});
        return;
    end

    if ~strcmpi(basis_order, 'kron')
        error('Unknown basis_order. Use ''kron'' or ''column-major''.');
    end

    idx0 = 0;

    for k = 1:length(dims)
        idx0 = idx0 * dims(k) + (sub(k) - 1);
    end

    idx = idx0 + 1;
end
