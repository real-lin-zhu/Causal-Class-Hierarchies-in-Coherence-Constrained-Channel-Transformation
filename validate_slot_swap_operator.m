function validate_slot_swap_operator(P, dim)
%VALIDATE_SLOT_SWAP_OPERATOR Basic checks for a slot-swap permutation matrix.

    if nargin < 2
        dim = size(P, 1);
    end

    if ~isequal(size(P), [dim, dim])
        error('Slot-swap operator must be a %d-by-%d matrix.', dim, dim);
    end

    tol = 1e-10;
    if norm(P * P' - eye(dim), 'fro') > tol || ...
            norm(P' * P - eye(dim), 'fro') > tol
        error('Slot-swap operator is not unitary/permutation-like.');
    end

    if norm(P * P - eye(dim), 'fro') > tol
        error('Slot-swap operator is not an involution.');
    end
end
