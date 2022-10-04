function varargout = compare_ix(ix_in, ix_ref)

% Creates a type of confusion matrix

ix_min = min(min(ix_in), min(ix_ref));
ix_max = max(max(ix_in), max(ix_ref));

comp_out = zeros(ix_max - ix_min + 1);

for kk = ix_min:ix_max
    for ll = ix_min:ix_max
        comp_out(kk, ll) = sum((ix_in == kk)&(ix_ref==ll));
    end
end

varargout{1} = comp_out;

end
