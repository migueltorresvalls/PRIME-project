function out = entrelazado(bits, Ncbps)
    bits_matrix = reshape(bits, length(bits)/Ncbps, Ncbps);
    out = reshape(bits_matrix', [], 1);
end
