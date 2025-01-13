function out = desentrelazado(bits, Ncbps)
    bits_matrix = reshape(bits, Ncbps, length(bits)/Ncbps);
    out = reshape(bits_matrix', [], 1);
end
