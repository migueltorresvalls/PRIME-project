function [out, piloto] = mod_ofdm_cp(Nfft, Nofdm, Nf, L_cp, x_mod)
    s = reshape(x_mod, Nf, Nofdm);
    X = zeros(Nfft, Nofdm);
    
    X(88:88+Nf-1, :) = s(1:Nf, :);

    X = X + flipud(conj(X));

    % Saco piloto
    piloto = X(88:88+Nf-1,1);
    
    y_ofdm = ifft(X, Nfft, 'symmetric');
    
    % Añadir prefijo cíclico
    y_ofdm_cp = [y_ofdm(end-L_cp+1:end,:); y_ofdm];
    out = reshape(y_ofdm_cp, [], 1);
end
