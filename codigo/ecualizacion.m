function out = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp)
x_rx = reshape(x_ruido, Nfft+Lcp, Nofdm);

% Eliminar prefijo cíclico
x_rx = x_rx(Lcp+1:end, :); 

% Hago DFT
x = fft(x_rx,Nfft);

% Estimación del canal en funcnión del primer símbolo recibido
% (correspondiente al piloto).
w_i = piloto./x(88:88+Nf-1,1);

% Concateno coeficientes del canal
w = zeros(Nfft,Nofdm);
w(88:88+Nf-1,:) = repmat(w_i,1,Nofdm);

% Ecualizo
x_eq = x.*w;

% Vector fila
x_mod = x_eq(88:88+Nf-1,:);
out = reshape(x_mod,[],1);
end