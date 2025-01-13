function out = demod_ofdm(x_rx,Nfft,Nf,Nofdm)
x_rx = reshape(x_rx,Nfft,Nofdm);
x = fft(x_rx,Nfft);
x_mod = x(88:88+Nf-1,:);
out = reshape(x_mod,[],1);

end