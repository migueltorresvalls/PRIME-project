function rx_bits = demod_todo(M, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,PrefijoCiclico, Lcp)
if PrefijoCiclico == 0
    y_ofdm = demod_ofdm(x_ruido,Nfft,Nf,Nofdm);
else 
    y_ofdm = demod_ofdm_cp(x_ruido,Nfft,Nf,Nofdm,Lcp);
end
if M == 2
    dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
    rx_bits_aleatorio = dbpsk_demod(y_ofdm);
elseif M == 4
    dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
    rx_bits_aleatorio = dqpsk_demod(y_ofdm);
elseif M == 8
    d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
    rx_bits_aleatorio = d8psk_demod(y_ofdm);
end

rx_bits = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion);