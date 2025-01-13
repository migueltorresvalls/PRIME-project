function [x_ofdm, piloto] = mod_todo(M, tx_bits, Nfft, Nofdm, Nf, PrefijoCiclico, Lcp)
if M == 2
    x_mod = mod_dbpsk(tx_bits);
elseif M == 4
    x_mod = mod_dqpsk(tx_bits);
elseif M == 8
    x_mod = mod_d8psk(tx_bits);
end
if PrefijoCiclico == 0
    [x_ofdm, piloto] = mod_ofdm(Nfft,Nofdm,Nf,x_mod);
else 
    [x_ofdm, piloto] = mod_ofdm_cp(Nfft,Nofdm,Nf,Lcp,x_mod);
end