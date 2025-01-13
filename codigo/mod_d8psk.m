function x_d8psk = mod_d8psk(tx_bits)
d8psk_mod = comm.DPSKModulator(8,0,'BitInput',true);

x_d8psk = d8psk_mod(tx_bits);
end