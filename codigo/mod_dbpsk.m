function x_dbpsk = mod_dbpsk(tx_bits)
dbpsk_mod = comm.DPSKModulator(2,0,'BitInput',true);

x_dbpsk = dbpsk_mod(tx_bits);
end