function x_dqpsk = mod_dqpsk(tx_bits)
dqpsk_mod = comm.DPSKModulator(4,0,'BitInput',true);

x_dqpsk = dqpsk_mod(tx_bits);
end