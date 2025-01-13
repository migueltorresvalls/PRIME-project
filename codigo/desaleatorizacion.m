function bits_out = desaleatorizacion(bits_in,vector_aleatorizacion)
bits_out = xor(bits_in, vector_aleatorizacion);
end