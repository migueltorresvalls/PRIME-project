clear all
clc
%% Modulaciones
% Número de símbolos ofdm. Igual para todas las modulaciones
M = 44; 

%d8psk
N_d8psk = 288*M;
% Fijamos número de tramas d8psk a 2
N_tramas_d8psk = 2;
N_bits_totales = N_d8psk*N_tramas_d8psk;

%dqpsk
N_dqpsk = 192*M;
N_tramas_dqpsk = ceil(N_bits_totales/N_dqpsk);

% dbpsk
N_dbpsk = 96*M;
N_tramas_dbpsk = ceil(N_bits_totales/N_dbpsk);

%% Cálculo BER
% Vectores auxiliares para bucle for
N_tramas = [N_tramas_dbpsk, N_tramas_dqpsk, N_tramas_d8psk];
L_tramas = [N_dbpsk, N_dqpsk, N_d8psk];

Nfft = 512; % Número de puntos fft. Número de subportadores disponibles
Nofdm = M; % Número de símbolos ofdm
Nf = 96; % Número de subportadoras

% Vector de SNRs
SNR_dB = -5:2:40;
offset_SNR = 10*log10(Nfft/(2*Nf));

% Vectores de BER calculada y BER teórica
BER = zeros(3, length(SNR_dB));
BER_teor = zeros(3, length(SNR_dB));

% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(2, tx_aleatorio, Nfft, Nofdm, Nf,0,0);
                
                x_ruido = awgn(x_ofdm, SNR_dB(snr)-offset_SNR,'measured');
                rx_bits = demod_todo(2, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
                
            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(4, tx_aleatorio, Nfft, Nofdm, Nf,0,0);
                
                x_ruido = awgn(x_ofdm, SNR_dB(snr)-offset_SNR,'measured');
                rx_bits = demod_todo(4, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(8, tx_aleatorio, Nfft, Nofdm, Nf,0,0);
                
                x_ruido = awgn(x_ofdm, SNR_dB(snr)-offset_SNR,'measured');
                rx_bits = demod_todo(8, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        
        % Cálculo BER teórica para todas las tramas de la modulación
        if i == 1
            BER_t = DBPSK_BER(SNR_dB(snr));
        elseif i == 2
            BER_t = DQPSK_BER(SNR_dB(snr));
        elseif i == 3
            BER_t = D8PSK_BER(SNR_dB(snr));
        end
        BER_t(BER_t<1e-5)=NaN;
        BER_teor(i,snr) = BER_t;
    end
end

%% Representación BER vs SNR
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR canal ideal')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR canal ideal')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR canal ideal')


%% Expresión algebraica de la señal recibida incluyendo únicamente el efecto del canal
% La señal recibida se calcula como:
% r(t) = s(t) * h(t)
% Donde h(t) es la respuesta al impulso del canal.

% Respuesta al impulso del canal
h = [-0.1 0.3 -0.5 0.7 -0.9 0.7 -0.5 0.3 -0.1];

%% Representación gráfica de la función de transferencia del canal en frecuencia, en dB
Nfft = 512; % Número de puntos para la FFT
H = fft(h, Nfft); % Transformada en frecuencia del canal
H_dB = 20 * log10(abs(H)); % Magnitud en dB
fs = 250e3; % A partir de número de muestras FFT y periodo de símbolo OFDM
f = ((0:length(H)-1)/length(H)-0.5)*fs;
%f = fs*(Nfft-1);

% Gráfica de la función de transferencia del canal
figure
plot(f, H_dB)
grid
title('Función de transferencia del canal en frecuencia')
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')

%% Inclusión del efecto del canal en el código y obtención de curvas BER sin prefijo cíclico ni ecualización
% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(2, tx_aleatorio, Nfft, Nofdm, Nf,0,0);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Demodulación
                rx_bits = demod_todo(2, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
                
            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(4, tx_aleatorio, Nfft, Nofdm, Nf,0,0);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Demodulación
                rx_bits = demod_todo(4, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(8, tx_aleatorio, Nfft, Nofdm, Nf,0,0);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Demodulación
                rx_bits = demod_todo(8, x_ruido, Nfft, Nofdm, Nf, vector_aleatorizacion,0,0);
            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
    end
end

%% Representación BER vs SNR con canal
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal')

%% Inclusión de prefijo cíclico y ecualizador en recepción
% Según especificaciones PRIME
Lcp = 48;

% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(2, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_bits_aleatorio = dbpsk_demod(x_eq);
                rx_bits = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion);  

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(4, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_bits_aleatorio = dqpsk_demod(x_eq);
                rx_bits = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion);      
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                [x_ofdm, piloto] = mod_todo(8, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_bits_aleatorio = d8psk_demod(x_eq);
                rx_bits = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion);   

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end


%% Representación BER vs SNR con canal y ecualizacion
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ecualizado')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ecualizado')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ecualizado')

%% Inclusión de entrelazado en transmisor y receptor
% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);

                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(2, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_entrelazado = dbpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion);  

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);

                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(4, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_entrelazado = dqpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion);  
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(8, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_entrelazado = d8psk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion);  

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end

%% Representación BER vs SNR con canal sin distorisón, ecualizacion y entrelazado
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y entrelazado')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y entrelazado')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y entrelazado')

%% Modificación de los parámetros de simulación para soportar códigos convolucionales
% Memoria del codificador
constraint_length = 7; 

% Bits vaciado
flush_bits = zeros(constraint_length,1);

%d8psk
N_d8psk = 144*M;
% Fijamos número de tramas d8psk a 2
N_tramas_d8psk = 2;
N_bits_totales = N_d8psk*N_tramas_d8psk;

%dqpsk
N_dqpsk = 96*M;
N_tramas_dqpsk = ceil(N_bits_totales/N_dqpsk);

% dbpsk
N_dbpsk = 48*M;
N_tramas_dbpsk = ceil(N_bits_totales/N_dbpsk);

% Vector dimensiones
N_tramas = [N_tramas_dbpsk, N_tramas_dqpsk, N_tramas_d8psk];
L_tramas = [N_dbpsk, N_dqpsk, N_d8psk];

%% Inclusión de FEC a canal sin ruido
% 171 y 133 son números en base octal correspondiente a 1111001 y 1011011
% respectivamente. 
trellis = poly2trellis(constraint_length,[171 133]);
% Según sistema PRIME
traceback_length = 42;

% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);

                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(2, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_bits_aleatorio = dbpsk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(4, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_bits_aleatorio = dqpsk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(8, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_canal,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_bits_aleatorio = d8psk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end

%% Representación BER vs SNR con canal sin ruido, ecualizacion y FEC
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y FEC')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y FEC')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal ideal y FEC')

%% Inclusión de FEC a canal con ruido
% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);
                
                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(2, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_bits_aleatorio = dbpsk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);
                
                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(4, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_bits_aleatorio = dqpsk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                [x_ofdm, piloto] = mod_todo(8, tx_aleatorio, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_bits_aleatorio = d8psk_demod(x_eq);
                rx_codec = desaleatorizacion(rx_bits_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end

%% Representación BER vs SNR con canal, ecualizacion y FEC
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal, ecualizador y FEC')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal, ecualizador y FEC')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR con canal, ecualizador y FEC')

%% Todos los bloques del sistema: Codificador, Aleatorización, Entrelazado, Ecualización
constraint_length = 7; 
trellis = poly2trellis(constraint_length,[171 133]);
traceback_length = 42;

% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);
                
                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);

                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(2, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_entrelazado = dbpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_codec = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);
                
                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(4, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_entrelazado = dqpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_codec = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Codifico bits
                tx_codec = convenc([tx_bits; flush_bits], trellis);
                tx_codec_flush = tx_codec(end-2*length(flush_bits)+1:end);
                tx_codec = tx_codec(1:end-2*length(flush_bits));

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_codec);
                
                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(8, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_entrelazado = d8psk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_codec = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 

                % Decodifico
                rx_bits = vitdec([rx_codec; tx_codec_flush],trellis,traceback_length,'trunc', 'hard');
                rx_bits = rx_bits(1:end-length(flush_bits));

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end

%% Representación BER vs SNR sistema completo
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo')

%% Todos los bloques del sistema menos FEC: Aleatorización, Entrelazado, Ecualización
% Ajuste parámetros sin FEC
%d8psk
N_d8psk = 288*M;
% Fijamos número de tramas d8psk a 2
N_tramas_d8psk = 2;
N_bits_totales = N_d8psk*N_tramas_d8psk;

%dqpsk
N_dqpsk = 192*M;
N_tramas_dqpsk = ceil(N_bits_totales/N_dqpsk);

% dbpsk
N_dbpsk = 96*M;
N_tramas_dbpsk = ceil(N_bits_totales/N_dbpsk);

% Vector dimensiones
N_tramas = [N_tramas_dbpsk, N_tramas_dqpsk, N_tramas_d8psk];
L_tramas = [N_dbpsk, N_dqpsk, N_d8psk];

% Bucle que itera para todas las SNR
for snr=1:length(SNR_dB)
    % Bucle que itera para todas las modulaciones
    for i=1:3
        l_trama = L_tramas(i);
        n_tramas = N_tramas(i);
    
        bits_recibidos = zeros(n_tramas,l_trama);
        bits_transmitidos = zeros(n_tramas,l_trama);
        % Bucle que itera para todas las tramas
        for j=1:n_tramas
            % Si modulación de subportadora es DBPSK
            if i == 1
                tx_bits = randi([0,1],l_trama,1);
                
                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);

                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(2, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                % Piloto del 87 al 87 + 96
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dbpsk_demod = comm.DPSKDemodulator(2,0,'BitOutput',true);
                rx_entrelazado = dbpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 

            % Si la modulación de subportadora es DQPSK
            elseif i == 2
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(4, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                dqpsk_demod = comm.DPSKDemodulator(4,0,'BitOutput',true);
                rx_entrelazado = dqpsk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 
                
            % Si la modulación de subportadora es D8PSK
            elseif i == 3
                tx_bits = randi([0,1],l_trama,1);

                % Aleatorización
                [tx_aleatorio, vector_aleatorizacion] = aleatorizacion(tx_bits);
                
                % Entrelazado
                tx_entrelazado = entrelazado(tx_aleatorio,2*l_trama/M);
                
                [x_ofdm, piloto] = mod_todo(8, tx_entrelazado, Nfft, Nofdm, Nf,1,Lcp);

                % Aplicación del canal
                x_canal = filter(h, 1, x_ofdm); 
                x_ruido = awgn(x_canal, SNR_dB(snr)-offset_SNR,'measured');

                % Ecualizacion
                x_eq = ecualizacion(piloto,Nfft,Nofdm,Nf,x_ruido,Lcp);

                % Demodulación
                d8psk_demod = comm.DPSKDemodulator(8,0,'BitOutput',true);
                rx_entrelazado = d8psk_demod(x_eq);
                rx_aleatorio = desentrelazado(rx_entrelazado,2*l_trama/M);
                rx_bits = desaleatorizacion(rx_aleatorio, vector_aleatorizacion); 

            end
            % Guardo tramas recibidas en una matriz
            bits_recibidos(j,:) = rx_bits;
            bits_transmitidos(j,:) = tx_bits;
        end
        % Concateno filas de la matriz
        bits_recibidos_concatenados = reshape(bits_recibidos,1,N_bits_totales);
        bits_transmitidos_concatenados = reshape(bits_transmitidos,1,N_bits_totales);

        % Cálculo BER para todas las tramas de la modulación
        BER(i,snr) = sum(abs(bits_recibidos_concatenados-bits_transmitidos_concatenados))/N_bits_totales;
        %BER(BER<1e-5) = NaN;
    end
end

%% Representación BER vs SNR sistema completo menos FEC
figure
semilogy(SNR_dB, BER(1,:))
hold on
semilogy(SNR_dB, BER_teor(1,:))
grid
legend('DBPSK','DBPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo menos FEC')

figure
semilogy(SNR_dB, BER(2,:))
hold on
semilogy(SNR_dB, BER_teor(2,:))
grid
legend('DQPSK','DQPSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo menos FEC')

figure
semilogy(SNR_dB, BER(3,:))
hold on
semilogy(SNR_dB, BER_teor(3,:))
grid
legend('D8PSK','D8PSK_{teorica}')
ylabel('BER (Bit Error Rate)')
xlabel('SNR (Signal to Noise Relation)')
title('Curvas BER vs SNR sistema completo menos FEC')
