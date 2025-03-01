close all; clear; clc;
%% 參數設定
M = 4;                  % 4-PAM
k = log2(M);            % 每個符號所含的 bit 數 (4-PAM -> 2 bits/symbol)
EsN0_dB = 0:12;       % 模擬的 Es/N0 (dB) 範圍

trellis = poly2trellis(3, [5 7]);
k_tb = 16;
%% 產生隨機位元
txBits = load('dataIn.asv', '-ascii');
txBits = txBits(1:16000);
txBits = transpose(txBits);

%% 不同長度組的conv_data

traceback = [16 32 128]; 
conv_msg = cell(1, length(traceback)); 

for io = 1:length(traceback)
    conv_msg{io} = []; % 存入空陣
end

for conv_rest_length = 1:length(traceback)
    count = 1 ;
    j_st = length(txBits) / traceback(conv_rest_length);
    for cj = 1:j_st
        conv_code = conv_hardware_213(txBits(count:count+traceback(conv_rest_length)-1));
        conv_msg{conv_rest_length} = [conv_msg{conv_rest_length}, conv_code];   
        count = count + traceback(conv_rest_length);
    end
end

symIdxTx_Uncoded_data = bi2de(reshape(txBits, k, []).','left-msb');
txSymbols_Uncoded = pammod(symIdxTx_Uncoded_data, M, 0, 'gray');  

conv_data = conv_hardware_213(txBits);
symIdxTx_data = bi2de(reshape(conv_data, k, []).','left-msb');
txSymbols_hardware_behavior = pammod(symIdxTx_data, M, 0, 'gray');  

%% viterbi_v5.1 ber
symIdxTx_data = bi2de(reshape(conv_msg{1}, k, []).','left-msb');
txSymbols_16 = pammod(symIdxTx_data, M, 0, 'gray');  

symIdxTx_data = bi2de(reshape(conv_msg{2}, k, []).','left-msb');
txSymbols_32 = pammod(symIdxTx_data, M, 0, 'gray');  

symIdxTx_data = bi2de(reshape(conv_msg{3}, k, []).','left-msb');
txSymbols_128 = pammod(symIdxTx_data, M, 0, 'gray');  

% 每個SNR 測試10次平均
num_trials = 10; % 每個 SNR 測試 10 次

for ii = 1:length(EsN0_dB)
    ber_16_sum = 0;
    ber_errors_16_sum = 0;
    for trial = 1:num_trials
    rxSymbolBlock_16 = awgn(txSymbols_16,ii,'measured');
    % 解調 (帶灰度解調)
    symIdxRx_16 = pamdemod(rxSymbolBlock_16, M, 0, 'gray');
    recovered_symIdxRx_16 = de2bi(symIdxRx_16, k, 'left-msb');
    recovered_bits_16 = reshape(recovered_symIdxRx_16.', 1, []);
    
    decoded_hardware_behavior_V51_16 = viterbi_hardwareV5_1(recovered_bits_16,16);
    
    [bit_errors_16, ber_16] = biterr(txBits, decoded_hardware_behavior_V51_16);
    ber_16_sum = ber_16_sum + ber_16;
    ber_errors_16_sum = ber_errors_16_sum + bit_errors_16;
    end
    BER_16(ii) = ber_16_sum / num_trials;
    BER_errors_16(ii) = bit_errors_16 / num_trials;
end

for ii = 1:length(EsN0_dB)
    ber_32_sum = 0;
    ber_errors_32_sum = 0;
    for trial = 1:num_trials
        rxSymbolBlock_32 = awgn(txSymbols_32,ii,'measured');
        % 解調 (帶灰度解調)
        symIdxRx_32 = pamdemod(rxSymbolBlock_32, M, 0, 'gray');
        recovered_symIdxRx_32 = de2bi(symIdxRx_32, k, 'left-msb');
        recovered_bits_32 = reshape(recovered_symIdxRx_32.', 1, []);
        
        decoded_hardware_behavior_V51_32 = viterbi_hardwareV5_1(recovered_bits_32,32);
        
        [bit_errors_32, ber_32] = biterr(txBits, decoded_hardware_behavior_V51_32);
        ber_32_sum = ber_32_sum + ber_32;
        ber_errors_32_sum = ber_errors_32_sum + bit_errors_32;
    end
    BER_32(ii) = ber_32_sum / num_trials;
    BER_errors_32(ii) = bit_errors_32 / num_trials;
end

for ii = 1:length(EsN0_dB)
    ber_128_sum = 0;
    ber_errors_128_sum = 0;
    for trial = 1:num_trials    
    rxSymbolBlock_128 = awgn(txSymbols_128,ii,'measured');
    % 解調 (帶灰度解調)
    symIdxRx_128 = pamdemod(rxSymbolBlock_128, M, 0, 'gray');
    recovered_symIdxRx_128 = de2bi(symIdxRx_128, k, 'left-msb');
    recovered_bits_128 = reshape(recovered_symIdxRx_128.', 1, []);
    
    decoded_hardware_behavior_V51_128 = viterbi_hardwareV5_1(recovered_bits_128,128);
    
    [bit_errors_128, ber_128] = biterr(txBits, decoded_hardware_behavior_V51_128);
    ber_128_sum = ber_128_sum + ber_128;
    ber_errors_128_sum = ber_errors_128_sum + bit_errors_128;
    end
    BER_128(ii) = ber_128_sum / num_trials;
    BER_errors_128(ii) = ber_errors_128_sum / num_trials;
end

%% viterbi_V5
for ii = 1:length(EsN0_dB)
    % 初始化累加變數
    matlab_ber_sum = 0;
    matlab_errors_sum = 0;
    hardware_ber_sum = 0;
    hardware_errors_sum = 0;
    for trial = 1:num_trials
        rxSymbolBlock = awgn(txSymbols_hardware_behavior,ii,'measured');
        % 解調 (帶灰度解調)
        symIdxRx = pamdemod(rxSymbolBlock, M, 0, 'gray');
        recovered_symIdxRx_data = de2bi(symIdxRx, k, 'left-msb');
        recovered_bits = reshape(recovered_symIdxRx_data.', 1, []);

        % MATLAB Viterbi 解碼
        decoded_matlab = vitdec(recovered_bits, trellis, k_tb, 'term', 'hard');
        decoded_hardware_behavior = viterbi_hardwareV5(recovered_bits,k_tb);

        % MATLAB Viterbi BER 計算
        [ber_MATLAB_errors, ber_MATLAB_trial] = biterr(txBits, decoded_matlab);
        matlab_ber_sum = matlab_ber_sum + ber_MATLAB_trial;
        matlab_errors_sum = matlab_errors_sum + ber_MATLAB_errors;

        % MATLAB Hardware Behavior Viterbi BER 計算
        [ber_hardware_errors, ber_hardware_trial] = biterr(txBits, decoded_hardware_behavior);
        hardware_ber_sum = hardware_ber_sum + ber_hardware_trial;
        hardware_errors_sum = hardware_errors_sum + ber_hardware_errors;
    end

    % 計算 10 次的平均 BER 和錯誤數
    BER_MATLAB(ii) = matlab_ber_sum / num_trials;
    BER_MATLAB_ERROR(ii) = matlab_errors_sum / num_trials;
    BER_HARDWARE_BEHAVIOR(ii)  = hardware_ber_sum / num_trials;
    BER_HARDWARE_BEHAVIOR_ERROR(ii) = hardware_errors_sum / num_trials;
end

%% 計算Uncoded ber
for ii = 1:length(EsN0_dB)
    ber_Uncoded_sum = 0;
    ber_Uncoded_error_sum = 0;
    for trial = 1:num_trials 
    rxSymbolBlock_Uncoded = awgn(txSymbols_Uncoded,ii,'measured');
    % 解調 (帶灰度解調)
    symIdxRx_Uncoded = pamdemod(rxSymbolBlock_Uncoded, M, 0, 'gray');

    [bit_Uncoded_errors, ber_Uncoded] = biterr(symIdxTx_Uncoded_data, symIdxRx_Uncoded);
    
    ber_Uncoded_sum = ber_Uncoded_sum + ber_Uncoded;
    ber_Uncoded_error_sum = ber_Uncoded_error_sum + bit_Uncoded_errors;
    end
    BER_Uncoded(ii) = ber_Uncoded_sum / num_trials;
    BER_Uncoded_errors(ii) = ber_Uncoded_error_sum / num_trials;
end


%% **繪圖 (BER / SER)** 軟+硬
figure;
hold on; 
markers = {'-*', '-d', '-^', '-v', '-<', '->', '-p', '-h'}; % 8 個標記
styles = {'r--s','k--x', 'm--^'};  % 增加样式

semilogy(EsN0_dB, BER_Uncoded, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Uncoded BER');
semilogy(EsN0_dB, BER_16, styles{1}, 'LineWidth', 1.5, 'DisplayName', 'viterbi_v5.1_16');
semilogy(EsN0_dB, BER_32, styles{2}, 'LineWidth', 1.5, 'DisplayName', 'viterbi_v5.1_32');
semilogy(EsN0_dB, BER_128, styles{3}, 'LineWidth', 1.5, 'DisplayName', 'viterbi_v5.1_128');
semilogy(EsN0_dB, BER_MATLAB, markers{1}, 'LineWidth', 1.5, 'DisplayName', 'MATLAB_Viterbi');
semilogy(EsN0_dB, BER_HARDWARE_BEHAVIOR, markers{2}, 'LineWidth', 1.5, 'DisplayName', 'Hardware behavior ViterbiV51 TB=16');

grid on;
set(gca, 'YScale', 'log');
set(gca, 'XTick', EsN0_dB); 
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show; 
title('BER Performance Comparison (Viterbi & PAM4)');
hold off; 

%---------VITERBI_213_HARDWARE_V5.1-----------
function [decoded_msg] = viterbi_hardwareV5_1(conv_code,D)
    survivors = cell(4,1);       % 儲存各狀態的存活路徑
    new_survivors = cell(4,1);   % 暫存新生成的存活路徑
    path_metrics = [0;3;3;3];    % 各狀態的路徑度量值
    decoded_msg = [];            % 儲存逐步輸出的解碼結果
    
    for step = 1:length(conv_code)/2
        idx = 2*step - 1;
        received_bits = conv_code(idx:idx+1);
        new_metrics = inf(4,1);    % 初始化新度量為無窮大
        new_survivors = cell(4,1); % 重置暫存器
        
        for current_state = 0:3
            for input_bit = 0:1
    
                next_state = viterbi_next_state(current_state, input_bit);
                output_dec = viterbi_outputs(current_state, input_bit);
                expected_bits = de2bi(output_dec, 2, 'left-msb');
                
                hamming_dist = sum(received_bits ~= expected_bits);
                
                candidate_metric = path_metrics(current_state+1) + hamming_dist;
                
                if candidate_metric < new_metrics(next_state+1)
                    new_metrics(next_state+1) = candidate_metric;
                    new_survivors{next_state+1} = [survivors{current_state+1}, input_bit];
                end
            end
        end
        
        path_metrics = new_metrics;
        survivors = new_survivors;
    
        row_lengths = cellfun(@length, survivors);
        if any(row_lengths >= D)
            [~, best_state] = min(path_metrics);
            best_path = survivors{best_state};
            decoded_msg = [decoded_msg, best_path];
            survivors = cell(4,1);
            path_metrics = [0;3;3;3];
        end
    end
    
    [~, final_state] = min(path_metrics);
    remaining_bits = survivors{final_state};
    decoded_msg = [decoded_msg, remaining_bits];
end

%---------VITERBI_213_HARDWARE_V5-----------
function [decoded_msg] = viterbi_hardwareV5(conv_code,D)
    survivors = cell(4,1);       % 儲存各狀態的存活路徑
    new_survivors = cell(4,1);   % 暫存新生成的存活路徑
    path_metrics = [0;3;3;3];    % 各狀態的路徑度量值
    decoded_msg = [];            % 儲存逐步輸出的解碼結果
    
    for step = 1:length(conv_code)/2
        idx = 2*step - 1;
        received_bits = conv_code(idx:idx+1);
        new_metrics = inf(4,1);    % 初始化新度量為無窮大
        new_survivors = cell(4,1); % 重置暫存器
        
        for current_state = 0:3
            for input_bit = 0:1
    
                next_state = viterbi_next_state(current_state, input_bit);
                output_dec = viterbi_outputs(current_state, input_bit);
                expected_bits = de2bi(output_dec, 2, 'left-msb');
                
                hamming_dist = sum(received_bits ~= expected_bits);
                
                candidate_metric = path_metrics(current_state+1) + hamming_dist;
                
                if candidate_metric < new_metrics(next_state+1)
                    new_metrics(next_state+1) = candidate_metric;
                    new_survivors{next_state+1} = [survivors{current_state+1}, input_bit];
                end
            end
        end
        
        path_metrics = new_metrics;
        survivors = new_survivors;
    
        row_lengths = cellfun(@length, survivors);
        if any(row_lengths >= D)
            [~, best_state] = min(path_metrics);
            best_path = survivors{best_state};
            decoded_msg = [decoded_msg, best_path];
            survivors = cell(4,1);
        end
    end
    
    [~, final_state] = min(path_metrics);
    remaining_bits = survivors{final_state};
    decoded_msg = [decoded_msg, remaining_bits];
end

%---------VITERBI_HARDWARE_TABLE-----------
function nextState = viterbi_next_state(currentState,inptBits)
    switch currentState
        case 0 
            if(inptBits == 0)
                nextState = 0;
            else
                nextState = 2;
            end
         case 1 
            if(inptBits == 0)
                nextState = 0;
            else
                nextState = 2;
            end
         case 2 
            if(inptBits == 0)
                nextState = 1;
            else
                nextState = 3;
            end
         case 3 
            if(inptBits == 0)
                nextState = 1;
            else
                nextState = 3;
            end
    end
end

function outputs = viterbi_outputs(currentState,inptBits)
    switch currentState
        case 0 
            if(inptBits == 0)
                outputs = 0;
            else
                outputs = 3;
            end
         case 1 
            if(inptBits == 0)
                outputs = 3;
            else
                outputs = 0;
            end
         case 2 
            if(inptBits == 0)
                outputs = 1;
            else
                outputs = 2;
            end
         case 3 
            if(inptBits == 0)
                outputs = 2;
            else
                outputs = 1;
            end
    end
end

%---------hardware_conv213_function-----------
function codeword = conv_hardware_213(msg_source)
    s1 = 0;
    s2 = 0;
    bit_string_length  = length(msg_source);
    codeword = zeros(1, bit_string_length * 2);
    for i = 1:bit_string_length
        u0 = xor(msg_source(i), s2);
        u1 = xor(xor(msg_source(i), s1), s2);
        s2 = s1;
        s1 = msg_source(i);
        codeword(2*i-1) = u0;
        codeword(2*i) = u1;
    end
end

%--------awgn_measured----------
function [noisy_signal,signal_power] = awgn_measured(input_signal, snr_dB)
    % 測量輸入信號的平均功率
    signal_power = sum(abs(input_signal(:)).^2)/numel(input_signal);
    
    % 將 SNR 從 dB 轉換為線性比例
    snr_linear = 10^(snr_dB/10);
    
    % 計算所需的噪聲功率
    noise_power = signal_power / snr_linear;
    
    % 產生高斯白噪聲，均值為 0，變異數為 noise_power
    noise = sqrt(noise_power) * randn(size(input_signal));
    
    % 將噪聲添加到原始信號
    noisy_signal = input_signal + noise;
end
