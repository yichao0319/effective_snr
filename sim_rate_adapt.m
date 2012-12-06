%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%   input_sym_file: file that contains received symbols
%     format: <index> <magnitude> <phase> <real> <image>
%
% Output:
%
%
% Example:
%   sim_rate_adapt('RXDATA.s96-2.pdat')
%

function sim_rate_adapt(input_sym_file, method_ber2prr)
    

    %% ----------------------------------
    % constants
    PREDICTION_ORACLE = 'ORACLE';
    PREDICTION_EWMA = 'EWMA';

    GET_BER_ACTUAL_SNR = 'ACTUAL_SNR';
    GET_BER_GUESSED_EVM = 'GUESSED_EVM';
    GET_BER_ACTUAL_EVM = 'ACTUAL_EVM';
    GET_BER_AVERAGE_SNR = 'AVERAGE_SNR';

    PACKET_BER_PREAMBLE = 'PREAMBLE';
    PACKET_BER_ENTIRE = 'ENTIRE';
    PACKET_BER_ACTUAL = 'ACTUAL_BER';

    MODULATION_BPSK = 'BPSK';
    MODULATION_QPSK = 'QPSK';
    MODULATION_16QAM = '16QAM';
    MODULATION_64QAM = '64QAM';

    SNR2BER_FORMULA = 'FORMULA';
    SNR2BER_THRESHOLD = 'THRESHOLD';

    FRAME_ERR_FROM_BER = 'FROM_BER';
    FRAME_ERR_GILBERT = 'GILBERT';

    FEC_NO = 'NO_FEC';
    FEC_RS = 'RS';
    FEC_DIVERSITY_GAIN = 'DIVERSITY_GAIN';

    REMOVE_OUTLIER = 0;
    REMOVE_RAWOFDM_EMPTY_PERIOD = 1;
    READ_SNRDATA = 0;
    DEBUG_EXCEL_CHECK = 0;
    DEBUG_EVALUATE_FORMULA = 0;

    
    
    %% ----------------------------------
    % global variables
    input_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/rawofdm.4modulation/';
    tx_data_file = 'tx_syms_plane.dat';

    figure_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/figures_sim/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';

    num_subcarriers = 48;
    num_ofdm_symbol_per_pkt = 96;
    num_ofdm_symbol_valid = 95;
    num_ofdm_symbol_BPSK = 24;
    num_ofdm_symbol_QPSK = 24;
    num_ofdm_symbol_16QAM = 24;
    num_ofdm_symbol_64QAM = 23;

    % snr_shift = 0;
    % if strcmp(method_snr2ber, SNR2BER_FORMULA)
    %     snr_shift = 1.5;
    % end
    % if strcmp(method_FEC, FEC_DIVERSITY_GAIN)
    %     snr_shift = snr_shift + 3;
    % end
    num_preamble_ofdm_sym = 2;

    max_snr = 100;
    min_snr = -50;
    gran_snr = 0.5;
    max_evm = 100;
    min_evm = 0;
    gran_evm = 0.05;

    num_schemes = 10;
    num_modulation = 4;

    throughput = zeros(num_schemes, 1);
    wrong_selection = zeros(num_schemes, 1);
    throughput_table = [1 2 4 6];
    fid = fopen([output_dir input_sym_file '.' method_ber2prr '.out'], 'w');

    selected_modulations = ones(num_schemes, 1);


    %% ----------------------------------
    % initinalization
    % [qpsk_table, num_bit_per_sym] = mod_table(MODULATION_QPSK);
    % num_bits_per_frame = num_subcarriers * num_ofdm_symbol_per_pkt * num_bit_per_sym;
    % code_rate = code_rate_k / code_rate_n;
    % fprintf('code rate = %f\n', code_rate);
    fprintf(fid, 'pkt, actual BER,,,, mod, pred mod, effSNR BER (formula),,,, effSNR PRR,,,, mod, pred mod, effSNR BER (threshold),,,, effSNR PRR,,,, mod, pred mod, entire SNR BER (formula),,,, entire SNR PRR,,,, mod, pred mod, entire SNR BER (threshold),,,, entire SNR PRR,,,, mod, pred mod\n');
    fprintf(fid, ', BPSK, QPSK, 16QAM, 64QAM, mod, pred mod, ');
    for scheme_i = 1:(num_schemes-1)/2
        fprintf(fid, 'BPSK, QPSK, 16QAM, 64QAM, BPSK, QPSK, 16QAM, 64QAM, mod, pred mod, ');
    end
    fprintf(fid, '\n');
    % fprintf(fid, 'pkt, actual BER,,,, modulation, effSNR BER (formula),,,, effSNR prr (formula),,,, modulation; effSNR BER (threshold),,,, effSNR prr (threshold),,,, modulation; entire BER (formula),,,, entire prr (formula),,,, modulation; entire BER (threshold),,,, entire prr (threshold),,,, modulation;\npkt_i: BPSK, QPSK, 16QAM, 64QAM, selected_ind; BPSK, QPSK, 16QAM, 64QAM, BPSK, QPSK, 16QAM, 64QAM, selected_ind; BPSK, QPSK, 16QAM, 64QAM, BPSK, QPSK, 16QAM, 64QAM, selected_ind; BPSK, QPSK, 16QAM, 64QAM, BPSK, QPSK, 16QAM, 64QAM, selected_ind; BPSK, QPSK, 16QAM, 64QAM, BPSK, QPSK, 16QAM, 64QAM, selected_ind\n')

    
    %% ----------------------------------
    % main

    %% ----------------------------------
    % load tx data
    %   format: <real> <image>
    tx_raw = load([input_dir tx_data_file]);
    tx_syms = tx_raw(:, 1) + tx_raw(:, 2) * i;
    if length(tx_syms) ~= num_subcarriers * num_ofdm_symbol_per_pkt
        error('wrong number of tx file');
    end
    [tx_syms_BPSK, tx_syms_QPSK, tx_syms_16QAM, tx_syms_64QAM] = partition_one_pkt(tx_syms, num_subcarriers*num_ofdm_symbol_BPSK, num_subcarriers*num_ofdm_symbol_QPSK, num_subcarriers*num_ofdm_symbol_16QAM, num_subcarriers*num_ofdm_symbol_64QAM);
    [tx_slice_syms_BPSK, guessed_BPSK_evm] = slice('BPSK', tx_syms_BPSK);
    [tx_slice_syms_QPSK, guessed_QPSK_evm] = slice('QPSK', tx_syms_QPSK);
    [tx_slice_syms_16QAM, guessed_16QAM_evm] = slice('16QAM', tx_syms_16QAM);
    [tx_slice_syms_64QAM, guessed_64QAM_evm] = slice('64QAM', tx_syms_64QAM);
    tx_bits_BPSK = reshape(symbol2bit('BPSK', tx_slice_syms_BPSK)', [], 1);
    tx_bits_QPSK = reshape(symbol2bit('QPSK', tx_slice_syms_QPSK)', [], 1);
    tx_bits_16QAM = reshape(symbol2bit('16QAM', tx_slice_syms_16QAM)', [], 1);
    tx_bits_64QAM = reshape(symbol2bit('64QAM', tx_slice_syms_64QAM)', [], 1);
    tx_bits = [tx_bits_BPSK; tx_bits_QPSK; tx_bits_16QAM; tx_bits_64QAM];

    % tx_syms_grid = reshape(tx_syms, num_subcarriers, num_ofdm_symbol_per_pkt);
    tx_syms_BPSK_grid = reshape(tx_slice_syms_BPSK, num_subcarriers, num_ofdm_symbol_BPSK);
    tx_syms_QPSK_grid = reshape(tx_slice_syms_QPSK, num_subcarriers, num_ofdm_symbol_QPSK);
    tx_syms_16QAM_grid = reshape(tx_slice_syms_16QAM, num_subcarriers, num_ofdm_symbol_16QAM);
    tx_syms_64QAM_grid = reshape(tx_slice_syms_64QAM, num_subcarriers, num_ofdm_symbol_64QAM);
    tx_bits_BPSK_grid = symbol2bit('BPSK', tx_syms_BPSK_grid);
    tx_bits_QPSK_grid = symbol2bit('QPSK', tx_syms_QPSK_grid);
    tx_bits_16QAM_grid = symbol2bit('16QAM', tx_syms_16QAM_grid);
    tx_bits_64QAM_grid = symbol2bit('64QAM', tx_syms_64QAM_grid);



    %% ----------------------------------
    % load rx data
    %   format: <real> <image>
    rx_raw = load([input_dir input_sym_file]);
    rx_syms = rx_raw(:, 1) + rx_raw(:, 2) * i;
    num_pkts = floor(length(rx_syms) / num_subcarriers / num_ofdm_symbol_per_pkt);
    fprintf('num pkt: %d\n', num_pkts);
    % if num_pkts ~= 499
    %     error('wrong number of input pkt');
    % end



    for pkt_i = 1:num_pkts
        pkt_std_ind = (pkt_i - 1) * num_subcarriers * num_ofdm_symbol_per_pkt + 1;
        pkt_end_ind = pkt_i * num_subcarriers * num_ofdm_symbol_per_pkt;
        this_rx_syms = rx_syms(pkt_std_ind:pkt_end_ind, 1);
        if length(this_rx_syms) ~= num_subcarriers * num_ofdm_symbol_per_pkt
            error('wrong rx pkt size');
        end

        [rx_syms_BPSK, rx_syms_QPSK, rx_syms_16QAM, rx_syms_64QAM] = partition_one_pkt(this_rx_syms, num_subcarriers*num_ofdm_symbol_BPSK, num_subcarriers*num_ofdm_symbol_QPSK, num_subcarriers*num_ofdm_symbol_16QAM, num_subcarriers*num_ofdm_symbol_64QAM);
        [rx_slice_syms_BPSK, guessed_BPSK_evm] = slice('BPSK', rx_syms_BPSK);
        [rx_slice_syms_QPSK, guessed_QPSK_evm] = slice('QPSK', rx_syms_QPSK);
        [rx_slice_syms_16QAM, guessed_16QAM_evm] = slice('16QAM', rx_syms_16QAM);
        [rx_slice_syms_64QAM, guessed_64QAM_evm] = slice('64QAM', rx_syms_64QAM);

        rx_bits_BPSK = reshape(symbol2bit('BPSK', rx_slice_syms_BPSK)', [], 1);
        rx_bits_QPSK = reshape(symbol2bit('QPSK', rx_slice_syms_QPSK)', [], 1);
        rx_bits_16QAM = reshape(symbol2bit('16QAM', rx_slice_syms_16QAM)', [], 1);
        rx_bits_64QAM = reshape(symbol2bit('64QAM', rx_slice_syms_64QAM)', [], 1);
        rx_bits = [rx_bits_BPSK; rx_bits_QPSK; rx_bits_16QAM; rx_bits_64QAM];

        rx_syms_BPSK_grid = reshape(rx_syms_BPSK, num_subcarriers, num_ofdm_symbol_BPSK);
        rx_syms_QPSK_grid = reshape(rx_syms_QPSK, num_subcarriers, num_ofdm_symbol_QPSK);
        rx_syms_16QAM_grid = reshape(rx_syms_16QAM, num_subcarriers, num_ofdm_symbol_16QAM);
        rx_syms_64QAM_grid = reshape(rx_syms_64QAM, num_subcarriers, num_ofdm_symbol_64QAM);

        rx_slice_syms_BPSK_grid = reshape(rx_slice_syms_BPSK, num_subcarriers, num_ofdm_symbol_BPSK);
        rx_slice_syms_QPSK_grid = reshape(rx_slice_syms_QPSK, num_subcarriers, num_ofdm_symbol_QPSK);
        rx_slice_syms_16QAM_grid = reshape(rx_slice_syms_16QAM, num_subcarriers, num_ofdm_symbol_16QAM);
        rx_slice_syms_64QAM_grid = reshape(rx_slice_syms_64QAM, num_subcarriers, num_ofdm_symbol_64QAM);
        rx_bits_BPSK_grid = symbol2bit('BPSK', rx_slice_syms_BPSK_grid);
        rx_bits_QPSK_grid = symbol2bit('QPSK', rx_slice_syms_QPSK_grid);
        rx_bits_16QAM_grid = symbol2bit('16QAM', rx_slice_syms_16QAM_grid);
        rx_bits_64QAM_grid = symbol2bit('64QAM', rx_slice_syms_64QAM_grid);


        %% -------------------------
        % actual BER
        actual_ber = get_bit_error_rate(rx_bits', tx_bits');
        % fprintf('pkt %d: actual BER = %f\n', pkt_i, actual_ber);

        actual_ber_BPSK = get_bit_error_rate(rx_bits_BPSK', tx_bits_BPSK');
        actual_ber_QPSK = get_bit_error_rate(rx_bits_QPSK', tx_bits_QPSK');
        actual_ber_16QAM = get_bit_error_rate(rx_bits_16QAM', tx_bits_16QAM');
        actual_ber_64QAM = get_bit_error_rate(rx_bits_64QAM', tx_bits_64QAM');
        actual_bers = [actual_ber_BPSK, actual_ber_QPSK, actual_ber_16QAM, actual_ber_64QAM];

        % oracle throughput
        scheme_index = 1;
        min_ber = min(actual_bers);
        min_ber_ind = find(actual_bers == min_ber);
        selected_ind = min_ber_ind(end);
        fprintf('pkt %d: oracle=(%d, %f), ', pkt_i, selected_ind, min_ber);
        fprintf(fid, '%d: %f, %f, %f, %f, %d; ', pkt_i, actual_ber_BPSK, actual_ber_QPSK, actual_ber_16QAM, actual_ber_64QAM, selected_ind);
        % fprintf('-- a) %d: %f, %f, %f, %f, %d; \n', pkt_i, actual_ber_BPSK, actual_ber_QPSK, actual_ber_16QAM, actual_ber_64QAM, selected_ind);
        if min_ber == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end

        % oracle prediction throughput
        scheme_index = scheme_index + 1;
        new_selected_modulation = selected_ind;
        selected_ind = selected_modulations(scheme_index);
        fprintf('predicted oracle=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%d; ', selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end
        selected_modulations(scheme_index) = new_selected_modulation;




        %% -------------------------
        % actual SNR
        actual_snr_BPSK = calculate_actual_SNR(tx_syms_BPSK_grid, rx_syms_BPSK_grid);
        actual_snr_QPSK = calculate_actual_SNR(tx_syms_QPSK_grid, rx_syms_QPSK_grid);
        actual_snr_16QAM = calculate_actual_SNR(tx_syms_16QAM_grid, rx_syms_16QAM_grid);
        actual_snr_64QAM = calculate_actual_SNR(tx_syms_64QAM_grid, rx_syms_64QAM_grid);


        %% -------------------------
        % BER from actual SNR
        ber_per_sym_from_snr_grid_formula_BPSK = SNR2BER('BPSK', actual_snr_BPSK, SNR2BER_FORMULA);
        ber_per_sym_from_snr_grid_formula_QPSK = SNR2BER('QPSK', actual_snr_QPSK, SNR2BER_FORMULA);
        ber_per_sym_from_snr_grid_formula_16QAM = SNR2BER('16QAM', actual_snr_16QAM, SNR2BER_FORMULA);
        ber_per_sym_from_snr_grid_formula_64QAM = SNR2BER('64QAM', actual_snr_64QAM, SNR2BER_FORMULA);

        ber_per_sym_from_snr_grid_threshold_BPSK = SNR2BER('BPSK', actual_snr_BPSK, SNR2BER_THRESHOLD);
        ber_per_sym_from_snr_grid_threshold_QPSK = SNR2BER('QPSK', actual_snr_QPSK, SNR2BER_THRESHOLD);
        ber_per_sym_from_snr_grid_threshold_16QAM = SNR2BER('16QAM', actual_snr_16QAM, SNR2BER_THRESHOLD);
        ber_per_sym_from_snr_grid_threshold_64QAM = SNR2BER('64QAM', actual_snr_64QAM, SNR2BER_THRESHOLD);


        %% -------------------------
        % BER for effective SNR
        ber_eff_snr_formula_BPSK = mean(ber_per_sym_from_snr_grid_formula_BPSK(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_formula_QPSK = mean(ber_per_sym_from_snr_grid_formula_QPSK(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_formula_16QAM = mean(ber_per_sym_from_snr_grid_formula_16QAM(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_formula_64QAM = mean(ber_per_sym_from_snr_grid_formula_64QAM(:, 1:num_preamble_ofdm_sym), 2);

        ber_eff_snr_threshold_BPSK = mean(ber_per_sym_from_snr_grid_threshold_BPSK(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_threshold_QPSK = mean(ber_per_sym_from_snr_grid_threshold_QPSK(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_threshold_16QAM = mean(ber_per_sym_from_snr_grid_threshold_16QAM(:, 1:num_preamble_ofdm_sym), 2);
        ber_eff_snr_threshold_64QAM = mean(ber_per_sym_from_snr_grid_threshold_64QAM(:, 1:num_preamble_ofdm_sym), 2);
        if length(ber_eff_snr_formula_BPSK) ~= num_subcarriers
            size(ber_eff_snr_formula_BPSK)
            error('wrong ber length');
        end

        % packet reception rate (PRR) for effective SNR -- without FEC
        prr_eff_snr_formula_BPSK = ber2prr(method_ber2prr, 'BPSK', mean(ber_eff_snr_formula_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_formula_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_formula_QPSK = ber2prr(method_ber2prr, 'QPSK', mean(ber_eff_snr_formula_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_formula_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_formula_16QAM = ber2prr(method_ber2prr, '16QAM', mean(ber_eff_snr_formula_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_formula_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_formula_64QAM = ber2prr(method_ber2prr, '64QAM', mean(ber_eff_snr_formula_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_formula_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt);
        if length(prr_eff_snr_formula_BPSK) ~= 1
            size(prr_eff_snr_formula_BPSK)
            error('wrong ber length');
        end

        prr_eff_snr_threshold_BPSK = ber2prr(method_ber2prr, 'BPSK', mean(ber_eff_snr_threshold_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_threshold_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_threshold_QPSK = ber2prr(method_ber2prr, 'QPSK', mean(ber_eff_snr_threshold_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_threshold_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_threshold_16QAM = ber2prr(method_ber2prr, '16QAM', mean(ber_eff_snr_threshold_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_threshold_16QAM),num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_eff_snr_threshold_64QAM = ber2prr(method_ber2prr, '64QAM', mean(ber_eff_snr_threshold_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_eff_snr_threshold_64QAM),num_subcarriers*num_ofdm_symbol_per_pkt);

        % effective SNR throughput
        scheme_index = scheme_index + 1;
        if prr_eff_snr_formula_64QAM >= 0.9
            selected_ind = 4;
        elseif prr_eff_snr_formula_16QAM >= 0.9
            selected_ind = 3;
        elseif prr_eff_snr_formula_QPSK >= 0.9
            selected_ind = 2;
        else
            selected_ind = 1;
        end
        fprintf('effSNR1=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%f, %f, %f, %f; %f, %f, %f, %f, %d; ', mean(ber_eff_snr_formula_BPSK), mean(ber_eff_snr_formula_QPSK), mean(ber_eff_snr_formula_16QAM), mean(ber_eff_snr_formula_64QAM), prr_eff_snr_formula_BPSK, prr_eff_snr_formula_QPSK, prr_eff_snr_formula_16QAM, prr_eff_snr_formula_64QAM, selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end

        % effective SNR throughput with prediction
        scheme_index = scheme_index + 1;
        new_selected_modulation = selected_ind;
        selected_ind = selected_modulations(scheme_index);
        fprintf('predicted effSNR1=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%d; ', selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end
        selected_modulations(scheme_index) = new_selected_modulation;


        % effective SNR throughput (threshold)
        scheme_index = scheme_index + 1;
        if prr_eff_snr_threshold_64QAM >= 0.9
            selected_ind = 4;
        elseif prr_eff_snr_threshold_16QAM >= 0.9
            selected_ind = 3;
        elseif prr_eff_snr_threshold_QPSK >= 0.9
            selected_ind = 2;
        else
            selected_ind = 1;
        end
        fprintf('effSNR2=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%f, %f, %f, %f; %f, %f, %f, %f, %d; ', mean(ber_eff_snr_threshold_BPSK), mean(ber_eff_snr_threshold_QPSK), mean(ber_eff_snr_threshold_16QAM), mean(ber_eff_snr_threshold_64QAM), prr_eff_snr_threshold_BPSK, prr_eff_snr_threshold_QPSK, prr_eff_snr_threshold_16QAM, prr_eff_snr_threshold_64QAM, selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end

        % effective SNR throughput with prediction (threshold)
        scheme_index = scheme_index + 1;
        new_selected_modulation = selected_ind;
        selected_ind = selected_modulations(scheme_index);
        fprintf('predicted effSNR2=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%d; ', selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end
        selected_modulations(scheme_index) = new_selected_modulation;


        %% -------------------------
        % BER for entire SNR
        ber_entire_snr_formula_BPSK = mean(ber_per_sym_from_snr_grid_formula_BPSK, 2);
        ber_entire_snr_formula_QPSK = mean(ber_per_sym_from_snr_grid_formula_QPSK, 2);
        ber_entire_snr_formula_16QAM = mean(ber_per_sym_from_snr_grid_formula_16QAM, 2);
        ber_entire_snr_formula_64QAM = mean(ber_per_sym_from_snr_grid_formula_64QAM, 2);

        ber_entire_snr_threshold_BPSK = mean(ber_per_sym_from_snr_grid_threshold_BPSK, 2);
        ber_entire_snr_threshold_QPSK = mean(ber_per_sym_from_snr_grid_threshold_QPSK, 2);
        ber_entire_snr_threshold_16QAM = mean(ber_per_sym_from_snr_grid_threshold_16QAM, 2);
        ber_entire_snr_threshold_64QAM = mean(ber_per_sym_from_snr_grid_threshold_64QAM, 2);
        if length(ber_entire_snr_formula_BPSK) ~= num_subcarriers
            error('wrong ber length');
        end

        % packet reception rate (PRR) for using entire SNR -- without FEC
        prr_entire_snr_formula_BPSK = ber2prr(method_ber2prr, 'BPSK', mean(ber_entire_snr_formula_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_formula_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_formula_QPSK = ber2prr(method_ber2prr, 'QPSK', mean(ber_entire_snr_formula_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_formula_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_formula_16QAM = ber2prr(method_ber2prr, '16QAM', mean(ber_entire_snr_formula_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_formula_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_formula_64QAM = ber2prr(method_ber2prr, '64QAM', mean(ber_entire_snr_formula_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_formula_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt);

        prr_entire_snr_threshold_BPSK = ber2prr(method_ber2prr, 'BPSK', mean(ber_entire_snr_threshold_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_threshold_BPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_threshold_QPSK = ber2prr(method_ber2prr, 'QPSK', mean(ber_entire_snr_threshold_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_threshold_QPSK), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_threshold_16QAM = ber2prr(method_ber2prr, '16QAM', mean(ber_entire_snr_threshold_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_threshold_16QAM), num_subcarriers*num_ofdm_symbol_per_pkt);
        prr_entire_snr_threshold_64QAM = ber2prr(method_ber2prr, '64QAM', mean(ber_entire_snr_threshold_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt); %power(1-mean(ber_entire_snr_threshold_64QAM), num_subcarriers*num_ofdm_symbol_per_pkt);

        % entire SNR throughput
        scheme_index = scheme_index + 1;
        if prr_entire_snr_formula_64QAM >= 0.9
            selected_ind = 4;
        elseif prr_entire_snr_formula_16QAM >= 0.9
            selected_ind = 3;
        elseif prr_entire_snr_formula_QPSK >= 0.9
            selected_ind = 2;
        else
            selected_ind = 1;
        end
        fprintf('entireSNR1=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%f, %f, %f, %f; %f, %f, %f, %f, %d; ', mean(ber_entire_snr_formula_BPSK), mean(ber_entire_snr_formula_QPSK), mean(ber_entire_snr_formula_16QAM), mean(ber_entire_snr_formula_64QAM), prr_entire_snr_formula_BPSK, prr_entire_snr_formula_QPSK, prr_entire_snr_formula_16QAM, prr_entire_snr_formula_64QAM, selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end

        % entire SNR throughput with prediction
        scheme_index = scheme_index + 1;
        new_selected_modulation = selected_ind;
        selected_ind = selected_modulations(scheme_index);
        fprintf('predicted entireSNR1=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%d; ', selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end
        selected_modulations(scheme_index) = new_selected_modulation;


        % entire SNR throughput (threshold)
        scheme_index = scheme_index + 1;
        if prr_entire_snr_threshold_64QAM >= 0.9
            selected_ind = 4;
        elseif prr_entire_snr_threshold_16QAM >= 0.9
            selected_ind = 3;
        elseif prr_entire_snr_threshold_QPSK >= 0.9
            selected_ind = 2;
        else
            selected_ind = 1;
        end
        fprintf('entireSNR2=(%d, %f)\n', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%f, %f, %f, %f; %f, %f, %f, %f, %d; ', mean(ber_entire_snr_threshold_BPSK), mean(ber_entire_snr_threshold_QPSK), mean(ber_entire_snr_threshold_16QAM), mean(ber_entire_snr_threshold_64QAM), prr_entire_snr_threshold_BPSK, prr_entire_snr_threshold_QPSK, prr_entire_snr_threshold_16QAM, prr_entire_snr_threshold_64QAM, selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end

        % entire SNR throughput with prediction (threshold)
        scheme_index = scheme_index + 1;
        new_selected_modulation = selected_ind;
        selected_ind = selected_modulations(scheme_index);
        fprintf('predicted entireSNR2=(%d, %f), ', selected_ind, actual_bers(selected_ind));
        fprintf(fid, '%d; ', selected_ind);
        if actual_bers(selected_ind) == 0
            throughput(scheme_index) = throughput(scheme_index) + throughput_table(selected_ind);
        else
            wrong_selection(scheme_index) = wrong_selection(scheme_index) + 1;
        end
        selected_modulations(scheme_index) = new_selected_modulation;
        

        if scheme_index ~= num_schemes
            error('wrong number of schemes');
        end
        fprintf(fid, '\n');
    end


    wrong_selection
    throughput


    fclose(fid);

end


%% ----------------------------------
% functions

%% get_snr_of_packet_preamble: function description
function [estimated_snr] = get_snr_of_packet_preamble(snr, num_preamble_ofdm_sym)
    estimated_snr = mean(snr(:, 1:num_preamble_ofdm_sym), 2);
end

%% get_snr_of_packet_entire: function description
function [estimated_snr] = get_snr_of_packet_entire(snr)
   estimated_snr = mean(snr, 2);
end


%% get_ber_of_packet_preamble: function description
function [estimated_ber] = get_ber_of_packet_preamble(ber, method, num_preamble_ofdm_sym)
    if strcmp(method, 'mean')
        estimated_ber = mean(ber(:, 1:num_preamble_ofdm_sym), 2);
    elseif strcmp(method, 'calculation')
        estimated_ber = 1-prod(1-ber(:, 1:num_preamble_ofdm_sym), 2);
    else
        error('Unimplemented method');
    end 
end

%% get_ber_of_packet_entire: function description
function [estimated_ber] = get_ber_of_packet_entire(ber, method)
    if strcmp(method, 'mean')
        estimated_ber = mean(ber, 2);
    elseif strcmp(method, 'calculation')
        estimated_ber = 1-prod(1-ber, 2);
    else
        error('Unimplemented method');
    end 
end


%% prediction_error:real_ts, estimate description
function [err] = prediction_error(real_ts, estimated_ts)
   err = zeros(size(real_ts));
   real_len = size(real_ts, 2);

   err = abs(real_ts - estimated_ts(:, 1:real_len) );
end


%% get_symbol_error_rate: function description
function [ber] = get_symbol_error_rate(rx_grid, tx_grid)
    if size(rx_grid) ~= size(tx_grid) 
        disp('wrong pkt size');
    end

    ber = zeros(size(rx_grid, 1), 1);
    % rx_grid(1, 1)
    % tx_grid(1, 1)
    % abs(rx_grid(1, 1) - tx_grid(1, 1))

    for sc_i = 1:size(rx_grid, 1)
        wrong_ind = find(abs(rx_grid(sc_i, :) - tx_grid(sc_i, :)) > 0.2);
        num_err = length(wrong_ind);

        ber(sc_i, 1) = num_err / size(rx_grid, 2);
    end
end


%% get_bit_error_rate: function description
function [ber] = get_bit_error_rate(rx_bit_grid, tx_bit_grid)
    if size(rx_bit_grid) ~= size(tx_bit_grid) 
        disp('wrong pkt size');
    end

    [num_subcarriers, num_bits_per_sc] = size(rx_bit_grid);
    ber = zeros(num_subcarriers, 1);

    for sc_i = 1:num_subcarriers
        wrong_ind = find(rx_bit_grid(sc_i, :) ~= tx_bit_grid(sc_i, :) );
        num_err = length(wrong_ind);

        ber(sc_i, 1) = num_err / num_bits_per_sc;
    end
end


%% calculate_actual_EVM: function description
function [evm] = calculate_actual_EVM(rx_grid, tx_grid)
    num_pkts = floor(size(rx_grid, 2) / size(tx_grid, 2));
    remain_syms = size(rx_grid, 2) - num_pkts * size(tx_grid, 2);
    tx_grid_rep = [repmat(tx_grid, 1, num_pkts) tx_grid(:, 1:remain_syms)];

    evm = abs(rx_grid - tx_grid_rep) / 1;
end


%% calculate_actual_SNR: function description
function [snr] = calculate_actual_SNR(rx_grid, tx_grid)
    num_pkts = floor(size(rx_grid, 2) / size(tx_grid, 2));
    remain_syms = size(rx_grid, 2) - num_pkts * size(tx_grid, 2);
    tx_grid_rep = [repmat(tx_grid, 1, num_pkts) tx_grid(:, 1:remain_syms)];

    snr = power(abs(tx_grid_rep), 2) ./ power(abs(rx_grid - tx_grid_rep), 2);
    snr = 10 * log(snr) / log(10);
end


%% cal_actual_SNR2BER: function description
function [snr_err_rate, snr_cnt] = cal_actual_SNR2BER(snr, rx_bit_grid, tx_bit_grid, min_snr, max_snr, gran)
    if(size(snr) ~= [size(rx_bit_grid, 1) size(rx_bit_grid, 2)*2]) 
        disp('wrong number of snr and rx_grid');
        return
    end


    num_pkts = floor(size(rx_bit_grid, 2) / size(tx_bit_grid, 2));
    remain_bits = size(rx_bit_grid, 2) - num_pkts * size(tx_bit_grid, 2);
    tx_bit_grid_rep = [repmat(tx_bit_grid, 1, num_pkts) tx_bit_grid(:, 1:remain_bits)];


    snr_err_cnt = zeros((max_snr - min_snr)/gran + 1, 2);
    snr_err_rate = zeros((max_snr - min_snr)/gran + 1, 1);


    [row col] = size(snr);
    for row_i = 1:row
        for col_i = 1:col
            this_snr = floor((snr(row_i, col_i) - min_snr) / gran ) + 1;
            if(rx_bit_grid(row_i, 2*(col_i-1)+1) == tx_bit_grid_rep(row_i, 2*(col_i-1)+1))
                snr_err_cnt(this_snr, 2) = snr_err_cnt(this_snr, 2) + 1;
            else
                snr_err_cnt(this_snr, 2) = snr_err_cnt(this_snr, 2) + 1;
                snr_err_cnt(this_snr, 1) = snr_err_cnt(this_snr, 1) + 1;
            end

            if(rx_bit_grid(row_i, 2*col_i) == tx_bit_grid_rep(row_i, 2*col_i))
                snr_err_cnt(this_snr, 2) = snr_err_cnt(this_snr, 2) + 1;
            else
                snr_err_cnt(this_snr, 2) = snr_err_cnt(this_snr, 2) + 1;
                snr_err_cnt(this_snr, 1) = snr_err_cnt(this_snr, 1) + 1;
            end
        end
    end

    % disp('actual SNR 2 BER');
    for snr_i = 1:(max_snr - min_snr) / gran + 1
        % fprintf('%d, %d\n', snr_err_cnt(snr_i, 1), snr_err_cnt(snr_i, 2));
        if(snr_err_cnt(snr_i, 1) > 0)
            snr_err_rate(snr_i, 1) = snr_err_cnt(snr_i, 1) / snr_err_cnt(snr_i, 2);
        end
    end

    snr_cnt = snr_err_cnt(:, 2);
end


%% cal_actual_EVM2BER: function description
function [evm_err_rate, evm_cnt] = cal_actual_EVM2BER(evm, rx_bit_grid, tx_bit_grid, min_evm, max_evm, gran)
    if(size(evm) ~= [size(rx_bit_grid, 1) size(rx_bit_grid, 2)*2]) 
        disp('wrong number of evm and rx_grid');
        return
    end


    num_pkts = floor(size(rx_bit_grid, 2) / size(tx_bit_grid, 2));
    remain_bits = size(rx_bit_grid, 2) - num_pkts * size(tx_bit_grid, 2);
    tx_bit_grid_rep = [repmat(tx_bit_grid, 1, num_pkts) tx_bit_grid(:, 1:remain_bits)];


    evm_err_cnt = zeros((max_evm - min_evm)/gran + 1, 2);
    evm_err_rate = zeros((max_evm - min_evm)/gran + 1, 1);


    [row col] = size(evm);
    for row_i = 1:row
        for col_i = 1:col
            this_evm = floor((evm(row_i, col_i) - min_evm) / gran ) + 1;
            if(rx_bit_grid(row_i, 2*(col_i-1)+1) == tx_bit_grid_rep(row_i, 2*(col_i-1)+1))
                evm_err_cnt(this_evm, 2) = evm_err_cnt(this_evm, 2) + 1;
            else
                evm_err_cnt(this_evm, 2) = evm_err_cnt(this_evm, 2) + 1;
                evm_err_cnt(this_evm, 1) = evm_err_cnt(this_evm, 1) + 1;
            end

            if(rx_bit_grid(row_i, 2*col_i) == tx_bit_grid_rep(row_i, 2*col_i))
                evm_err_cnt(this_evm, 2) = evm_err_cnt(this_evm, 2) + 1;
            else
                evm_err_cnt(this_evm, 2) = evm_err_cnt(this_evm, 2) + 1;
                evm_err_cnt(this_evm, 1) = evm_err_cnt(this_evm, 1) + 1;
            end
        end
    end

    for evm_i = 1:(max_evm - min_evm) / gran + 1
        if(evm_err_cnt(evm_i, 1) > 0)
            evm_err_rate(evm_i, 1) = evm_err_cnt(evm_i, 1) / evm_err_cnt(evm_i, 2);
        end
    end

    evm_cnt = evm_err_cnt(:, 2);
end



%% partition_one_pkt: function description
function [syms_BPSK, syms_QPSK, syms_16QAM, syms_64QAM] = partition_one_pkt(pkt_syms, a, b, c, d)
    syms_BPSK  = pkt_syms(1:a, 1);
    syms_QPSK  = pkt_syms(1+a:a+b, 1);
    syms_16QAM = pkt_syms(1+a+b:a+b+c, 1);
    syms_64QAM = pkt_syms(1+a+b+c:a+b+c+d, 1);
end


%% ber2prr: function description
function [prr] = ber2prr(method, modulation, ber, num_bits)
    if strcmp(method, 'FORMULA')
        prr = power(1-ber, num_bits);
    elseif strcmp(method, 'THRESHOLD')
        if strcmp(modulation, 'BPSK')
            prr = power(1-ber, num_bits);
        elseif strcmp(modulation, 'QPSK')
            if ber > 0.001
                prr = 0;
            else
                prr = 1;
            end
        elseif strcmp(modulation, '16QAM')
            if ber > 0.000207
                prr = 0;
            else
                prr = 1;
            end
        elseif strcmp(modulation, '64QAM')
            if ber > 0.0026631
                prr = 0;
            else
                prr = 1;
            end
        end
    end
end
    






