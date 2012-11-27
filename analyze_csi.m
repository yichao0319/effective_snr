%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%   input_sym_file: file that contains received symbols
%     format: <index> <magnitude> <phase> <real> <image>
%   method_prediction: how to predict BER of next packet
%     [ORACLE | EWMA]
%   method_get_ber: the way to get BER of each bit
%     [ACTUAL_SNR | GUESSED_EVM | ACTUAL_EVM | AVERAGE_SNR]
%   method_pkt_ber: use which part of BER of the packet to get effective SNR
%     [PREAMBLE | ENTIRE]
%   method_snr2ber: the function used to convert SNR to BER
%     [FORMULA | THRESHOLD]
%   num_ofdm_symbol_per_pkt:
%   modulation:
%     [BPSK | QPSK | 16QAM | 64QAM]
%   alpha: the parameter for EWMA
%
% Output:
%
%
% Example:
%   analyze_csi('RXDATA.s96-2.pdat', 'ORACLE', 'ACTUAL_SNR', 'PREAMBLE', 90, 'QPSK', 1)
%

function [avg_actual_ber, avg_predit_ber] = analyze_csi(input_sym_file, ...
                     method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, ...
                     num_ofdm_symbol_per_pkt, modulation, ...
                     alpha)
    

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

    REMOVE_OUTLIER = 0;
    REMOVE_RAWOFDM_EMPTY_PERIOD = 1;
    READ_SNRDATA = 0;
    DEBUG_EXCEL_CHECK = 0;
    DEBUG_EVALUATE_FORMULA = 1;

    
    
    %% ----------------------------------
    % global variables
    input_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/PARSEDDATA/';
    tx_data_file = 'tx-data-s240.pdat';

    figure_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/figures/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT/';

    num_subcarriers = 48;
    % num_ofdm_symbol_per_pkt = 90;
    % modulation = 'QPSK';
    snr_shift = 0;
    if strcmp(method_snr2ber, SNR2BER_FORMULA)
        snr_shift = 1.5;
    end
    num_preamble_ofdm_sym = 2;

    max_snr = 100;
    min_snr = -50;
    gran_snr = 0.5;
    max_evm = 100;
    min_evm = 0;
    gran_evm = 0.05;


    %% ----------------------------------
    % initinalization
    [qpsk_table, num_bit_per_sym] = mod_table(MODULATION_QPSK);


    
    %% ----------------------------------
    % main

    %% ----------------------------------
    % load tx data
    %   format: <magnitude> <phase> <real> <image>
    tx_data = load([input_dir tx_data_file]);
    tx_data_cpx = tx_data(:, 3) + tx_data(:, 4) * i;    
    tx_data_cpx = tx_data_cpx(1:num_subcarriers*num_ofdm_symbol_per_pkt);   % e.g. (4320 * 1)
    tx_grid = zeros(num_subcarriers, num_ofdm_symbol_per_pkt);
    tx_grid = reshape(tx_data_cpx, num_subcarriers, []);                    % e.g. (48 * 90)
    tx_bit_grid = symbol2bit(modulation, tx_grid);                          % e.g. (48 * 180)
    disp(sprintf('load tx data size = %d, %d', size(tx_data_cpx) ) );
    if size(tx_grid, 2) ~= num_ofdm_symbol_per_pkt
        error('wrong number of tx grid');
    end

    %% DEBUG
    %{
    tmp_ind = find(tx_data_cpx == 0);
    fprintf('number of 0 in tx data: %d\n', length(tmp_ind));
    %}


    %% ----------------------------------
    % load SNR reported from rawOFDM
    %   format: <index> <SNR>
    %{
    if READ_SNRDATA == 1
        raw_snr_data = load([input_dir input_snr_file]);
        snr_data = raw_snr_data(:, 2);
    end
    %}


    %% ----------------------------------
    % load rx symbols to rx_grid
    %   format: <index> <magnitude> <phase> <real> <image>
    data = load([input_dir input_sym_file]);
    data_cpx = data(:, 4) + data(:, 5) * i;                             % e.g. (4320x * 1)
    num_symbols = size(data, 1);
    disp(sprintf('load data size = %d, %d', size(data) ) );


    %% ----------------------------------
    % remove outlier
    if REMOVE_OUTLIER == 1
        outlier_index = find(abs(data_cpx > 5) );
        data_cpx(outlier_index) = 0.000001;
    end


    %% ----------------------------------
    % put rx symbols into grid
    rx_grid = zeros(num_subcarriers, num_symbols/num_subcarriers);
    rx_grid = reshape(data_cpx, num_subcarriers, []);                   % e.g. (48 * 90x)
    disp(sprintf('grid size: %d, %d', size(rx_grid) ));


    %% ----------------------------------
    % For some reason, rawOFDM generates 0+0i symbols in the end of packets.
    % The codes here recognize the initial 7 symbols of each packet and 
    % only take the specified number of OFDM symbols per packet
    if REMOVE_RAWOFDM_EMPTY_PERIOD == 1
        num_pkts = 1;
        pkt_start_ind_ary = [];
        sym_i = 1;
        while sym_i < num_symbols - num_subcarriers*num_ofdm_symbol_per_pkt
            while ~( ( real(data_cpx(sym_i  )) < 0 & imag(data_cpx(sym_i  )) < 0 ) & ...
                     ( real(data_cpx(sym_i+1)) < 0 & imag(data_cpx(sym_i+1)) > 0 ) & ...
                     ( real(data_cpx(sym_i+2)) > 0 & imag(data_cpx(sym_i+2)) < 0 ) & ...
                     ( real(data_cpx(sym_i+3)) < 0 & imag(data_cpx(sym_i+3)) > 0 ) & ...
                     ( real(data_cpx(sym_i+4)) < 0 & imag(data_cpx(sym_i+4)) < 0 ) & ...
                     ( real(data_cpx(sym_i+5)) > 0 & imag(data_cpx(sym_i+5)) > 0 ) & ...
                     ( real(data_cpx(sym_i+6)) > 0 & imag(data_cpx(sym_i+6)) > 0 ) )

                sym_i = sym_i + 1;
            end
            pkt_start_ind_ary = [pkt_start_ind_ary sym_i];
            sym_i = sym_i + num_subcarriers*num_ofdm_symbol_per_pkt;

            % if length(pkt_start_ind_ary) > num_pkts
            %     break;
            % end
        end
        % pkt_start_ind_ary
        % return


        rx_grid = [];
        for pkt_i = 1:length(pkt_start_ind_ary)-1
            start_ind = pkt_start_ind_ary(pkt_i);
            end_ind = start_ind + num_subcarriers * num_ofdm_symbol_per_pkt - 1;

            rx_grid = [rx_grid reshape(data_cpx(start_ind:end_ind), num_subcarriers, [])];  % e.g. (48 * 90x)
        end
        disp(sprintf('grid size: %d, %d', size(rx_grid) ));
    end
    
    num_pkts = floor(size(rx_grid, 2) / num_ofdm_symbol_per_pkt);
    fprintf('number of packets: %d\n', num_pkts);
    

    %% DEBUG
    %{
    tmp_ind = find(data_cpx == 0);
    fprintf('number of 0 in rx data: %d\n', length(tmp_ind));
    %}



    %% ----------------------------------
    % DEBUG
    %{
    tmp_tx = tx_data_cpx';
    tmp_tx_bit = symbol2bit(modulation, tmp_tx);
    size(tmp_tx_bit)
    tmp_rx = data_cpx(1:48*96)';
    tmp_rx_bit = symbol2bit(modulation, tmp_rx);
    size(tmp_rx_bit)
    tmp_err = get_bit_error_rate(tmp_tx_bit, tmp_rx_bit)
    wrong_ind = find(tmp_tx_bit ~= tmp_rx_bit);
    dlmwrite([output_dir input_sym_file '.tmp_tx.txt'], [real(tmp_tx)' imag(tmp_tx)']);
    dlmwrite([output_dir input_sym_file '.tmp_rx.txt'], [real(tmp_rx)' imag(tmp_rx)']);
    dlmwrite([output_dir input_sym_file '.tmp_tx_bit.txt'], reshape(tmp_tx_bit, 2, [])');
    dlmwrite([output_dir input_sym_file '.tmp_rx_bit.txt'], reshape(tmp_rx_bit, 2, [])');
    dlmwrite([output_dir input_sym_file '.tmp_err_ind.txt'], wrong_ind');
    dlmwrite([output_dir input_sym_file '.tmp_snr_data.txt'], snr_data(1:48*96)');
    tmp_snr2ber = SNR2BER(modulation, snr_data(1:48*96), method_snr2ber);
    dlmwrite([output_dir input_sym_file '.tmp_snr2ber.txt'], tmp_snr2ber');
    tmp_avg_snr = mean(snr_data(1:48*96));
    fprintf('avg snr: %.10f\n', tmp_avg_snr);
    tmp_avg_snr2ber = SNR2BER(modulation, tmp_avg_snr, method_snr2ber);
    fprintf('avg estimated ber: %.10f\n', tmp_avg_snr2ber);
    tmp_avg_snr2ber2 = mean(tmp_snr2ber);
    fprintf('avg estimated ber: %.10f\n', tmp_avg_snr2ber2);
    tmp_num_err1 = 48*96*tmp_avg_snr2ber;
    fprintf('number of estimated bit error: %.10f\n', tmp_num_err1);
    tmp_num_err2 = 48*96*tmp_avg_snr2ber2;
    fprintf('number of estimated bit error: %.10f\n', tmp_num_err2);
    tmp_overall_ber = 1;
    for tmp_i = 1:48*96
        tmp_overall_ber = tmp_overall_ber * (1-tmp_snr2ber(tmp_i));
    end
    tmp_overall_ber = 1-tmp_overall_ber
    48*96*tmp_overall_ber
    return
    %}


    %% ----------------------------------
    % DEBUG: 
    %   pick one packet, calculate SNR, BER, ... and 
    %   put them to excel file (excel_check.xlsx) for manual examination
    %
    if DEBUG_EXCEL_CHECK == 1
        tmp_pkt = 1;
        tmp_pkt_std_ind = (tmp_pkt-1) * num_ofdm_symbol_per_pkt + 1;
        tmp_pkt_end_ind = tmp_pkt * num_ofdm_symbol_per_pkt;
        tmp_tx_sym = reshape(tx_data_cpx, 1, []);
        tmp_tx_bit = symbol2bit(modulation, tmp_tx_sym);
        tmp_rx_sym = reshape(rx_grid(:, tmp_pkt_std_ind:tmp_pkt_end_ind), [], 1);
        tmp_rx_bit = symbol2bit(modulation, reshape(tmp_rx_sym, 1, []));
        dlmwrite([output_dir input_sym_file '.tmp_rx_sym.txt'], [real(tmp_rx_sym) imag(tmp_rx_sym)]);
        dlmwrite([output_dir input_sym_file '.tmp_rx_bit.txt'], reshape(tmp_rx_bit, 2, [])');
        tmp_actual_snr = calculate_actual_SNR(reshape(tmp_rx_sym, 1, []), reshape(tx_data_cpx, 1, [])) + snr_shift;
        dlmwrite([output_dir input_sym_file '.tmp_actual_snr.txt'], tmp_actual_snr');
        tmp_ber = SNR2BER(modulation, tmp_actual_snr, method_snr2ber);
        dlmwrite([output_dir input_sym_file '.tmp_ber.txt'], tmp_ber');
        tmp_cal_ber = get_bit_error_rate(tmp_rx_bit, tmp_tx_bit);
        dlmwrite([output_dir input_sym_file '.tmp_cal_ber.txt'], tmp_cal_ber');
        % return;
    end
    %



    %% ----------------------------------
    % 1. EVM and rx bits
    %   i) get EVM from closest constellation points
    [rx_slice_grid, guessed_evm] = slice(modulation, rx_grid);  % e.g. rx_slice_grid: (48 * 90x)
                                                                % e.g. guessed_evm: (48 * 90x)
    rx_bit_grid = symbol2bit(modulation, rx_slice_grid);        % e.g. rx_bit_grid: (48 * 180x)
    %   ii) get EVM by assuming we know the correct constellation points 
    % actual_evm = calculate_actual_EVM(rx_grid, tx_grid);              % e.g. actual_evm: (48 * 90x)
    

    %% ----------------------------------
    % 2. SNR and BER
    snr_per_sym = [];
    ber_per_sym = [];
    if strcmp(method_get_ber, GET_BER_ACTUAL_SNR)
        % knowing the tx and rx symbols, we can calculate the exact SNRs
        snr_per_sym = calculate_actual_SNR(rx_grid, tx_grid) + snr_shift;
        ber_per_sym = SNR2BER(modulation, snr_per_sym, method_snr2ber);

    elseif strcmp(method_get_ber, GET_BER_ACTUAL_EVM)
        % get EVM by assuming we know the correct constellation points 
        actual_evm = calculate_actual_EVM(rx_grid, tx_grid);              % e.g. actual_evm: (48 * 90x)
        % converted actual EVM to SNR 
        snr_per_sym = EVM2SNR(actual_evm) + snr_shift;
        ber_per_sym = SNR2BER(modulation, snr_per_sym, method_snr2ber);

    elseif strcmp(method_get_ber, GET_BER_GUESSED_EVM)
        % converted from EVM calculated from closest constellation points
        snr_per_sym = EVM2SNR(guessed_evm) + snr_shift;
        ber_per_sym = SNR2BER(modulation, snr_per_sym, method_snr2ber);

    elseif strcmp(method_get_ber, GET_BER_AVERAGE_SNR)
        snr_per_sym = calculate_actual_SNR(rx_grid, tx_grid) + snr_shift;

    else
        error('unknown method to get BER');
    end


    %% ----------------------------------
    % DEBUG:
    %   evaluate SNR to BER, EVM to BER curve
    if DEBUG_EVALUATE_FORMULA == 1
        actual_snr = calculate_actual_SNR(rx_grid, tx_grid) + snr_shift;
        actual_evm = calculate_actual_EVM(rx_grid, tx_grid);
        [snr_err_rate, snr_cnt] = cal_actual_SNR2BER(actual_snr, rx_bit_grid, tx_bit_grid, min_snr, max_snr, gran_snr);
        [evm_err_rate, evm_cnt] = cal_actual_EVM2BER(actual_evm, rx_bit_grid, tx_bit_grid, min_evm, max_evm, gran_evm);
        f9 = figure;
        plot(min_snr:gran_snr:max_snr, snr_err_rate, 'r*', ...
             min_snr:gran_snr:max_snr, snr_cnt/max(snr_cnt), 'g');
        hold on;
        x_values = min_snr:0.1:max_snr;
        y_values_qpsk = SNR2BER('QPSK', x_values, method_snr2ber);
        plot(x_values, y_values_qpsk);
        xlabel('SNR (dB)');
        ylabel('BER');
        legend('actual', ...
               ['count (/' int2str(max(snr_cnt)) ')'], ...
               'formula');
        axis([min_snr max_snr 0 1]);
        print(f9, '-dpsc', [figure_dir input_sym_file '.' method_snr2ber '.snr2ber.ps']);

        f10 = figure;
        plot(min_evm:gran_evm:max_evm, evm_err_rate, 'r*', ...
             min_evm:gran_evm:max_evm, evm_cnt/max(evm_cnt), 'g');
        hold on;
        x_values = min_evm:0.01:max_evm;
        y_values_qpsk = SNR2BER('QPSK', EVM2SNR(x_values), method_snr2ber);
        y_values_qpsk_shift = SNR2BER('QPSK', EVM2SNR(x_values+1), method_snr2ber);
        plot(x_values, y_values_qpsk);
        xlabel('EVM');
        ylabel('BER');
        legend('actual', ...
               ['count (/' int2str(max(evm_cnt)) ')'], ...
               'formula');
        % axis([-10 15 0 0.5]);
        print(f10, '-dpsc', [figure_dir input_sym_file '.' method_snr2ber '.evm2ber.ps']);
        % return;
    end
    



    %% ----------------------------------
    % 3. get packet BER
    actual_pkt_ber = zeros(num_subcarriers, num_pkts);          % e.g. actual_pkt_ber: (48 * x)
    estimated_pkt_ber = zeros(num_subcarriers, num_pkts);       % e.g. estimated_pkt_ber: (48 * x)

    for pkt_i = 1:num_pkts
        % actual packet BER
        pkt_str_bit_ind = (pkt_i-1) * num_bit_per_sym * num_ofdm_symbol_per_pkt + 1;
        pkt_end_bit_ind = pkt_i     * num_bit_per_sym * num_ofdm_symbol_per_pkt;
        this_rx_bit_grid = rx_bit_grid(:, pkt_str_bit_ind:pkt_end_bit_ind);
        actual_pkt_ber(:, pkt_i) = get_bit_error_rate(this_rx_bit_grid, tx_bit_grid);
        

        % estimated packet BER
        pkt_str_ofdm_ind = (pkt_i-1) * num_ofdm_symbol_per_pkt + 1;
        pkt_end_ofdm_ind = pkt_i     * num_ofdm_symbol_per_pkt;
        

        if strcmp(method_get_ber, GET_BER_AVERAGE_SNR)
            this_snr_grid = snr_per_sym(:, pkt_str_ofdm_ind:pkt_end_ofdm_ind);
            if strcmp(method_pkt_ber, PACKET_BER_PREAMBLE)
                this_pkt_snr = get_snr_of_packet_preamble(this_snr_grid, num_preamble_ofdm_sym);
                estimated_pkt_ber(:, pkt_i) = SNR2BER(modulation, this_pkt_snr, method_snr2ber);

            elseif strcmp(method_pkt_ber, PACKET_BER_ENTIRE)
                this_pkt_snr = get_snr_of_packet_entire(this_snr_grid);
                estimated_pkt_ber(:, pkt_i) = SNR2BER(modulation, this_pkt_snr, method_snr2ber);

            else
                error('unknown method to get packet BER');
            end            

        else
            this_ber_grid = ber_per_sym(:, pkt_str_ofdm_ind:pkt_end_ofdm_ind);
            if strcmp(method_pkt_ber, PACKET_BER_PREAMBLE)
                estimated_pkt_ber(:, pkt_i) = get_ber_of_packet_preamble(this_ber_grid, 'mean', num_preamble_ofdm_sym);
            elseif strcmp(method_pkt_ber, PACKET_BER_ENTIRE)
                estimated_pkt_ber(:, pkt_i) = get_ber_of_packet_entire(this_ber_grid, 'mean');
            else
                error('unknown method to get packet BER');
            end

        end
    end

        
    %% ----------------------------------
    % 4. prediction
    predicted_estimated_pkt_ber = [];
    if strcmp(method_prediction, PREDICTION_ORACLE)
        predicted_estimated_pkt_ber = estimated_pkt_ber;

    elseif strcmp(method_prediction, PREDICTION_EWMA)
        predicted_estimated_pkt_ber = ewma(estimated_pkt_ber, alpha);
        
        num_prediction = size(predicted_estimated_pkt_ber, 2);
        predicted_estimated_pkt_ber = predicted_estimated_pkt_ber(:, 1:num_prediction-1);
        if(size(predicted_estimated_pkt_ber) ~= size(actual_pkt_ber))
            error('wrong number of ewma prediction');
        end

    else
        error('unknown method for prediction');
    end
        
        
    %% ----------------------------------
    % 5. evaluation
    %   i) get error per pkt
    prediction_err_per_pkt = abs(mean(actual_pkt_ber(:, 2:end), 1) - mean(predicted_estimated_pkt_ber(:, 2:end), 1) );
    prediction_err_per_pkt_ratio = prediction_err_per_pkt ./ mean(actual_pkt_ber(:, 2:end), 1);
    dlmwrite([output_dir input_sym_file '.' method_prediction '.' method_get_ber '.' method_pkt_ber '.' method_snr2ber '.pkt.eval.txt'], [prediction_err_per_pkt' mean(actual_pkt_ber(:, 2:end), 1)' prediction_err_per_pkt_ratio']);
    % fprintf('error per pkt = \n');
    % fprintf('%f,', prediction_err_per_pkt_ratio);
    % fprintf('\n');
    
    f5a = figure;
    plot(1:num_pkts-1, mean(actual_pkt_ber(:, 2:end), 1), ...
         1:num_pkts-1, mean(predicted_estimated_pkt_ber(:, 2:end), 1), ...
         1:num_pkts-1, prediction_err_per_pkt);
    xlabel('pkt');
    ylabel('BER');
    legend('actual', 'predicted', 'error');
    print(f5a, '-dpsc', [figure_dir input_sym_file '.' method_prediction '.' method_get_ber '.' method_pkt_ber '.' method_snr2ber '.pkt.eval.ps']);
    

    %   ii) get error per subcarrier
    prediction_err_per_sc = abs(mean(actual_pkt_ber(:, 2:end), 2) - mean(predicted_estimated_pkt_ber(:, 2:end), 2) );
    prediction_err_per_sc_ratio = prediction_err_per_sc ./ mean(actual_pkt_ber(:, 2:end), 2);
    dlmwrite([output_dir input_sym_file '.' method_prediction '.' method_get_ber '.' method_pkt_ber '.' method_snr2ber '.sc.eval.txt'], [prediction_err_per_sc mean(actual_pkt_ber(:, 2:end), 2) prediction_err_per_sc_ratio]);
    % fprintf('error per subcarrier = \n');
    % fprintf('%f,', prediction_err_per_sc_ratio);
    % fprintf('\n');
    
    f5b = figure;
    plot(1:num_subcarriers, mean(actual_pkt_ber(:, 2:end), 2), ...
         1:num_subcarriers, mean(predicted_estimated_pkt_ber(:, 2:end), 2), ...
         1:num_subcarriers, prediction_err_per_sc);
    xlabel('subcarrier');
    ylabel('BER');
    legend('actual', 'predicted', 'error');
    print(f5b, '-dpsc', [figure_dir input_sym_file '.' method_prediction '.' method_get_ber '.' method_pkt_ber '.' method_snr2ber '.sc.eval.ps']);


    %   iii) get error as a whole
    avg_actual_ber = mean(mean(actual_pkt_ber(:, 2:end) ) );
    avg_predit_ber = mean(mean(predicted_estimated_pkt_ber(:, 2:end) ) );

    prediction_err_whole = abs(avg_actual_ber - avg_predit_ber);
    prediction_err_whole_ratio = prediction_err_whole / avg_actual_ber;
    dlmwrite([output_dir input_sym_file '.' method_prediction '.' method_get_ber '.' method_pkt_ber '.' method_snr2ber '.whole.eval.txt'], [avg_predit_ber avg_actual_ber prediction_err_whole prediction_err_whole_ratio]);
    % fprintf('error as a whole = |%f - %f| / %f = %f\n', avg_actual_ber, avg_predit_ber, avg_actual_ber, prediction_err_whole_ratio);
    



    % %   i) knowing the tx and rx symbols, we can calculate the exact SNRs
    % actual_snr = calculate_actual_SNR(rx_grid, tx_grid) + snr_shift;
    % %   ii) converted from EVM calculated from closest constellation points
    % snr_from_guessed_evm = EVM2SNR(guessed_evm) + snr_shift;
    % %   iii) converted from EVM knowing the correct constellation points
    % snr_from_actual_evm = EVM2SNR(actual_evm) + snr_shift;


    % %% ----------------------------------
    % % 3. BER per symbols (from SNR)
    % ber_from_actual_snr = SNR2BER(modulation, actual_snr, method_snr2ber);
    % ber_from_guessed_evm = SNR2BER(modulation, snr_from_guessed_evm, method_snr2ber);
    % ber_from_actual_evm = SNR2BER(modulation, snr_from_actual_evm, method_snr2ber);
    % disp(sprintf('ber_from_actual_snr  size: %d, %d', size(ber_from_actual_snr) ));
    % disp(sprintf('ber_from_guessed_evm size: %d, %d', size(ber_from_guessed_evm) ));
    % disp(sprintf('ber_from_actual_evm size: %d, %d', size(ber_from_actual_evm) ));


    % %% ----------------------------------
    % % 4. BER per packets
    % num_pkts = floor(size(rx_grid, 2) / num_ofdm_symbol_per_pkt);
    % fprintf('number of packets: %d\n', num_pkts);
    
    % % i) from tx and rx bit
    % pkt_ber_from_rxtx_bit = zeros(num_subcarriers, num_pkts);

    % % ii) from preamble
    % pkt_ber_preamble_from_actual_snr_mean = zeros(num_subcarriers, num_pkts);
    % pkt_ber_preamble_from_guessed_evm_mean = zeros(num_subcarriers, num_pkts);
    % pkt_ber_preamble_from_actual_evm_mean = zeros(num_subcarriers, num_pkts);

    % pkt_ber_preamble_from_actual_snr_calculation = zeros(num_subcarriers, num_pkts);
    % pkt_ber_preamble_from_guessed_evm_calculation = zeros(num_subcarriers, num_pkts);
    % pkt_ber_preamble_from_actual_evm_calculation = zeros(num_subcarriers, num_pkts);

    % % iii) from entire packet
    % pkt_ber_entire_from_actual_snr_mean = zeros(num_subcarriers, num_pkts);
    % pkt_ber_entire_from_guessed_evm_mean = zeros(num_subcarriers, num_pkts);
    % pkt_ber_entire_from_actual_evm_mean = zeros(num_subcarriers, num_pkts);

    % pkt_ber_entire_from_actual_snr_calculation = zeros(num_subcarriers, num_pkts);
    % pkt_ber_entire_from_guessed_evm_calculation = zeros(num_subcarriers, num_pkts);
    % pkt_ber_entire_from_actual_evm_calculation = zeros(num_subcarriers, num_pkts);

    % for pkt_i = 1:num_pkts
    % % for pkt_i = 1:1
    %     pkt_start_ind = (pkt_i-1) * num_ofdm_symbol_per_pkt + 1;
    %     pkt_end_ind   = pkt_i * num_ofdm_symbol_per_pkt;
    %     % fprintf('pkt: %d, %d~%d\n', pkt_i, pkt_start_ind, pkt_end_ind);


    %     % i) from tx and rx bit
    %     pkt_bit_start_ind = (pkt_i-1) * 2 * num_ofdm_symbol_per_pkt + 1;
    %     pkt_bit_end_ind   = pkt_i * 2 * num_ofdm_symbol_per_pkt;
    %     rx_bit_grid = rx_bit_grid(:, pkt_bit_start_ind:pkt_bit_end_ind);
    %     pkt_ber_from_rxtx_bit(:, pkt_i) = get_bit_error_rate(rx_bit_grid, tx_bit_grid);
    %     % fprintf('actual BER: ')
    %     % fprintf('%f, ', pkt_ber_from_rxtx_bit(:, pkt_i));
    %     % fprintf('\n');


    %     % ii) from preamble
    %     pkt_ber_preamble_from_actual_snr_mean(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_actual_snr(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     pkt_ber_preamble_from_guessed_evm_mean(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_guessed_evm(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     pkt_ber_preamble_from_actual_evm_mean(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_actual_evm(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     % fprintf('estimated BER from preamble: ')
    %     % fprintf('%f, ', pkt_ber_preamble_from_actual_snr_mean(:, pkt_i));
    %     % fprintf('\n');

    %     pkt_ber_preamble_from_actual_snr_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_actual_snr(:, pkt_start_ind:pkt_end_ind), 'calculation');
    %     pkt_ber_preamble_from_guessed_evm_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_guessed_evm(:, pkt_start_ind:pkt_end_ind), 'calculation');
    %     pkt_ber_preamble_from_actual_evm_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_preamble(ber_from_actual_evm(:, pkt_start_ind:pkt_end_ind), 'calculation');
        

    %     % iii) from entire packet
    %     pkt_ber_entire_from_actual_snr_mean(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_actual_snr(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     pkt_ber_entire_from_guessed_evm_mean(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_guessed_evm(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     pkt_ber_entire_from_actual_evm_mean(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_actual_evm(:, pkt_start_ind:pkt_end_ind), 'mean');
    %     % fprintf('estimated BER from entire packet: ')
    %     % fprintf('%f, ', pkt_ber_entire_from_actual_snr_mean(:, pkt_i));
    %     % fprintf('\n');


    %     pkt_ber_entire_from_actual_snr_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_actual_snr(:, pkt_start_ind:pkt_end_ind), 'calculation');
    %     pkt_ber_entire_from_guessed_evm_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_guessed_evm(:, pkt_start_ind:pkt_end_ind), 'calculation');
    %     pkt_ber_entire_from_actual_evm_calculation(:, pkt_i) = ...
    %         get_ber_of_packet_entire(ber_from_actual_evm(:, pkt_start_ind:pkt_end_ind), 'calculation');
    % end


    % %% ----------------------------------
    % % 5. evaluate 
    % fprintf('SC actual BER: %f\n', mean(pkt_ber_from_rxtx_bit, 2));

    % err_sc_pkt_ber_preamble_from_actual_snr_mean = abs(mean(pkt_ber_preamble_from_actual_snr_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_preamble_from_guessed_evm_mean = abs(mean(pkt_ber_preamble_from_guessed_evm_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_preamble_from_actual_evm_mean = abs(mean(pkt_ber_preamble_from_actual_evm_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % fprintf('error sc preamble_from_actual_snr_mean: %f\n', err_sc_pkt_ber_preamble_from_actual_snr_mean);
    % fprintf('error sc preamble_from_guessed_evm_mean: %f\n', err_sc_pkt_ber_preamble_from_guessed_evm_mean);
    % fprintf('error sc preamble_from_actual_evm_mean: %f\n', err_sc_pkt_ber_preamble_from_actual_evm_mean);

    % err_sc_pkt_ber_preamble_from_actual_snr_calculation = abs(mean(pkt_ber_preamble_from_actual_snr_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_preamble_from_guessed_evm_calculation = abs(mean(pkt_ber_preamble_from_guessed_evm_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_preamble_from_actual_evm_calculation = abs(mean(pkt_ber_preamble_from_actual_evm_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % fprintf('error sc preamble_from_actual_snr_calculation: %f\n', err_sc_pkt_ber_preamble_from_actual_snr_calculation);
    % fprintf('error sc preamble_from_guessed_evm_calculation: %f\n', err_sc_pkt_ber_preamble_from_guessed_evm_calculation);
    % fprintf('error sc preamble_from_actual_evm_calculation: %f\n', err_sc_pkt_ber_preamble_from_actual_evm_calculation);

    % err_sc_pkt_ber_entire_from_actual_snr_mean = abs(mean(pkt_ber_entire_from_actual_snr_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_entire_from_guessed_evm_mean = abs(mean(pkt_ber_entire_from_guessed_evm_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_entire_from_actual_evm_mean = abs(mean(pkt_ber_entire_from_actual_evm_mean, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % fprintf('error sc entire_from_actual_snr_mean: %f\n', err_sc_pkt_ber_entire_from_actual_snr_mean);
    % fprintf('error sc entire_from_guessed_evm_mean: %f\n', err_sc_pkt_ber_entire_from_guessed_evm_mean);
    % fprintf('error sc entire_from_actual_evm_mean: %f\n', err_sc_pkt_ber_entire_from_actual_evm_mean);

    % err_sc_pkt_ber_entire_from_actual_snr_calculation = abs(mean(pkt_ber_entire_from_actual_snr_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_entire_from_guessed_evm_calculation = abs(mean(pkt_ber_entire_from_guessed_evm_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % err_sc_pkt_ber_entire_from_actual_evm_calculation = abs(mean(pkt_ber_entire_from_actual_evm_calculation, 2) - mean(pkt_ber_from_rxtx_bit, 2) );
    % fprintf('error sc entire_from_actual_snr_calculation: %f\n', err_sc_pkt_ber_entire_from_actual_snr_calculation);
    % fprintf('error sc entire_from_guessed_evm_calculation: %f\n', err_sc_pkt_ber_entire_from_guessed_evm_calculation);
    % fprintf('error sc entire_from_actual_evm_calculation: %f\n', err_sc_pkt_ber_entire_from_actual_evm_calculation);


    % fprintf('all actual BER: %f\n', mean(mean(pkt_ber_from_rxtx_bit)));
    % fprintf('error all preamble_from_actual_snr_mean: %f\n', abs(mean(mean(pkt_ber_preamble_from_actual_snr_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) ) ;
    % fprintf('error all preamble_from_guessed_evm_mean: %f\n', abs(mean(mean(pkt_ber_preamble_from_guessed_evm_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );
    % fprintf('error all preamble_from_actual_evm_mean: %f\n', abs(mean(mean(pkt_ber_preamble_from_actual_evm_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );

    % fprintf('error all preamble_from_actual_snr_calculation: %f\n', abs(mean(mean(pkt_ber_preamble_from_actual_snr_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) ) ;
    % fprintf('error all preamble_from_guessed_evm_calculation: %f\n', abs(mean(mean(pkt_ber_preamble_from_guessed_evm_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );
    % fprintf('error all preamble_from_actual_evm_calculation: %f\n', abs(mean(mean(pkt_ber_preamble_from_actual_evm_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );

    % fprintf('error all entire_from_actual_snr_mean: %f\n', abs(mean(mean(pkt_ber_entire_from_actual_snr_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) ) ;
    % fprintf('error all entire_from_guessed_evm_mean: %f\n', abs(mean(mean(pkt_ber_entire_from_guessed_evm_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );
    % fprintf('error all entire_from_actual_evm_mean: %f\n', abs(mean(mean(pkt_ber_entire_from_actual_evm_mean) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );

    % fprintf('error all entire_from_actual_snr_calculation: %f\n', abs(mean(mean(pkt_ber_entire_from_actual_snr_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) ) ;
    % fprintf('error all entire_from_guessed_evm_calculation: %f\n', abs(mean(mean(pkt_ber_entire_from_guessed_evm_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );
    % fprintf('error all entire_from_actual_evm_calculation: %f\n', abs(mean(mean(pkt_ber_entire_from_actual_evm_calculation) ) - mean(mean(pkt_ber_from_rxtx_bit) ) ) );


    
    %


    % %% ----------------------------------
    % % calculate effective SNR for each packet
    % num_pkts = floor(size(rx_grid, 2) / num_ofdm_symbol_per_pkt);
    % fprintf('number of packets: %d\n', num_pkts);
    % pkt_snr_preamble = zeros(num_subcarriers, num_pkts);
    % pkt_snr_entire   = zeros(num_subcarriers, num_pkts);
    % pkt_bit_error_rate = zeros(num_subcarriers, num_pkts);
    
    % for pkt_i = 1:num_pkts
    %     pkt_start_ind = (pkt_i-1) * num_ofdm_symbol_per_pkt + 1;
    %     pkt_end_ind   = pkt_i * num_ofdm_symbol_per_pkt;

    %     snr_preamble = get_snr_of_packet_preamble(snr(:, pkt_start_ind:pkt_end_ind));
    %     pkt_snr_preamble(:, pkt_i) = snr_preamble;
    %     if length(snr_preamble) ~= num_subcarriers
    %         disp('wrong number of estimated SNR');
    %     end

    %     snr_entire = get_snr_of_packet_entire(snr(:, pkt_start_ind:pkt_end_ind));
    %     pkt_snr_entire(:, pkt_i) = snr_entire;
    %     if length(snr_entire) ~= num_subcarriers
    %         disp('wrong number of estimated SNR');
    %     end


    %     % rx_bit_grid = symbol2bit(modulation, rx_slice_grid(:, pkt_start_ind:pkt_end_ind));
    %     pkt_start_ind2 = (pkt_i-1) * 2 * num_ofdm_symbol_per_pkt + 1;
    %     pkt_end_ind2   = pkt_i * 2 * num_ofdm_symbol_per_pkt;
    %     rx_bit_grid = rx_bit_grid(:, pkt_start_ind2:pkt_end_ind2);
    %     pkt_bit_error_rate(:, pkt_i) = get_bit_error_rate(rx_bit_grid, tx_bit_grid);

    % end


    % %% ----------------------------------
    % % SNR prediction
    % alpha = 1;
    % ewma_preamble = ewma(pkt_snr_preamble, alpha);
    % ewma_entire   = ewma(pkt_snr_entire, alpha);
    % ewma_err_preamble = prediction_error(pkt_snr_preamble, ewma_preamble);
    % ewma_err_entire   = prediction_error(pkt_snr_entire, ewma_entire);


    % %% ----------------------------------
    % % best SNR shift that minimize BER estimation error
    % best_snr_shift_preamble = 0;
    % best_diff = 1;
    % for shift_i = -5:0.01:15
    %     this_diff = mean(abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - shift_i, method_snr2ber) - mean(pkt_bit_error_rate, 2) ) );
    %     if length(this_diff) > 1
    %         fprintf('wrong usage in getting diff: %d, %d\n', size(this_diff));
    %     end

    %     if this_diff < best_diff
    %         best_snr_shift_preamble = shift_i;
    %         best_diff = this_diff;
    %     end
    % end
    
    % best_snr_shift_entire = 0;
    % best_diff = 1;
    % for shift_i = -5:0.01:15
    %     this_diff = mean(abs(SNR2BER(modulation, mean(pkt_snr_entire, 2) - shift_i, method_snr2ber) - mean(pkt_bit_error_rate, 2) ) );
    %     if length(this_diff) > 1
    %         fprintf('wrong usage in getting diff: %d, %d\n', size(this_diff));
    %     end

    %     if this_diff < best_diff
    %         best_snr_shift_entire = shift_i;
    %         best_diff = this_diff;
    %     end
    % end


    % %% ----------------------------------
    % % write to file
    % ber_diff_preamble = abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2), method_snr2ber) - mean(pkt_bit_error_rate, 2));
    % ber_diff_entire   = abs(SNR2BER(modulation, mean(pkt_snr_entire, 2), method_snr2ber)   - mean(pkt_bit_error_rate, 2));
    % ber_diff_preamble_shift = abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - best_snr_shift_preamble, method_snr2ber) - mean(pkt_bit_error_rate, 2));
    % ber_diff_entire_shift   = abs(SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire, method_snr2ber) - mean(pkt_bit_error_rate, 2));
    
        
    % fh = fopen([output_dir input_sym_file ewma_error_file], 'w');
    % fprintf(fh, '# <preable prediction error>, <entire prediction error>, <actuall ber>, <preamble: BER diff> <entire: BER diff> <preamble shift: BER diff> <entire shift: BER diff>\n');
    % for sc_i = 1:num_subcarriers
    %     fprintf(fh, '%f, %f, %f, %f, %f, %f, %f\n', ...
    %         mean(ewma_err_preamble(sc_i, :), 2), ...
    %         mean(ewma_err_entire(sc_i, :), 2), ...
    %         mean(pkt_bit_error_rate(sc_i, :), 2), ... 
    %         ber_diff_preamble(sc_i, :), ...
    %         ber_diff_entire(sc_i, :), ...
    %         ber_diff_preamble_shift(sc_i, :), ...
    %         ber_diff_entire_shift(sc_i, :) );
    % end
    % R = corrcoef(SNR2BER(modulation, mean(pkt_snr_preamble, 2), method_snr2ber), mean(pkt_bit_error_rate, 2));
    % fprintf(fh, 'corrcoef between SNR_preamble and actual BER: %f\n', R(1, 2));
    % R = corrcoef(SNR2BER(modulation, mean(pkt_snr_entire, 2), method_snr2ber), mean(pkt_bit_error_rate, 2));
    % fprintf(fh, 'corrcoef between SNR_entire and actual BER: %f\n', R(1, 2));

    % R = corrcoef(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - best_snr_shift_preamble, method_snr2ber), mean(pkt_bit_error_rate, 2));
    % fprintf(fh, 'corrcoef between SNR_preamble (shift) and actual BER: %f\n', R(1, 2));
    % R = corrcoef(SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire, method_snr2ber), mean(pkt_bit_error_rate, 2));
    % fprintf(fh, 'corrcoef between SNR_entire (shift) and actual BER: %f\n', R(1, 2));

    % fprintf(fh, 'best SNR shift for preamble: %f\n', best_snr_shift_preamble);
    % fprintf(fh, 'best SNR shift for entire: %f\n', best_snr_shift_entire);

    % fclose(fh);


    % %% ----------------------------------
    % % actual BER for all packets
    % dlmwrite([output_dir input_sym_file '.evm.txt'], evm');
    % dlmwrite([output_dir input_sym_file '.estimated_snr.txt'], snr');
    % dlmwrite([output_dir input_sym_file '.actual_snr.txt'], actual_snr');
    % dlmwrite([output_dir input_sym_file '.ber.txt'], [pkt_bit_error_rate'; mean(pkt_bit_error_rate')]);
    % dlmwrite([output_dir input_sym_file '.actual_snr2ber.txt'], snr_err_rate);
    




    % %% ----------------------------------
    % % plot
    % plot_subcarrier_start = 1;
    % plot_subcarrier_end   = 48;
    % plot_symbol_start = 1;
    % plot_symbol_end   = 1000;
    
    % %% ----------------------------------
    % %   fig 0. csi over time
    % %
    % f0 = figure;
    % plot(abs(data_cpx(1:end)));
    % xlabel('symbols');
    % ylabel('magnitude');
    % print(f0, '-dpng', [figure_dir input_sym_file '.allmag.png']);
    % %


    % %% ----------------------------------
    % % DEBUG
    % %{
    % disp('pkt size');
    % dbg_index = find(abs(data_cpx(5000:end)) > 0);
    % size(dbg_index)
    % disp('pkt size 2');
    % dbg_index = find(abs(data_cpx(15000:end)) > 0);
    % size(dbg_index)
    % disp('pkt interval');
    % dbg_index = find(abs(data_cpx(2000:20000)) == 0);
    % size(dbg_index)
    % disp('pkt interval 2');
    % dbg_index = find(abs(data_cpx(10000:20000)) == 0);
    % size(dbg_index)
    % %}


    % %% ----------------------------------
    % %   fig 1. subcarrier csi over time
    % f1 = figure;
    % plot(abs(rx_grid)');
    % xlabel('symbols');
    % ylabel('magnitude');
    % % print(f1, '-dpsc', [figure_dir input_sym_file '.mag.ps']);
    % print(f1, '-dpng', [figure_dir input_sym_file '.magnitude.png']);

    % f1a = figure;
    % plot(abs(rx_grid([2, 21], 1:1000))');
    % xlabel('symbols');
    % ylabel('magnitude');
    % legend('subcarrier 1', 'subcarrier 21');
    % print(f1a, '-dpsc', [figure_dir input_sym_file '.magnitude_part.ps']);
    % % print(f1a, '-dpng', [figure_dir input_sym_file '.magnitude_part.png']);

    

    % %% ----------------------------------
    % % plot
    % %   fig 2. IQ plane
    % f2 = figure;
    % plot(real(rx_grid(1, 1:500)), imag(rx_grid(1, 1:500)), 'b.', ...
    %      real(qpsk_table), imag(qpsk_table), 'r*' );
    % xlabel('I');
    % ylabel('Q');
    % print(f2, '-dpsc', [figure_dir input_sym_file '.IQplane.ps']);
    % % print(f2, '-dpng', [figure_dir input_sym_file '.IQplane.png']);


    % %% ----------------------------------
    % % plot
    % %   fig 3. EVM
    % f3 = figure;
    % plot(evm(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    % xlabel('symbols');
    % ylabel('EVM');
    % print(f3, '-dpng', [figure_dir input_sym_file '.evm.png']);

    % f3a = figure;
    % plot(evm([1, 21], 1:200)' );
    % xlabel('symbols');
    % ylabel('EVM');
    % print(f3a, '-dpsc', [figure_dir input_sym_file '.evm_part.ps']);



    % %% ----------------------------------
    % % plot
    % %   fig 4. SNR
    % f4 = figure;
    % plot(snr(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    % xlabel('symbols');
    % ylabel('SNR');
    % print(f4, '-dpng', [figure_dir input_sym_file '.snr.png']);

    % f4a = figure;
    % plot(1:3*num_ofdm_symbol_per_pkt, snr([1], 1:3*num_ofdm_symbol_per_pkt)', 'bo-', ...
    %      [num_ofdm_symbol_per_pkt, num_ofdm_symbol_per_pkt], [0, 50], 'r-', ...
    %      [2*num_ofdm_symbol_per_pkt, 2*num_ofdm_symbol_per_pkt], [0, 50], 'r-');
    % xlabel('symbols');
    % ylabel('SNR');
    % print(f4a, '-dpsc', [figure_dir input_sym_file '.snr_part.ps']);

    % f4b = figure;
    % plot(snr([1 21], 1:200)');
    % xlabel('symbols');
    % ylabel('SNR (dB)');
    % legend('subcarrier 1', 'subcarrier 21')
    % print(f4b, '-dpsc', [figure_dir input_sym_file '.snr_part2.ps']);


    % %% ----------------------------------
    % % plot
    % %   fig 5. BER
    % f5 = figure;
    % plot(ber(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    % xlabel('symbols');
    % ylabel('BER');
    % print(f5, '-dpng', [figure_dir input_sym_file '.ber.png']);

    % f5a = figure;
    % plot(1:3*num_ofdm_symbol_per_pkt, ber([1], 1:3*num_ofdm_symbol_per_pkt)', 'bo-', ...
    %      [num_ofdm_symbol_per_pkt, num_ofdm_symbol_per_pkt], [0, 1], 'r-', ...
    %      [2*num_ofdm_symbol_per_pkt, 2*num_ofdm_symbol_per_pkt], [0, 1], 'r-');
    % xlabel('symbols');
    % ylabel('BER');
    % print(f5a, '-dpsc', [figure_dir input_sym_file '.ber_part.ps']);


    % %% ----------------------------------
    % % plot
    % %   fig 6. ewma preamble
    % f6 = figure;
    % plot(1:size(pkt_snr_preamble, 2), pkt_snr_preamble(1, :), 'b', ...
    %      1:size(pkt_snr_preamble, 2), ewma_err_preamble(1, 1:size(pkt_snr_preamble, 2) ), 'r');
    %      % 1:size(pkt_snr_preamble, 2), ewma_preamble(1, 1:size(pkt_snr_preamble, 2) ), ...
    % xlabel('packet');
    % ylabel('SNR');
    % legend('SNR', 'error');
    % print(f6, '-dpsc', [figure_dir input_sym_file '.ewma_preamble.ps']);

    % f6a = figure;
    % % plot(1:num_subcarriers, pkt_snr_preamble(:, 100), 'b', ...
    % %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    % plot(1:num_subcarriers, mean(pkt_snr_preamble, 2), 'b', ...
    %      1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_preamble, 2), method_snr2ber) * 2000, 'r', ...
    %      1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    % xlabel('subcarrier');
    % ylabel('SNR');
    % legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    % print(f6a, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_preamble.ps']);

    % f6b = figure;
    % % plot(1:num_subcarriers, pkt_snr_preamble(:, 100), 'b', ...
    % %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    % plot(1:num_subcarriers, mean(pkt_snr_preamble, 2), 'b', ...
    %      1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_preamble, 2)-best_snr_shift_preamble, method_snr2ber) * 2000, 'r', ...
    %      1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    % xlabel('subcarrier');
    % ylabel('SNR');
    % legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    % print(f6b, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_preamble_shift.ps']);

    

    % %% ----------------------------------
    % % plot
    % %   fig 7. ewma entire
    % f7 = figure;
    % plot(1:size(pkt_snr_entire, 2), pkt_snr_entire(1, :), 'b', ...
    %      1:size(pkt_snr_entire, 2), ewma_err_entire(1, 1:size(pkt_snr_entire, 2) ), 'r');
    %      % 1:size(pkt_snr_entire, 2), ewma_entire(1, 1:size(pkt_snr_entire, 2) ), ...
    % xlabel('packet');
    % ylabel('SNR');
    % legend('SNR', 'error');
    % print(f7, '-dpsc', [figure_dir input_sym_file '.ewma_entire.ps']);

    % f7a = figure;
    % % plot(1:num_subcarriers, pkt_snr_entire(:, 100), 'b', ...
    % %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    % plot(1:num_subcarriers, mean(pkt_snr_entire, 2), 'b', ...
    %      1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_entire, 2), method_snr2ber) * 2000, 'r', ...
    %      1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    % xlabel('subcarrier');
    % ylabel('SNR');
    % legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    % print(f7a, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_entire.ps']);

    % f7b = figure;
    % % plot(1:num_subcarriers, pkt_snr_entire(:, 100), 'b', ...
    % %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    % plot(1:num_subcarriers, mean(pkt_snr_entire, 2), 'b', ...
    %      1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire, method_snr2ber) * 2000, 'r', ...
    %      1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    % xlabel('subcarrier');
    % ylabel('SNR');
    % legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    % print(f7b, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_entire_shift.ps']);

    

    % %% ----------------------------------
    % % plot
    % %   fig 8. actual pkt error rate
    % f8 = figure;
    % plot(pkt_bit_error_rate([1 21], :)' );
    % xlabel('packet');
    % ylabel('actual error rate');
    % print(f8, '-dpsc', [figure_dir input_sym_file '.actual_ber.ps']);

    % f8a = figure;
    % plot(mean(pkt_bit_error_rate, 2));
    % xlabel('subcarrier');
    % ylabel('actual error rate');
    % print(f8a, '-dpsc', [figure_dir input_sym_file '.subcarrier_actual_ber.ps']);
    

    % f8b = figure;
    % plot(1:num_subcarriers, mean(ewma_err_preamble, 2), ...
    %      1:num_subcarriers, mean(ewma_err_entire, 2) );
    % xlabel('subcarrier');
    % ylabel('SNR prediction error');
    % legend('preamble', 'entire');
    % print(f8b, '-dpsc', [figure_dir input_sym_file '.subcarrier_ewma_error.ps']);

    

    % %% ----------------------------------
    % % plot
    % %   fig 9.
    % f9 = figure;
    % plot(min_snr:gran_snr:max_snr, snr_err_rate, 'r*');
    % hold on;
    % x_values = -10:0.1:50;
    % y_values_bpsk = SNR2BER('BPSK', x_values, method_snr2ber);
    % y_values_qpsk = SNR2BER('QPSK', x_values, method_snr2ber);
    % y_values_16qam = SNR2BER('16QAM', x_values, method_snr2ber);
    % y_values_qpsk_shift = SNR2BER(modulation, x_values+1, method_snr2ber);
    % plot(x_values, y_values_qpsk, ...
    %      x_values, y_values_qpsk_shift);
    % xlabel('SNR (dB)');
    % ylabel('BER');
    % legend('actual', ...
    %        'QPSK', ...
    %        'QPSK shift');
    % axis([-10 15 0 0.5]);
    % print(f9, '-dpsc', [figure_dir input_sym_file '.actual_snr2ber.ps']);


    % %% ----------------------------------
    % % plot
    % %   fig 10.
    % f10 = figure;
    % plot(min_evm:gran_evm:max_evm, evm_err_rate, 'r*');
    % hold on;
    % x_values = 0:0.01:2;
    % y_values_qpsk = SNR2BER('QPSK', EVM2SNR(x_values), method_snr2ber);
    % y_values_qpsk_shift = SNR2BER('QPSK', EVM2SNR(x_values+1), method_snr2ber);
    % plot(x_values, y_values_qpsk);
    % xlabel('EVM');
    % ylabel('BER');
    % legend('actual', ...
    %        'QPSK');
    % % axis([-10 15 0 0.5]);
    % print(f10, '-dpsc', [figure_dir input_sym_file '.actual_evm2ber.ps']);

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










