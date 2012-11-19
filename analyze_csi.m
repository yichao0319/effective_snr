
function analyze_csi(input_sym_file, input_snr_file)
    % addpath('../matlab');


    REMOVE_OUTLIER = 0;
    REMOVE_RAWOFDM_EMPTY_PERIOD = 0;

    
    %% ----------------------------------
    % global variables
    input_dir = './PARSEDDATA/';
    % input_sym_file = 'RXDATA.static.d0.pdat';
    % input_sym_file = 'rxdata_tmp.pdat';
    % input_sym_file = 'tmp_rxdata_s16.short.pdat';
    % input_sym_file = 'tmp_rxdata_s32.short.pdat';
    % input_sym_file = 'tmp_rxdata_s240.short.pdat';
    % input_sym_file = 'tmp_rxdata_s128.short.pdat';
    % input_sym_file = 'tmp_rxdata_s200.short.pdat';
    % input_sym_file = 'tmp_rxdata_s400.short.pdat';
    % input_sym_file = 'RXDATA.mobile.1.short.pdat';
    % input_sym_file = 'RXDATA.static.s16.pdat';
    % input_sym_file = 'RXDATA.static.s32.short.pdat';
    % input_sym_file = 'RXDATA.static.s64.short.pdat';

    tx_data_file = 'tx-data-s240.pdat';

    figure_dir = './figures/';
    output_dir = './OUTPUT/';
    ewma_error_file = '.ewma_err.txt';

    num_subcarriers = 48;
    num_ofdm_symbol_per_pkt = 96;
    modulation = 'QPSK';
    
    qpsk_table = mod_table(modulation);

    max_snr = 50;
    min_snr = -20;
    gran_snr = 0.5;
    max_evm = 2;
    min_evm = 0;
    gran_evm = 0.05;


    
    %% ----------------------------------
    % main

    %% ----------------------------------
    % load tx data
    %   format: <magnitude> <phase> <real> <image>
    tx_data = load([input_dir tx_data_file]);
    tx_data_cpx = tx_data(:, 3) + tx_data(:, 4) * i;
    tx_data_cpx = tx_data_cpx(1:num_subcarriers*num_ofdm_symbol_per_pkt);
    tx_grid = zeros(num_subcarriers, num_ofdm_symbol_per_pkt);
    tx_grid = reshape(tx_data_cpx(1:num_subcarriers*num_ofdm_symbol_per_pkt), num_subcarriers, []);
    tx_bit_grid = symbol2bit(modulation, tx_grid);
    disp(sprintf('load tx data size = %d, %d', size(tx_data_cpx) ) );
    if size(tx_grid, 2) ~= num_ofdm_symbol_per_pkt
        disp('wrong number of tx grid');
        return;
    end


    %% ----------------------------------
    % load SNR reported from rawOFDM
    %   format: <index> <SNR>
    raw_snr_data = load([input_dir input_snr_file]);
    snr_data = raw_snr_data(:, 2);
    

    %% ----------------------------------
    % load data to rx_grid
    %   format: <index> <magnitude> <phase> <real> <image>
    data = load([input_dir input_sym_file]);
    data_cpx = data(:, 4) + data(:, 5) * i;
    num_symbols = size(data, 1);
    disp(sprintf('load data size = %d, %d', size(data) ) );
    % if num_symbols ~= length(snr_data)
    %     fprintf('problem in size of RX symbols (%d) and SNRs (%d) !!!!!!!!!!!!!\n', num_symbols, length(snr_data));
    %     return;
    % end


    %% ----------------------------------
    % remove outlier
    if REMOVE_OUTLIER == 1
        outlier_index = find(abs(data_cpx > 5) );
        data_cpx(outlier_index) = 0.000001;
    end


    %% ----------------------------------
    % put rx symbols into grid
    rx_grid = zeros(num_subcarriers, num_symbols/num_subcarriers);
    rx_grid = reshape(data_cpx, num_subcarriers, []);
    disp(sprintf('grid size: %d, %d', size(rx_grid) ));


    %% ----------------------------------
    % remove empty generated by rawOFDM
    if REMOVE_RAWOFDM_EMPTY_PERIOD == 1
        sum_mag_ofdm_sym = sum(rx_grid);
        if(length(sum_mag_ofdm_sym) ~= num_symbols/num_subcarriers) 
            disp('sum_mag_ofdm_sym length error');
            return;
        end
        non_empty_ofdm_sym_ind = find(sum_mag_ofdm_sym ~= 0);
        rx_grid2 = rx_grid(:, non_empty_ofdm_sym_ind);
        rx_grid = rx_grid2;
        disp(sprintf('grid size after remove empty ofdm symbols: %d, %d', size(rx_grid) ));
    end


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
    tmp_snr2ber = SNR2BER(modulation, snr_data(1:48*96));
    dlmwrite([output_dir input_sym_file '.tmp_snr2ber.txt'], tmp_snr2ber');
    tmp_avg_snr = mean(snr_data(1:48*96));
    fprintf('avg snr: %.10f\n', tmp_avg_snr);
    tmp_avg_snr2ber = SNR2BER(modulation, tmp_avg_snr);
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
    %
    % 1. get EVM
    [out_grid, evm] = slice(modulation, rx_grid);
    out_bit_grid = symbol2bit(modulation, out_grid);
    evm = calculate_EVM(rx_grid, tx_grid);
    actual_snr = calculate_SNR(rx_grid, tx_grid);
    % 2. Map EVM to SNR
    snr = EVM2SNR(evm);
    % 3. Map SNR to BER
    ber = SNR2BER(modulation, snr);
    disp(sprintf('BER size: %d, %d', size(ber) ));

    snr_err_rate = cal_actual_SNR2BER(snr, out_bit_grid, tx_bit_grid, min_snr, max_snr, gran_snr);
    evm_err_rate = cal_actual_EVM2BER(evm, out_bit_grid, tx_bit_grid, min_evm, max_evm, gran_evm);
    %


    %% ----------------------------------
    % calculate effective SNR for each packet
    num_pkts = floor(size(rx_grid, 2) / num_ofdm_symbol_per_pkt);
    fprintf('number of packets: %d\n', num_pkts);
    pkt_snr_preamble = zeros(num_subcarriers, num_pkts);
    pkt_snr_entire   = zeros(num_subcarriers, num_pkts);
    pkt_bit_error_rate = zeros(num_subcarriers, num_pkts);
    
    for pkt_i = 1:num_pkts
        pkt_start_ind = (pkt_i-1) * num_ofdm_symbol_per_pkt + 1;
        pkt_end_ind   = pkt_i * num_ofdm_symbol_per_pkt;

        snr_preamble = get_snr_of_packet_preamble(snr(:, pkt_start_ind:pkt_end_ind));
        pkt_snr_preamble(:, pkt_i) = snr_preamble;
        if length(snr_preamble) ~= num_subcarriers
            disp('wrong number of estimated SNR');
        end

        snr_entire = get_snr_of_packet_entire(snr(:, pkt_start_ind:pkt_end_ind));
        pkt_snr_entire(:, pkt_i) = snr_entire;
        if length(snr_entire) ~= num_subcarriers
            disp('wrong number of estimated SNR');
        end


        % rx_bit_grid = symbol2bit(modulation, out_grid(:, pkt_start_ind:pkt_end_ind));
        pkt_start_ind2 = (pkt_i-1) * 2 * num_ofdm_symbol_per_pkt + 1;
        pkt_end_ind2   = pkt_i * 2 * num_ofdm_symbol_per_pkt;
        rx_bit_grid = out_bit_grid(:, pkt_start_ind2:pkt_end_ind2);
        pkt_bit_error_rate(:, pkt_i) = get_bit_error_rate(rx_bit_grid, tx_bit_grid);

    end


    %% ----------------------------------
    % SNR prediction
    alpha = 1;
    ewma_preamble = ewma(pkt_snr_preamble, alpha);
    ewma_entire   = ewma(pkt_snr_entire, alpha);
    ewma_err_preamble = prediction_error(pkt_snr_preamble, ewma_preamble);
    ewma_err_entire   = prediction_error(pkt_snr_entire, ewma_entire);


    %% ----------------------------------
    % best SNR shift that minimize BER estimation error
    best_snr_shift_preamble = 0;
    best_diff = 1;
    for shift_i = -5:0.01:15
        this_diff = mean(abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - shift_i) - mean(pkt_bit_error_rate, 2) ) );
        if length(this_diff) > 1
            fprintf('wrong usage in getting diff: %d, %d\n', size(this_diff));
        end

        if this_diff < best_diff
            best_snr_shift_preamble = shift_i;
            best_diff = this_diff;
        end
    end
    
    best_snr_shift_entire = 0;
    best_diff = 1;
    for shift_i = -5:0.01:15
        this_diff = mean(abs(SNR2BER(modulation, mean(pkt_snr_entire, 2) - shift_i) - mean(pkt_bit_error_rate, 2) ) );
        if length(this_diff) > 1
            fprintf('wrong usage in getting diff: %d, %d\n', size(this_diff));
        end

        if this_diff < best_diff
            best_snr_shift_entire = shift_i;
            best_diff = this_diff;
        end
    end


    %% ----------------------------------
    % write to file
    ber_diff_preamble = abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2) ) - mean(pkt_bit_error_rate, 2));
    ber_diff_entire   = abs(SNR2BER(modulation, mean(pkt_snr_entire, 2) )   - mean(pkt_bit_error_rate, 2));
    ber_diff_preamble_shift = abs(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - best_snr_shift_preamble ) - mean(pkt_bit_error_rate, 2));
    ber_diff_entire_shift   = abs(SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire ) - mean(pkt_bit_error_rate, 2));
    
        
    fh = fopen([output_dir input_sym_file ewma_error_file], 'w');
    fprintf(fh, '# <preable prediction error>, <entire prediction error>, <actuall ber>, <preamble: BER diff> <entire: BER diff> <preamble shift: BER diff> <entire shift: BER diff>\n');
    for sc_i = 1:num_subcarriers
        fprintf(fh, '%f, %f, %f, %f, %f, %f, %f\n', ...
            mean(ewma_err_preamble(sc_i, :), 2), ...
            mean(ewma_err_entire(sc_i, :), 2), ...
            mean(pkt_bit_error_rate(sc_i, :), 2), ... 
            ber_diff_preamble(sc_i, :), ...
            ber_diff_entire(sc_i, :), ...
            ber_diff_preamble_shift(sc_i, :), ...
            ber_diff_entire_shift(sc_i, :) );
    end
    R = corrcoef(SNR2BER(modulation, mean(pkt_snr_preamble, 2) ), mean(pkt_bit_error_rate, 2));
    fprintf(fh, 'corrcoef between SNR_preamble and actual BER: %f\n', R(1, 2));
    R = corrcoef(SNR2BER(modulation, mean(pkt_snr_entire, 2) ), mean(pkt_bit_error_rate, 2));
    fprintf(fh, 'corrcoef between SNR_entire and actual BER: %f\n', R(1, 2));

    R = corrcoef(SNR2BER(modulation, mean(pkt_snr_preamble, 2) - best_snr_shift_preamble ), mean(pkt_bit_error_rate, 2));
    fprintf(fh, 'corrcoef between SNR_preamble (shift) and actual BER: %f\n', R(1, 2));
    R = corrcoef(SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire ), mean(pkt_bit_error_rate, 2));
    fprintf(fh, 'corrcoef between SNR_entire (shift) and actual BER: %f\n', R(1, 2));

    fprintf(fh, 'best SNR shift for preamble: %f\n', best_snr_shift_preamble);
    fprintf(fh, 'best SNR shift for entire: %f\n', best_snr_shift_entire);

    fclose(fh);


    %% ----------------------------------
    % actual BER for all packets
    dlmwrite([output_dir input_sym_file '.evm.txt'], evm');
    dlmwrite([output_dir input_sym_file '.estimated_snr.txt'], snr');
    dlmwrite([output_dir input_sym_file '.actual_snr.txt'], actual_snr');
    dlmwrite([output_dir input_sym_file '.ber.txt'], [pkt_bit_error_rate'; mean(pkt_bit_error_rate')]);
    dlmwrite([output_dir input_sym_file '.actual_snr2ber.txt'], snr_err_rate);
    




    %% ----------------------------------
    % plot
    plot_subcarrier_start = 1;
    plot_subcarrier_end   = 48;
    plot_symbol_start = 1;
    plot_symbol_end   = 1000;
    
    %% ----------------------------------
    %   fig 0. csi over time
    %
    f0 = figure;
    plot(abs(data_cpx(1:end)));
    xlabel('symbols');
    ylabel('magnitude');
    print(f0, '-dpng', [figure_dir input_sym_file '.allmag.png']);
    %


    %% ----------------------------------
    % DEBUG
    %{
    disp('pkt size');
    dbg_index = find(abs(data_cpx(5000:end)) > 0);
    size(dbg_index)
    disp('pkt size 2');
    dbg_index = find(abs(data_cpx(15000:end)) > 0);
    size(dbg_index)
    disp('pkt interval');
    dbg_index = find(abs(data_cpx(2000:20000)) == 0);
    size(dbg_index)
    disp('pkt interval 2');
    dbg_index = find(abs(data_cpx(10000:20000)) == 0);
    size(dbg_index)
    %}


    %% ----------------------------------
    %   fig 1. subcarrier csi over time
    f1 = figure;
    plot(abs(rx_grid)');
    xlabel('symbols');
    ylabel('magnitude');
    % print(f1, '-dpsc', [figure_dir input_sym_file '.mag.ps']);
    print(f1, '-dpng', [figure_dir input_sym_file '.magnitude.png']);

    f1a = figure;
    plot(abs(rx_grid([2, 21], 1:1000))');
    xlabel('symbols');
    ylabel('magnitude');
    legend('subcarrier 1', 'subcarrier 21');
    print(f1a, '-dpsc', [figure_dir input_sym_file '.magnitude_part.ps']);
    % print(f1a, '-dpng', [figure_dir input_sym_file '.magnitude_part.png']);

    

    %% ----------------------------------
    % plot
    %   fig 2. IQ plane
    f2 = figure;
    plot(real(rx_grid(1, 1:500)), imag(rx_grid(1, 1:500)), 'b.', ...
         real(qpsk_table), imag(qpsk_table), 'r*' );
    xlabel('I');
    ylabel('Q');
    print(f2, '-dpsc', [figure_dir input_sym_file '.IQplane.ps']);
    % print(f2, '-dpng', [figure_dir input_sym_file '.IQplane.png']);


    %% ----------------------------------
    % plot
    %   fig 3. EVM
    f3 = figure;
    plot(evm(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    xlabel('symbols');
    ylabel('EVM');
    print(f3, '-dpng', [figure_dir input_sym_file '.evm.png']);

    f3a = figure;
    plot(evm([1, 21], 1:200)' );
    xlabel('symbols');
    ylabel('EVM');
    print(f3a, '-dpsc', [figure_dir input_sym_file '.evm_part.ps']);



    %% ----------------------------------
    % plot
    %   fig 4. SNR
    f4 = figure;
    plot(snr(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    xlabel('symbols');
    ylabel('SNR');
    print(f4, '-dpng', [figure_dir input_sym_file '.snr.png']);

    f4a = figure;
    plot(1:3*num_ofdm_symbol_per_pkt, snr([1], 1:3*num_ofdm_symbol_per_pkt)', 'bo-', ...
         [num_ofdm_symbol_per_pkt, num_ofdm_symbol_per_pkt], [0, 50], 'r-', ...
         [2*num_ofdm_symbol_per_pkt, 2*num_ofdm_symbol_per_pkt], [0, 50], 'r-');
    xlabel('symbols');
    ylabel('SNR');
    print(f4a, '-dpsc', [figure_dir input_sym_file '.snr_part.ps']);

    f4b = figure;
    plot(snr([1 21], 1:200)');
    xlabel('symbols');
    ylabel('SNR (dB)');
    legend('subcarrier 1', 'subcarrier 21')
    print(f4b, '-dpsc', [figure_dir input_sym_file '.snr_part2.ps']);


    %% ----------------------------------
    % plot
    %   fig 5. BER
    f5 = figure;
    plot(ber(plot_subcarrier_start:plot_subcarrier_end, plot_symbol_start:plot_symbol_end)' );
    xlabel('symbols');
    ylabel('BER');
    print(f5, '-dpng', [figure_dir input_sym_file '.ber.png']);

    f5a = figure;
    plot(1:3*num_ofdm_symbol_per_pkt, ber([1], 1:3*num_ofdm_symbol_per_pkt)', 'bo-', ...
         [num_ofdm_symbol_per_pkt, num_ofdm_symbol_per_pkt], [0, 1], 'r-', ...
         [2*num_ofdm_symbol_per_pkt, 2*num_ofdm_symbol_per_pkt], [0, 1], 'r-');
    xlabel('symbols');
    ylabel('BER');
    print(f5a, '-dpsc', [figure_dir input_sym_file '.ber_part.ps']);


    %% ----------------------------------
    % plot
    %   fig 6. ewma preamble
    f6 = figure;
    plot(1:size(pkt_snr_preamble, 2), pkt_snr_preamble(1, :), 'b', ...
         1:size(pkt_snr_preamble, 2), ewma_err_preamble(1, 1:size(pkt_snr_preamble, 2) ), 'r');
         % 1:size(pkt_snr_preamble, 2), ewma_preamble(1, 1:size(pkt_snr_preamble, 2) ), ...
    xlabel('packet');
    ylabel('SNR');
    legend('SNR', 'error');
    print(f6, '-dpsc', [figure_dir input_sym_file '.ewma_preamble.ps']);

    f6a = figure;
    % plot(1:num_subcarriers, pkt_snr_preamble(:, 100), 'b', ...
    %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    plot(1:num_subcarriers, mean(pkt_snr_preamble, 2), 'b', ...
         1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_preamble, 2)) * 2000, 'r', ...
         1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    xlabel('subcarrier');
    ylabel('SNR');
    legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    print(f6a, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_preamble.ps']);

    f6b = figure;
    % plot(1:num_subcarriers, pkt_snr_preamble(:, 100), 'b', ...
    %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    plot(1:num_subcarriers, mean(pkt_snr_preamble, 2), 'b', ...
         1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_preamble, 2)-best_snr_shift_preamble) * 2000, 'r', ...
         1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    xlabel('subcarrier');
    ylabel('SNR');
    legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    print(f6b, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_preamble_shift.ps']);

    

    %% ----------------------------------
    % plot
    %   fig 7. ewma entire
    f7 = figure;
    plot(1:size(pkt_snr_entire, 2), pkt_snr_entire(1, :), 'b', ...
         1:size(pkt_snr_entire, 2), ewma_err_entire(1, 1:size(pkt_snr_entire, 2) ), 'r');
         % 1:size(pkt_snr_entire, 2), ewma_entire(1, 1:size(pkt_snr_entire, 2) ), ...
    xlabel('packet');
    ylabel('SNR');
    legend('SNR', 'error');
    print(f7, '-dpsc', [figure_dir input_sym_file '.ewma_entire.ps']);

    f7a = figure;
    % plot(1:num_subcarriers, pkt_snr_entire(:, 100), 'b', ...
    %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    plot(1:num_subcarriers, mean(pkt_snr_entire, 2), 'b', ...
         1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_entire, 2) ) * 2000, 'r', ...
         1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    xlabel('subcarrier');
    ylabel('SNR');
    legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    print(f7a, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_entire.ps']);

    f7b = figure;
    % plot(1:num_subcarriers, pkt_snr_entire(:, 100), 'b', ...
    %      1:num_subcarriers, pkt_bit_error_rate(:, 100) * 1000, 'r');
    plot(1:num_subcarriers, mean(pkt_snr_entire, 2), 'b', ...
         1:num_subcarriers, SNR2BER(modulation, mean(pkt_snr_entire, 2) - best_snr_shift_entire ) * 2000, 'r', ...
         1:num_subcarriers, mean(pkt_bit_error_rate, 2) * 2000, 'g');
    xlabel('subcarrier');
    ylabel('SNR');
    legend('SNR', 'estimated BER(*2000)', 'actual BER(*2000)');
    print(f7b, '-dpsc', [figure_dir input_sym_file '.subcarrier_snr_entire_shift.ps']);

    

    %% ----------------------------------
    % plot
    %   fig 8. actual pkt error rate
    f8 = figure;
    plot(pkt_bit_error_rate([1 21], :)' );
    xlabel('packet');
    ylabel('actual error rate');
    print(f8, '-dpsc', [figure_dir input_sym_file '.actual_ber.ps']);

    f8a = figure;
    plot(mean(pkt_bit_error_rate, 2));
    xlabel('subcarrier');
    ylabel('actual error rate');
    print(f8a, '-dpsc', [figure_dir input_sym_file '.subcarrier_actual_ber.ps']);
    

    f8b = figure;
    plot(1:num_subcarriers, mean(ewma_err_preamble, 2), ...
         1:num_subcarriers, mean(ewma_err_entire, 2) );
    xlabel('subcarrier');
    ylabel('SNR prediction error');
    legend('preamble', 'entire');
    print(f8b, '-dpsc', [figure_dir input_sym_file '.subcarrier_ewma_error.ps']);

    

    %% ----------------------------------
    % plot
    %   fig 9.
    f9 = figure;
    plot(min_snr:gran_snr:max_snr, snr_err_rate, 'r*');
    hold on;
    x_values = -10:0.1:50;
    y_values_bpsk = SNR2BER('BPSK', x_values);
    y_values_qpsk = SNR2BER('QPSK', x_values);
    y_values_16qam = SNR2BER('16QAM', x_values);
    y_values_qpsk_shift = SNR2BER(modulation, x_values+1);
    plot(x_values, y_values_qpsk, ...
         x_values, y_values_qpsk_shift);
    xlabel('SNR (dB)');
    ylabel('BER');
    legend('actual', ...
           'QPSK', ...
           'QPSK shift');
    axis([-10 15 0 0.5]);
    print(f9, '-dpsc', [figure_dir input_sym_file '.actual_snr2ber.ps']);


    %% ----------------------------------
    % plot
    %   fig 10.
    f10 = figure;
    plot(min_evm:gran_evm:max_evm, evm_err_rate, 'r*');
    hold on;
    x_values = 0:0.01:2;
    y_values_qpsk = SNR2BER('QPSK', EVM2SNR(x_values));
    y_values_qpsk_shift = SNR2BER('QPSK', EVM2SNR(x_values+1));
    plot(x_values, y_values_qpsk);
    xlabel('EVM');
    ylabel('BER');
    legend('actual', ...
           'QPSK');
    % axis([-10 15 0 0.5]);
    print(f10, '-dpsc', [figure_dir input_sym_file '.actual_evm2ber.ps']);

end


%% ----------------------------------
% functions

%% get_snr_of_packet_preamble: function description
function [estimated_snr] = get_snr_of_packet_preamble(snr)
    estimated_snr = mean(snr(:, 1:5), 2);
end

%% get_snr_of_packet_entire: function description
function [estimated_snr] = get_snr_of_packet_entire(snr)
   estimated_snr = mean(snr, 2);
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
function [ber] = get_bit_error_rate(rx_grid, tx_grid)
    if size(rx_grid) ~= size(tx_grid) 
        disp('wrong pkt size');
    end

    ber = zeros(size(rx_grid, 1), 1);

    for sc_i = 1:size(rx_grid, 1)
        wrong_ind = find(rx_grid(sc_i, :) ~= tx_grid(sc_i, :) );
        num_err = length(wrong_ind);

        ber(sc_i, 1) = num_err / size(rx_grid, 2);
    end
end


%% calculate_EVM: function description
function [evm] = calculate_EVM(rx_grid, tx_grid)
    num_pkts = floor(size(rx_grid, 2) / size(tx_grid, 2));
    remain_syms = size(rx_grid, 2) - num_pkts * size(tx_grid, 2);
    tx_grid_rep = [repmat(tx_grid, 1, num_pkts) tx_grid(:, 1:remain_syms)];

    evm = abs(rx_grid - tx_grid_rep) / 1;
end


%% calculate_SNR: function description
function [snr] = calculate_SNR(rx_grid, tx_grid)
    num_pkts = floor(size(rx_grid, 2) / size(tx_grid, 2));
    remain_syms = size(rx_grid, 2) - num_pkts * size(tx_grid, 2);
    tx_grid_rep = [repmat(tx_grid, 1, num_pkts) tx_grid(:, 1:remain_syms)];

    snr = power(abs(tx_grid_rep), 2) ./ power(abs(rx_grid - tx_grid_rep), 2);
    snr = 10 * log(snr) / log(10);
end


%% cal_actual_SNR2BER: function description
function [snr_err_rate] = cal_actual_SNR2BER(snr, rx_bit_grid, tx_bit_grid, min_snr, max_snr, gran)
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
end


%% cal_actual_EVM2BER: function description
function [evm_err_rate] = cal_actual_EVM2BER(evm, rx_bit_grid, tx_bit_grid, min_evm, max_evm, gran)
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
end









