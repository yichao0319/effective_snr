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

%% plot_evm_snr_ber_conversion: function description
function plot_evm_snr_ber_conversion()


    %% ----------------------------------
    % constants
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

    figure_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/figures2/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT/';

    input_sym_files = {'RXDATA.static.d0.pdat', 'RXDATA.static.d6.pdat', 'RXDATA.mobile.2.pdat', 'RXDATA.mobile.b1.pdat', 'RXDATA.mobile.b2.pdat', 'RXDATA.s96.pdat', 'RXDATA.s96-2.pdat', 'RXDATA.static.highSNR.pdat', 'RXDATA.static.highSNR.2.pdat', 'RXDATA.static.highSNR.d1.pdat', 'RXDATA.static.highSNR.d1.2.pdat', 'RXDATA.static.highSNR.d2.pdat', 'RXDATA.static.highSNR.d6.pdat', 'RXDATA.static.highSNR.d6.2.pdat', 'RXDATA.mobile.highSNR.pdat', 'RXDATA.mobile.highSNR.1.pdat', 'RXDATA.mobile.highSNR.b1.pdat', 'RXDATA.mobile.highSNR.b2.pdat'};
    % input_sym_files = {'RXDATA.static.d0.pdat', 'RXDATA.static.d6.pdat'};

    num_subcarriers = 48;
    num_ofdm_symbol_per_pkt = 50;
    modulation = 'QPSK';
    snr_shift = 0;
    % if strcmp(method_snr2ber, SNR2BER_FORMULA)
    %     snr_shift = 1.5;
    % end
    num_preamble_ofdm_sym = 2;

    max_snr = 100;
    min_snr = -50;
    gran_snr = 0.5;
    max_evm = 100;
    min_evm = 0;
    gran_evm = 0.05;

    snr_err_cnt_all = zeros((max_snr - min_snr)/gran_snr + 1, 2);
    snr_err_rate_all = zeros((max_snr - min_snr)/gran_snr + 1, 1);
    evm_err_cnt_all = zeros((max_evm - min_evm)/gran_evm + 1, 2);
    evm_err_rate_all = zeros((max_evm - min_evm)/gran_evm + 1, 1);


    %% ----------------------------------
    % initinalization
    
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


    
    %% ----------------------------------
    % main

    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));

        %% ----------------------------------
        % load rx symbols to rx_grid
        %   format: <index> <magnitude> <phase> <real> <image>
        data = load([input_dir input_sym_file]);
        data_cpx = data(:, 4) + data(:, 5) * i;                             % e.g. (4320x * 1)
        num_symbols = size(data, 1);
        disp(sprintf('load data size = %d, %d', size(data) ) );


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
        




        %% ----------------------------------
        % 1. EVM and rx bits
        %   i) get EVM from closest constellation points
        [rx_slice_grid, guessed_evm] = slice(modulation, rx_grid);  % e.g. rx_slice_grid: (48 * 90x)
                                                                    % e.g. guessed_evm: (48 * 90x)
        rx_bit_grid = symbol2bit(modulation, rx_slice_grid);        % e.g. rx_bit_grid: (48 * 180x)
        %   ii) get EVM by assuming we know the correct constellation points 
        % actual_evm = calculate_actual_EVM(rx_grid, tx_grid);              % e.g. actual_evm: (48 * 90x)
        

        %% ----------------------------------
        % evaluate SNR to BER, EVM to BER curve
        actual_snr = calculate_actual_SNR(rx_grid, tx_grid) + snr_shift;
        actual_evm = calculate_actual_EVM(rx_grid, tx_grid);
        [snr_err_rate, snr_cnt] = cal_actual_SNR2BER(actual_snr, rx_bit_grid, tx_bit_grid, min_snr, max_snr, gran_snr);
        [evm_err_rate, evm_cnt] = cal_actual_EVM2BER(actual_evm, rx_bit_grid, tx_bit_grid, min_evm, max_evm, gran_evm);
        snr_err_cnt_all(:, 2) = snr_err_cnt_all(:, 2) + snr_cnt;
        snr_err_cnt_all(:, 1) = snr_err_cnt_all(:, 1) + snr_err_rate .* snr_cnt;
        evm_err_cnt_all(:, 2) = evm_err_cnt_all(:, 2) + evm_cnt;
        evm_err_cnt_all(:, 1) = evm_err_cnt_all(:, 1) + evm_err_rate .* evm_cnt;
        
        f9 = figure;
        plot(min_snr:gran_snr:max_snr, snr_err_rate, 'r*', ...
             min_snr:gran_snr:max_snr, snr_cnt/max(snr_cnt), 'g');
        hold on;
        x_values = min_snr:0.1:max_snr;
        y_values_qpsk_formula = SNR2BER('QPSK', x_values, SNR2BER_FORMULA);
        y_values_qpsk_threshold = SNR2BER('QPSK', x_values, SNR2BER_THRESHOLD);
        plot(x_values, y_values_qpsk_formula, ...
             x_values, y_values_qpsk_threshold);
        xlabel('SNR (dB)');
        ylabel('BER');
        legend('actual', ...
               ['count (/' int2str(max(snr_cnt)) ')'], ...
               'formula', 'threshold');
        axis([min_snr max_snr 0 1]);
        print(f9, '-dpsc', [figure_dir input_sym_file '.snr2ber.ps']);

        f10 = figure;
        plot(min_evm:gran_evm:max_evm, evm_err_rate, 'r*', ...
             min_evm:gran_evm:max_evm, evm_cnt/max(evm_cnt), 'g');
        hold on;
        x_values = min_evm:0.01:max_evm;
        y_values_qpsk_formula = SNR2BER('QPSK', EVM2SNR(x_values), SNR2BER_FORMULA);
        y_values_qpsk_threshold = SNR2BER('QPSK', EVM2SNR(x_values), SNR2BER_THRESHOLD);
        plot(x_values, y_values_qpsk_formula, ...
             x_values, y_values_qpsk_threshold);
        xlabel('EVM');
        ylabel('BER');
        legend('actual', ...
               ['count (/' int2str(max(evm_cnt)) ')'], ...
               'formula', 'threshold');
        axis([0 15 0 1]);
        print(f10, '-dpsc', [figure_dir input_sym_file '.evm2ber.ps']);
        % return;
        

    end % end for all files






    f11 = figure;
    for snr_i = 1:(max_snr - min_snr) / gran_snr + 1
        if(snr_err_cnt_all(snr_i, 1) > 0)
            snr_err_rate_all(snr_i, 1) = snr_err_cnt_all(snr_i, 1) / snr_err_cnt_all(snr_i, 2);
        end
    end
    plot(min_snr:gran_snr:max_snr, snr_err_rate_all, 'r*', ...
         min_snr:gran_snr:max_snr, snr_err_cnt_all(:, 2)/max(snr_err_cnt_all(:, 2)), 'g');
    hold on;
    x_values = min_snr:0.1:max_snr;
    y_values_qpsk_formula = SNR2BER('QPSK', x_values, SNR2BER_FORMULA);
    y_values_qpsk_threshold = SNR2BER('QPSK', x_values, SNR2BER_THRESHOLD);
    plot(x_values, y_values_qpsk_formula, ...
         x_values, y_values_qpsk_threshold);
    xlabel('SNR (dB)');
    ylabel('BER');
    legend('actual', ...
           ['count (/' int2str(max(snr_err_cnt_all(:, 2))) ')'], ...
           'formula', 'threshold');
    axis([min_snr max_snr 0 1]);
    print(f11, '-dpsc', [figure_dir 'snr2ber.ps']);

    f12 = figure;
    for evm_i = 1:(max_evm - min_evm) / gran_evm + 1
        if(evm_err_cnt_all(evm_i, 1) > 0)
            evm_err_rate_all(evm_i, 1) = evm_err_cnt_all(evm_i, 1) / evm_err_cnt_all(evm_i, 2);
        end
    end
    plot(min_evm:gran_evm:max_evm, evm_err_rate_all, 'r*', ...
         min_evm:gran_evm:max_evm, evm_err_cnt_all(:, 2)/max(evm_err_cnt_all(:, 2)), 'g');
    hold on;
    x_values = min_evm:0.01:max_evm;
    y_values_qpsk_formula = SNR2BER('QPSK', EVM2SNR(x_values), SNR2BER_FORMULA);
    y_values_qpsk_threshold = SNR2BER('QPSK', EVM2SNR(x_values), SNR2BER_THRESHOLD);
    plot(x_values, y_values_qpsk_formula, ...
         x_values, y_values_qpsk_threshold);
    xlabel('EVM');
    ylabel('BER');
    legend('actual', ...
           ['count (/' int2str(max(evm_err_cnt_all(:, 2))) ')'], ...
           'formula', 'threshold');
    axis([0 15 0 1]);
    print(f12, '-dpsc', [figure_dir 'evm2ber.ps']);

end


%% ----------------------------------
% functions



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










