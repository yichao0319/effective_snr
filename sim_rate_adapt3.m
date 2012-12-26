%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%   file_ind: the index for the input file
%   method_SNR: the schemes to calculate SNR
%       preambleSNR
%       allSNR
%       allSNRoracle
%   method_RA: the schemes for Rate Adaptation
%       Oracle
%       Thresholding
%       Probability
%   method_PRED: the schemes for prediction
%       currPkt
%       prevPkt
%       EWMA
%       HW
%   threshold: when the PPR >= threshold, assume the packet can be received correctly
%
% Output:
%
%
% Example:
%   sim_rate_adapt3(1, 'preambleSNR', 'Thresholding', 'currPkt', 0.9)
%

function [tput] = sim_rate_adapt3(file_ind, method_SNR, method_RA, method_PRED, threshold)
    
    %% ----------------------------------
    % constants
    preambleSNR = 'preambleSNR';
    allSNR = 'allSNR';
    allSNRoracle = 'allSNRoracle';

    Oracle = 'Oracle';
    Thresholding = 'Thresholding';
    Probability_1 = 'Probability_1';
    Probability_2 = 'Probability_2';

    currPkt = 'currPkt';
    prevPkt = 'prevPkt';
    EWMA = 'EWMA';
    HW = 'HW';

    BPSK = 1;
    QPSK = 2;
    QAM16 = 3;
    QAM64 = 4;

    NUM_MOD = 4;
    NUM_EFF_SYM = 1;
    PRED_GRANULARITY = 0.1;

    %% configuration
    BLOCK_SIZE = 144;   % 144 bytes
    RATES_K = [1; 1; 3; 1; 3; 2; 3; 5];
    RATES_N = [2; 2; 4; 2; 4; 3; 4; 6];
    MOD = [BPSK; QPSK; QPSK; QAM16; QAM16; QAM64; QAM64; QAM64];
    TPUT_TABLE = [0.5; 1; 1.5; 2; 3; 4; 4.5; 5];
    

    CONST_BER_THRESHOLD_ALL = [0.001, 0.01, 0.004, 0.0013];
    CONST_BER_THRESHOLD_ALL2 = [0.001, 0.01, 0.004, 0.0025];
    CONST_BER_THRESHOLD_EFF = [0.001, 0.01, 0.005, 0.00025];
    CONST_BER_THRESHOLD_EVM = [0.001, 0.01, 0.0015, 0.0025];

    num_subcarriers = 48;
    num_ofdm_symbol_BPSK = 24;
    num_ofdm_symbol_QPSK = 24;
    num_ofdm_symbol_16QAM = 24;
    num_ofdm_symbol_64QAM = 23;
    num_ofdm_symbol = [24; 24; 24; 23];

    num_bits_per_pkt_BPSK = num_subcarriers * num_ofdm_symbol(1) * TPUT_TABLE(1);
    num_bits_per_pkt_QPSK = num_subcarriers * num_ofdm_symbol(2) * TPUT_TABLE(2);
    num_bits_per_pkt_16QAM = num_subcarriers * num_ofdm_symbol(3) * TPUT_TABLE(3);
    num_bits_per_pkt_64QAM = num_subcarriers * num_ofdm_symbol(4) * TPUT_TABLE(4);
    num_bits_per_pkt = [num_bits_per_pkt_BPSK; num_bits_per_pkt_QPSK; num_bits_per_pkt_16QAM; num_bits_per_pkt_64QAM];
    

    %% ----------------------------------
    % global variables
    %
    input_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/gen_traces/PARSEDDATA/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';
    actual_ber_file = ['rx_actual_bers_run' int2str(file_ind) '.dat'];
    
    ber_bpsk_file   = ['rx_bers_bpsk_run' int2str(file_ind) '.dat'];
    ber_qpsk_file   = ['rx_bers_qpsk_run' int2str(file_ind) '.dat'];
    ber_16qam_file  = ['rx_bers_16qam_run' int2str(file_ind) '.dat'];
    ber_64qam_file  = ['rx_bers_64qam_run' int2str(file_ind) '.dat'];

    evmber_bpsk_file   = ['rx_evmbers_bpsk_run' int2str(file_ind) '.dat'];
    evmber_qpsk_file   = ['rx_evmbers_qpsk_run' int2str(file_ind) '.dat'];
    evmber_16qam_file  = ['rx_evmbers_16qam_run' int2str(file_ind) '.dat'];
    evmber_64qam_file  = ['rx_evmbers_64qam_run' int2str(file_ind) '.dat'];


    %% ----------------------------------
    % initinalization
    %
    %  load actual ber for each packet: modulation * num_pkts
    actual_bers = load([input_dir actual_ber_file]);

    %  load actual SNR -> BERs of all symbols
    ber_bpsk_grid = load([input_dir ber_bpsk_file]);
    ber_qpsk_grid = load([input_dir ber_qpsk_file]);
    ber_16qam_grid = load([input_dir ber_16qam_file]);
    ber_64qam_grid = load([input_dir ber_64qam_file]);

    % load EVM SNR -> EVM BERs of all symbols
    evmber_bpsk_grid = load([input_dir evmber_bpsk_file]);
    evmber_qpsk_grid = load([input_dir evmber_qpsk_file]);
    evmber_16qam_grid = load([input_dir evmber_16qam_file]);
    evmber_64qam_grid = load([input_dir evmber_64qam_file]);


    %% ----------------------------------
    % main
    %

    %% ----------------------------------
    % get Effective SNR according to method_SNR
    %   @output: bers_BPSK, stdevs_BPSK,
    %            bers_QPSK, stdevs_QPSK, ...
    %
    % fprintf('- get Effective SNR according to method_SNR\n');
    if strcmp(method_SNR, preambleSNR) 
        [bers_BPSK,  stdevs_BPSK]  = cal_eff_ber(ber_bpsk_grid,  num_ofdm_symbol_BPSK,  NUM_EFF_SYM);
        [bers_QPSK,  stdevs_QPSK]  = cal_eff_ber(ber_qpsk_grid,  num_ofdm_symbol_QPSK,  NUM_EFF_SYM);
        [bers_16QAM, stdevs_16QAM] = cal_eff_ber(ber_16qam_grid, num_ofdm_symbol_16QAM, NUM_EFF_SYM);
        [bers_64QAM, stdevs_64QAM] = cal_eff_ber(ber_64qam_grid, num_ofdm_symbol_64QAM, NUM_EFF_SYM);
    elseif strcmp(method_SNR, allSNR) 
        [bers_BPSK,  stdevs_BPSK]  = cal_eff_ber(evmber_bpsk_grid,  num_ofdm_symbol_BPSK,  num_ofdm_symbol_BPSK);
        [bers_QPSK,  stdevs_QPSK]  = cal_eff_ber(evmber_qpsk_grid,  num_ofdm_symbol_QPSK,  num_ofdm_symbol_QPSK);
        [bers_16QAM, stdevs_16QAM] = cal_eff_ber(evmber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
        [bers_64QAM, stdevs_64QAM] = cal_eff_ber(evmber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    elseif strcmp(method_SNR, allSNRoracle) 
        [bers_BPSK,  stdevs_BPSK]  = cal_eff_ber(ber_bpsk_grid,  num_ofdm_symbol_BPSK,  num_ofdm_symbol_BPSK);
        [bers_QPSK,  stdevs_QPSK]  = cal_eff_ber(ber_qpsk_grid,  num_ofdm_symbol_QPSK,  num_ofdm_symbol_QPSK);
        [bers_16QAM, stdevs_16QAM] = cal_eff_ber(ber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
        [bers_64QAM, stdevs_64QAM] = cal_eff_ber(ber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    else
        error('wrong SNR method');
    end

    
    %% ----------------------------------
    % get the predicted BER mean and stdev
    %   @output: pred_bers_BPSK, pred_stdevs_BPSK,
    %            pred_bers_QPSK, pred_stdevs_QPSK, ...
    %
    % fprintf('- get the predicted BER mean and stdev\n');
    if strcmp(method_PRED, currPkt)
        pred_bers_BPSK  = bers_BPSK;
        pred_bers_QPSK  = bers_QPSK;
        pred_bers_16QAM = bers_16QAM;
        pred_bers_64QAM = bers_64QAM;

        pred_stdevs_BPSK  = stdevs_BPSK;
        pred_stdevs_QPSK  = stdevs_QPSK;
        pred_stdevs_16QAM = stdevs_16QAM;
        pred_stdevs_64QAM = stdevs_64QAM;

    elseif strcmp(method_PRED, prevPkt)
        pred_bers_BPSK  = [0, bers_BPSK(1, 1:end-1)];
        pred_bers_QPSK  = [0, bers_QPSK(1, 1:end-1)];
        pred_bers_16QAM = [0, bers_16QAM(1, 1:end-1)];
        pred_bers_64QAM = [0, bers_64QAM(1, 1:end-1)];

        pred_stdevs_BPSK  = [0, stdevs_BPSK(1, 1:end-1)];
        pred_stdevs_QPSK  = [0, stdevs_QPSK(1, 1:end-1)];
        pred_stdevs_16QAM = [0, stdevs_16QAM(1, 1:end-1)];
        pred_stdevs_64QAM = [0, stdevs_64QAM(1, 1:end-1)];

    elseif strcmp(method_PRED, EWMA)
        [alpha, pred_ts] = ewma_trial(bers_BPSK, bers_BPSK, PRED_GRANULARITY);
        pred_bers_BPSK = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(bers_QPSK, bers_QPSK, PRED_GRANULARITY);
        pred_bers_QPSK = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(bers_16QAM, bers_16QAM, PRED_GRANULARITY);
        pred_bers_16QAM = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(bers_64QAM, bers_64QAM, PRED_GRANULARITY);
        pred_bers_64QAM = pred_ts(1:end-1);

        [alpha, pred_ts] = ewma_trial(stdevs_BPSK, stdevs_BPSK, PRED_GRANULARITY);
        pred_stdevs_BPSK = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(stdevs_QPSK, stdevs_QPSK, PRED_GRANULARITY);
        pred_stdevs_QPSK = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(stdevs_16QAM, stdevs_16QAM, PRED_GRANULARITY);
        pred_stdevs_16QAM = pred_ts(1:end-1);
        [alpha, pred_ts] = ewma_trial(stdevs_64QAM, stdevs_64QAM, PRED_GRANULARITY);
        pred_stdevs_64QAM = pred_ts(1:end-1);

    elseif strcmp(method_PRED, HW)
        [alpha, beta, pred_ts] = hw_trial(bers_BPSK, bers_BPSK, PRED_GRANULARITY);
        pred_bers_BPSK = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(bers_QPSK, bers_QPSK, PRED_GRANULARITY);
        pred_bers_QPSK = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(bers_16QAM, bers_16QAM, PRED_GRANULARITY);
        pred_bers_16QAM = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(bers_64QAM, bers_64QAM, PRED_GRANULARITY);
        pred_bers_64QAM = pred_ts(1:end-1);

        [alpha, beta, pred_ts] = hw_trial(stdevs_BPSK, stdevs_BPSK, PRED_GRANULARITY);
        pred_stdevs_BPSK = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(stdevs_QPSK, stdevs_QPSK, PRED_GRANULARITY);
        pred_stdevs_QPSK = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(stdevs_16QAM, stdevs_16QAM, PRED_GRANULARITY);
        pred_stdevs_16QAM = pred_ts(1:end-1);
        [alpha, beta, pred_ts] = hw_trial(stdevs_64QAM, stdevs_64QAM, PRED_GRANULARITY);
        pred_stdevs_64QAM = pred_ts(1:end-1);

    else
        error('wrong prediction method');
    end
    pred_bers = [pred_bers_BPSK; pred_bers_QPSK; pred_bers_16QAM; pred_bers_64QAM];
    pred_stdevs = [pred_stdevs_BPSK; pred_stdevs_QPSK; pred_stdevs_16QAM; pred_stdevs_64QAM];



    %% ----------------------------------
    % Rate Adaptation according to method_RA
    %
    % fprintf('- Rate Adaptation according to method_RA\n');
    if strcmp(method_RA, Oracle)
        % tput = cal_oracle_tput(actual_bers, TPUT_TABLE);
        tput = cal_oracle_tput_fec(actual_bers, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);
    elseif strcmp(method_RA, Thresholding)

        %% ----------------------------------
        % Get thresholds
        %

        %% ----------------------------------
        % BER2PRR_conversion
        % threshold = 1;
        % threshold_BER = BER2PRR_conversion(pred_bers, actual_bers, 0.0000001, threshold);
        threshold_BER = BER2PRR_conversion_fec(pred_bers, actual_bers, 0.0000001, threshold, MOD);
        % threshold_BER = [0; 0.000000; 0.000041; 0.002853];
        % threshold_BER = [0.001, 0.01, 0.004, 0.0013];
        % fprintf('threshold: %f\n', threshold_BER);
        

        num_confs = length(MOD);
        num_pkts = size(pred_bers, 2);

        prrs = zeros(num_confs, num_pkts);
        for conf_i = 1:num_confs
            mod_i = MOD(conf_i);
            this_pred_bers = pred_bers(mod_i, :);

            %% not enough samples for BPSK & QPSK, so just accept it
            this_threshold_BER = 1;
            if mod_i > 2
                this_threshold_BER = threshold_BER(conf_i);
            end

            prrs(conf_i, :) = cal_prr_from_threshold(this_pred_bers, this_threshold_BER);
        end

        tput = cal_tput_fec(prrs, actual_bers, threshold, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);
    elseif strcmp(method_RA, Probability_1)
        prrs = cal_probability_fec_1(pred_bers, num_bits_per_pkt, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);

        tput = cal_tput_fec(prrs, actual_bers, threshold, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);
    elseif strcmp(method_RA, Probability_2)
        prrs = cal_probability_fec_2(pred_bers, pred_stdevs, num_bits_per_pkt, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);

        tput = cal_tput_fec(prrs, actual_bers, threshold, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N);
    else
        error('wrong RA method');
    end

    fprintf('trace %d, %s, %s, %s, tput = %f\n', file_ind, method_SNR, method_RA, method_PRED, tput);
end




%% -----------------------------------
% cal_oracle_tput: 
%  calculate throughput using oracle Rate Adaptation scheme
%  
% @input
%   - bers: modulation * num_pkts
%       effective BER of each packet
%   - TPUT_TABLE: modulation * 1
%       throughput gain for each modulation if the packet is correctly received
%
function [tput] = cal_oracle_tput(bers, TPUT_TABLE)
    tput = 0; 
    for pkt_i = 1:size(bers, 2)
        min_ber = min(bers(:, pkt_i));
        min_ber_ind = find(bers(:, pkt_i) == min_ber);
        selected_ind = min_ber_ind(end);

        if min_ber == 0
            tput = tput + TPUT_TABLE(selected_ind);
        end
        % fprintf('pkt %d: bers=(%f, %f, %f, %f), oracle=(%d, %f), tput=%d\n', pkt_i, bers(:, pkt_i), selected_ind, min_ber, tput);
    end
end


function [tput] = cal_oracle_tput_fec(bers, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N)
    num_confs = length(MOD);

    tput = 0;
    for pkt_i = 1:size(bers, 2)
        for conf_i = num_confs:-1:1
            tolerable_err = (RATES_N(conf_i) - RATES_K(conf_i) ) / 2 / RATES_N(conf_i);
            mod_i = MOD(conf_i);
            gain = TPUT_TABLE(conf_i);
            block_size = BLOCK_SIZE;

            this_tput = received_correctly(bers(mod_i, pkt_i), block_size, gain, mod_i, tolerable_err);
            if this_tput > 0
                tput = tput + this_tput;
                break;
            end
        end
    end
end


%% -----------------------------------
% received_correctly: 
%  decide if this packet is correctly received or not
%  XXX: ignore the block for now
%  
% @input
%   - ber: 1 * 1
%       effective BER of the packet
%   - block_size: 1 * 1
%       block_size of this configuration
%   - gain: 1 * 1
%       number of bits got if correctly received
%   - mod_i: 1 * 1
%       which modulation is used for this configuration
%   - tolerabe_err: 1 * 1
%       how many bit errors are torelable for this configuration
%
function [tput] = received_correctly(ber, block_size, gain, mod_i, tolerable_err)
    tput = 0;
    if ber <= tolerable_err
        tput = gain;
    end
end


%% -----------------------------------
% cal_tput: 
%  calculate throughput according PPR:
%  i.e. select the highest configuration which has PPR >= threshold
%  
% @input
%   - prrs: modulation * num_pkts
%       The Packet Reception Rate (PRR) of the trace
%   - actual_bers: modulation * num_pkts
%       actual ber of each frame
%   - threshold
%   - tput_table
%
function [tput] = cal_tput(prrs, actual_bers, threshold, tput_table)
    tput = 0; 
    for pkt_i = 1:size(prrs, 2)
        if prrs(4, pkt_i) >= threshold
            selected_ind = 4;
        elseif prrs(3, pkt_i) >= threshold
            selected_ind = 3;
        elseif prrs(2, pkt_i) >= threshold
            selected_ind = 2;
        else
            selected_ind = 1;
        end

        if actual_bers(selected_ind, pkt_i) == 0
            tput = tput + tput_table(selected_ind);
        end

        % fprintf('pkt %d: bers=(%f, %f, %f, %f), prrs=(%f, %f, %f, %f),  scheme=(%d, %f), tput=%d\n', pkt_i, bers(:, pkt_i), prrs(:, pkt_i), selected_ind, prrs(selected_ind, pkt_i), tput);
    end
end


function [tput] = cal_tput_fec(prrs, actual_bers, threshold, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N)
    num_confs = length(MOD);

    tput = 0; 
    for pkt_i = 1:size(prrs, 2)
        for conf_i = num_confs:-1:1
            if prrs(conf_i, pkt_i) >= threshold | conf_i == 1
                % select this configuration,
                % but need to check if it is actually received correctly

                tolerable_err = (RATES_N(conf_i) - RATES_K(conf_i) ) / 2 / RATES_N(conf_i);
                mod_i = MOD(conf_i);
                gain = TPUT_TABLE(conf_i);
                block_size = BLOCK_SIZE;

                this_tput = received_correctly(actual_bers(mod_i, pkt_i), block_size, gain, mod_i, tolerable_err);
                if this_tput > 0
                    tput = tput + this_tput;
                    break;
                end
            end
        end
    end
end


%% -------------------------------
% cal_prr_from_threshold
%   give the BER of a pkt, determine if the pkt is correctly received by thresholding the BER
%
% @input
%   - bers: 1 * num_pkt
%       the BERs of the whole trace
%   - threshold
%
function [prr] = cal_prr_from_threshold(bers, threshold)
    num_pkts = length(bers);
    prr = zeros(1, num_pkts);

    for pkt_i = 1:num_pkts
        ber = bers(pkt_i);

        if ber > threshold
            prr(1, pkt_i) = 0;
        else
            prr(1, pkt_i) = 1;
        end
    end
end


% @input
%   - bers: 1 * num_pkt
%       the BERs of the whole trace
%   - threshold
%
function [prr] = cal_prr_from_threshold_fec(bers, threshold)
    num_pkts = length(bers);
    prr = zeros(1, num_pkts);

    for pkt_i = 1:num_pkts
        ber = bers(pkt_i);

        if ber > threshold
            prr(1, pkt_i) = 0;
        else
            prr(1, pkt_i) = 1;
        end
    end
end


%% -------------------------------
% cal_eff_ber
%   use the first "num_eff_sym" OFDM symbols to calculate average and stdev SNR of each packet
%
% @input
%   - bers: num_subcarriers * num_ofdm_symbols
%       the BERs of the whole trace
%   - sym_per_pkt
%       number of OFDM symbols per packet
%   - num_eff_sym
%       number of OFDM symbols are used to calculate effective SNR
%
function [eff_bers, stdev_bers] = cal_eff_ber(bers, sym_per_pkt, num_eff_sym)
    [num_sc, num_sym] = size(bers);
    num_pkts = num_sym / sym_per_pkt;

    eff_bers   = zeros(1, num_pkts);
    stdev_bers = zeros(1, num_pkts);
    for pkt_i = 1:num_pkts
        std_ind = (pkt_i - 1) * sym_per_pkt + 1;
        end_ind = std_ind + num_eff_sym - 1;

        this_bers_vector = reshape(bers(:, std_ind:end_ind), [], 1);
        if length(this_bers_vector) ~= num_sc * num_eff_sym
            error('wrong number of symbols used to calculate SNR');
        end

        eff_bers(1, pkt_i)   = mean(this_bers_vector);
        stdev_bers(1, pkt_i) = std(this_bers_vector);
    end
end


%% ---------------------------------
% BER2PRR_conversion
%  i) to plot the figure: BER vs PRR
%  ii) to find the mapping from BER to PRR
%
% @input: 
%   - bers: modulation * num_pkts
%       the BERs of the whole trace
%   - actual_bers: modulation * num_pkts
%       actual ber of each frame
%   - granularity
%       the granularity to seperate BERs
%
function [threshold_BER] = BER2PRR_conversion(bers, actual_bers, granularity, threshold)
    figure_dir = '/u/yichao/mobile_streaming/mobile_streaming/effective_snr/figures_sim3/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';
    max_ber = 0.01;

    [num_mods, num_pkts] = size(bers);
    
    ber_sum = zeros(num_mods, floor(max_ber/granularity + 1));
    ber_cnt = zeros(num_mods, floor(max_ber/granularity + 1));
    for mod_i = 1:num_mods
        for pkt_i = 1:num_pkts
            if(bers(mod_i, pkt_i) > max_ber) 
                continue;
            end

            ber_ind = floor(bers(mod_i, pkt_i) / granularity) + 1;

            ber_sum(mod_i, ber_ind) = ber_sum(mod_i, ber_ind) + actual_bers(mod_i, pkt_i);
            ber_cnt(mod_i, ber_ind) = ber_cnt(mod_i, ber_ind) + 1;
        end
    end

    nonzero_ind = ber_cnt~=0;
    ber_sum(nonzero_ind) = ber_sum(nonzero_ind) ./ ber_cnt(nonzero_ind);
    zero_ind = ber_cnt<=1;
    ber_sum(zero_ind) = inf;


    %% -------------
    % DEBUG
    % tmp_ind = find(ber_cnt(4, :) >= 2);
    % ber_sum(4, tmp_ind)

    
    %% -------------------------
    % plot figure
    % f1 = figure;
    % plot(0:granularity:max_ber, 1-ber_sum(1, :), 'gx', ...
    %      0:granularity:max_ber, 1-ber_sum(2, :), 'y^', ...
    %      0:granularity:max_ber, 1-ber_sum(3, :), 'r*', ...
    %      0:granularity:max_ber, 1-ber_sum(4, :), 'b.');
    % legend('BPSK', 'QPSK', '16QAM', '64QAM');
    for mod_i = 1:num_mods
        fh = figure;
        tmp_ind = find(ber_sum(mod_i, :) <= 1);
        if length(tmp_ind) < 1
            continue;
        end

        plot((tmp_ind-1) * granularity, 1-ber_sum(mod_i, tmp_ind), '.');
        print(fh, '-dpsc', [figure_dir 'BER_2_PRR_conversion_mod' int2str(mod_i) '.eps']);
    end


    %% -------------------------
    % get BER that has 90% PRR
    threshold_BER = zeros(num_mods, 1);
    for mod_i = 1:num_mods
        good_bers = find(ber_sum(mod_i, :) < (1-threshold));

        if length(good_bers) > 0
            threshold_BER(mod_i, 1) = (good_bers(end) - 1) * granularity;
        else
            [val, ind] = min(ber_sum(mod_i, :));
            threshold_BER(mod_i, 1) = (ind(end) - 1) * granularity;
        end
    end


    %% ------------------------
    % DEBUG
    % dlmwrite([output_dir 'tmp_bers.dat'], bers);
    % dlmwrite([output_dir 'tmp_ber2prr.dat'], ber_sum);
    
end


% @input: 
%   - bers: modulation * num_pkts
%       the BERs of the whole trace
%   - actual_bers: modulation * num_pkts
%       actual ber of each frame
%   - granularity
%       the granularity to seperate BERs
function [threshold_BER] = BER2PRR_conversion_fec(bers, actual_bers, granularity, threshold, MOD)
    figure_dir = '/u/yichao/mobile_streaming/mobile_streaming/effective_snr/figures_sim3/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';
    max_ber = 0.01;

    [num_mods, num_pkts] = size(bers);
    num_confs = length(MOD);
    
    ber_sum = zeros(num_confs, floor(max_ber/granularity + 1));
    ber_cnt = zeros(num_confs, floor(max_ber/granularity + 1));
    for conf_i = 1:num_confs
        mod_i = MOD(conf_i);

        for pkt_i = 1:num_pkts
            if(bers(mod_i, pkt_i) > max_ber)
                continue;
            end

            if bers(mod_i, pkt_i) < 0
                bers(mod_i, pkt_i) = 0;
            end

            ber_ind = floor(bers(mod_i, pkt_i) / granularity) + 1;

            ber_sum(conf_i, ber_ind) = ber_sum(conf_i, ber_ind) + actual_bers(mod_i, pkt_i);
            ber_cnt(conf_i, ber_ind) = ber_cnt(conf_i, ber_ind) + 1;
        end
    end

    nonzero_ind = ber_cnt~=0;
    ber_sum(nonzero_ind) = ber_sum(nonzero_ind) ./ ber_cnt(nonzero_ind);
    zero_ind = ber_cnt<=1;
    ber_sum(zero_ind) = inf;


    %% -------------
    % DEBUG
    % tmp_ind = find(ber_cnt(4, :) >= 2);
    % ber_sum(4, tmp_ind)

    
    %% -------------------------
    % plot figure
    % f1 = figure;
    % plot(0:granularity:max_ber, 1-ber_sum(1, :), 'gx', ...
    %      0:granularity:max_ber, 1-ber_sum(2, :), 'y^', ...
    %      0:granularity:max_ber, 1-ber_sum(3, :), 'r*', ...
    %      0:granularity:max_ber, 1-ber_sum(4, :), 'b.');
    % legend('BPSK', 'QPSK', '16QAM', '64QAM');
    for conf_i = 1:num_confs
        fh = figure;
        tmp_ind = find(ber_sum(conf_i, :) <= 1);
        if length(tmp_ind) < 1
            continue;
        end

        plot((tmp_ind-1) * granularity, 1-ber_sum(conf_i, tmp_ind), '.');
        print(fh, '-dpsc', [figure_dir 'BER_2_PRR_conversion_mod' int2str(conf_i) '.eps']);
    end


    %% -------------------------
    % get BER that has 90% PRR
    threshold_BER = zeros(num_confs, 1);
    for conf_i = 1:num_confs
        good_bers = find(ber_sum(conf_i, :) < (1-threshold));

        if length(good_bers) > 0
            threshold_BER(conf_i, 1) = (good_bers(end) - 1) * granularity;
        else
            [val, ind] = min(ber_sum(conf_i, :));
            threshold_BER(conf_i, 1) = (ind(end) - 1) * granularity;
        end
    end


    %% ------------------------
    % DEBUG
    % dlmwrite([output_dir 'tmp_bers.dat'], bers);
    % dlmwrite([output_dir 'tmp_ber2prr.dat'], ber_sum);
    
end




%% --------------------------------
% cal_probability
%    given the mean & avg of BERs, calculate the probability that # of error bit < k
%
% @input:
%   - bers: modulation * num_pkts
%   - stdevs: modulation * num_pkts
%   - num_bits_per_pkt: modulation * 1
%   - k
%
function [prrs] = cal_probability(bers, stdevs, num_bits_per_pkt, k)
    output_dir = '/u/yichao/mobile_streaming/mobile_streaming/effective_snr/OUTPUT_sim_3/';
    [num_mods, num_pkts] = size(bers);

    prrs = zeros(num_mods, num_pkts);
    for mod_i = 1:num_mods
        this_num_bits_per_pkt = num_bits_per_pkt(mod_i, 1);
        this_k = k(mod_i, 1);

        for pkt_i = 1:num_pkts
            prob = 1/2 * (1 + erf((this_k/this_num_bits_per_pkt - bers(mod_i, pkt_i)) / sqrt(2*power(stdevs(mod_i, pkt_i), 2))));
            prrs(mod_i, pkt_i) = prob;
        end
    end


    %% ------------------------
    % DEBUG
    % dlmwrite([output_dir 'tmp_prrs.dat'], prrs);

end


function [prrs] = cal_probability_fec_1(bers, num_bits_per_pkt, TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N)
    output_dir = '/u/yichao/mobile_streaming/mobile_streaming/effective_snr/OUTPUT_sim_3/';
    [num_mods, num_pkts] = size(bers);
    num_confs = length(MOD);

    prrs = zeros(num_confs, num_pkts);
    for conf_i = 1:num_confs
        mod_i = MOD(conf_i);
        this_num_bits_per_pkt = num_bits_per_pkt(mod_i, 1);
        num_err_tolerable = floor(this_num_bits_per_pkt * (RATES_N(conf_i) - RATES_K(conf_i) ) / 2 / RATES_N(conf_i));
        block_size = BLOCK_SIZE;

        for pkt_i = 1:num_pkts
            this_ber = bers(mod_i, pkt_i);
            
            prrs(conf_i, pkt_i) = power(1-this_ber, this_num_bits_per_pkt - num_err_tolerable);
        end
    end


    %% ------------------------
    % DEBUG
    % dlmwrite([output_dir 'tmp_prob_1_prrs.dat'], prrs);

end


function [prrs] = cal_probability_fec_2(bers, stdevs, num_bits_per_pkt, ...
                                        TPUT_TABLE, BLOCK_SIZE, MOD, RATES_K, RATES_N)
    output_dir = '/u/yichao/mobile_streaming/mobile_streaming/effective_snr/OUTPUT_sim_3/';
    [num_mods, num_pkts] = size(bers);
    num_confs = length(MOD);

    prrs = zeros(num_confs, num_pkts);
    for conf_i = 1:num_confs
        mod_i = MOD(conf_i);
        this_num_bits_per_pkt = num_bits_per_pkt(mod_i, 1);
        tolerable_err = (RATES_N(conf_i) - RATES_K(conf_i) ) / 2 / RATES_N(conf_i);
        block_size = BLOCK_SIZE;

        for pkt_i = 1:num_pkts
            this_ber = bers(mod_i, pkt_i);
            this_std = stdevs(mod_i, pkt_i);
            % prob = 1/2 * (1 + erf((tolerable_err - this_ber) / sqrt(2*power(this_std, 2))));
            prob = normcdf(tolerable_err, this_ber, this_std);
            prrs(conf_i, pkt_i) = prob;
        end
    end


    %% ------------------------
    % DEBUG
    % dlmwrite([output_dir 'tmp_prob_2_prrs.dat'], prrs);

end
