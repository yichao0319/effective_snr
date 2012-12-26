%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%   input_file
%
% Output:
%
%
% Example:
%   sim_rate_adapt2(...)
%

function sim_rate_adapt2(file_ind)
    
    %% ----------------------------------
    % constants
    NUM_MOD = 4;
    NUM_EFF_SYM = 1;
    TPUT_TABLE = [1 2 4 6];
    CONST_BER_THRESHOLD_ALL = [0.001, 0.01, 0.004, 0.0013];
    CONST_BER_THRESHOLD_ALL2 = [0.001, 0.01, 0.004, 0.0025];
    CONST_BER_THRESHOLD_EFF = [0.001, 0.01, 0.005, 0.00025];
    CONST_BER_THRESHOLD_EVM = [0.001, 0.01, 0.0015, 0.0025];

    num_subcarriers = 48;
    num_ofdm_symbol_BPSK = 24;
    num_ofdm_symbol_QPSK = 24;
    num_ofdm_symbol_16QAM = 24;
    num_ofdm_symbol_64QAM = 23;

    num_bits_per_pkt_BPSK = num_subcarriers * num_ofdm_symbol_BPSK * TPUT_TABLE(1);
    num_bits_per_pkt_QPSK = num_subcarriers * num_ofdm_symbol_QPSK * TPUT_TABLE(2);
    num_bits_per_pkt_16QAM = num_subcarriers * num_ofdm_symbol_16QAM * TPUT_TABLE(3);
    num_bits_per_pkt_64QAM = num_subcarriers * num_ofdm_symbol_64QAM * TPUT_TABLE(4);
    

    %% ----------------------------------
    % global variables
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
    actual_bers = load([input_dir actual_ber_file]);
    ber_bpsk_grid = load([input_dir ber_bpsk_file]);
    ber_qpsk_grid = load([input_dir ber_qpsk_file]);
    ber_16qam_grid = load([input_dir ber_16qam_file]);
    ber_64qam_grid = load([input_dir ber_64qam_file]);
    evmber_bpsk_grid = load([input_dir evmber_bpsk_file]);
    evmber_qpsk_grid = load([input_dir evmber_qpsk_file]);
    evmber_16qam_grid = load([input_dir evmber_16qam_file]);
    evmber_64qam_grid = load([input_dir evmber_64qam_file]);

    %% -------------------
    eff_ber_bpsk = cal_eff_ber(ber_bpsk_grid, num_ofdm_symbol_BPSK, NUM_EFF_SYM);
    eff_ber_qpsk = cal_eff_ber(ber_qpsk_grid, num_ofdm_symbol_QPSK, NUM_EFF_SYM);
    eff_ber_16qam = cal_eff_ber(ber_16qam_grid, num_ofdm_symbol_16QAM, NUM_EFF_SYM);
    eff_ber_64qam = cal_eff_ber(ber_64qam_grid, num_ofdm_symbol_64QAM, NUM_EFF_SYM);
    % eff_ber = [eff_ber_bpsk; eff_ber_qpsk; eff_ber_16qam; eff_ber_64qam];

    eff_prr_bpsk = cal_prr(eff_ber_bpsk, num_bits_per_pkt_BPSK, 'BPSK', CONST_BER_THRESHOLD_ALL2);
    eff_prr_qpsk = cal_prr(eff_ber_qpsk, num_bits_per_pkt_QPSK, 'QPSK', CONST_BER_THRESHOLD_ALL2);
    eff_prr_16qam = cal_prr(eff_ber_16qam, num_bits_per_pkt_16QAM, '16QAM', CONST_BER_THRESHOLD_ALL2);
    eff_prr_64qam = cal_prr(eff_ber_64qam, num_bits_per_pkt_64QAM, '64QAM', CONST_BER_THRESHOLD_ALL2);
    eff_prrs = [eff_prr_bpsk; eff_prr_qpsk; eff_prr_16qam; eff_prr_64qam];

    %% -------------------
    ent_ber_bpsk = cal_eff_ber(ber_bpsk_grid, num_ofdm_symbol_BPSK, num_ofdm_symbol_BPSK);
    ent_ber_qpsk = cal_eff_ber(ber_qpsk_grid, num_ofdm_symbol_QPSK, num_ofdm_symbol_QPSK);
    ent_ber_16qam = cal_eff_ber(ber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
    ent_ber_64qam = cal_eff_ber(ber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    % ent_ber = [ent_ber_bpsk; ent_ber_qpsk; ent_ber_16qam; ent_ber_64qam];

    ent_prr_bpsk = cal_prr(ent_ber_bpsk, num_bits_per_pkt_BPSK, 'BPSK', CONST_BER_THRESHOLD_ALL);
    ent_prr_qpsk = cal_prr(ent_ber_qpsk, num_bits_per_pkt_QPSK, 'QPSK', CONST_BER_THRESHOLD_ALL);
    ent_prr_16qam = cal_prr(ent_ber_16qam, num_bits_per_pkt_16QAM, '16QAM', CONST_BER_THRESHOLD_ALL);
    ent_prr_64qam = cal_prr(ent_ber_64qam, num_bits_per_pkt_64QAM, '64QAM', CONST_BER_THRESHOLD_ALL);
    ent_prrs = [ent_prr_bpsk; ent_prr_qpsk; ent_prr_16qam; ent_prr_64qam];

    %% -------------------
    ent_evmber_bpsk  = cal_eff_ber(evmber_bpsk_grid,  num_ofdm_symbol_BPSK,  num_ofdm_symbol_BPSK);
    ent_evmber_qpsk  = cal_eff_ber(evmber_qpsk_grid,  num_ofdm_symbol_QPSK,  num_ofdm_symbol_QPSK);
    ent_evmber_16qam = cal_eff_ber(evmber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
    ent_evmber_64qam = cal_eff_ber(evmber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    % ent_ber = [ent_ber_bpsk; ent_ber_qpsk; ent_ber_16qam; ent_ber_64qam];

    ent_evmprr_bpsk  = cal_prr(ent_evmber_bpsk,  num_bits_per_pkt_BPSK, 'BPSK', CONST_BER_THRESHOLD_ALL);
    ent_evmprr_qpsk  = cal_prr(ent_evmber_qpsk,  num_bits_per_pkt_QPSK, 'QPSK', CONST_BER_THRESHOLD_ALL);
    ent_evmprr_16qam = cal_prr(ent_evmber_16qam, num_bits_per_pkt_16QAM, '16QAM', CONST_BER_THRESHOLD_ALL);
    ent_evmprr_64qam = cal_prr(ent_evmber_64qam, num_bits_per_pkt_64QAM, '64QAM', CONST_BER_THRESHOLD_ALL);
    ent_evmprrs = [ent_evmprr_bpsk; ent_evmprr_qpsk; ent_evmprr_16qam; ent_evmprr_64qam];


    %% ----------------------------------
    % main
    %  i) without prediction
    %     a) without BER shift
    %     b) with BER shift 1
    %     c) with BER shift 2
    %  ii) prediction 1: previous pkt
    %     a) ...
    %  iii) prediction 2: EWMA
    %     a) ...
    %  iv) prediction 3: Holt-Winters
    %     a) ...



    %% ----------------------------------
    %  i) without prediction
    oracle_tput = cal_actual_tput(actual_bers, TPUT_TABLE);
    eff_tput = cal_tput(eff_prrs, actual_bers, TPUT_TABLE);
    ent_tput = cal_tput(ent_prrs, actual_bers, TPUT_TABLE);
    fprintf('oracle=%d, EffSNR=%d(%f), EntireSNR=%d(%f)\n', oracle_tput, eff_tput, eff_tput/oracle_tput, ent_tput, ent_tput/oracle_tput);

    %% ----------------------------------
    %  ii) prediction 1: previous pkt
    this_actual_bers = actual_bers(:, 2:end);
    this_eff_prrs = eff_prrs(:, 1:end-1);
    this_ent_prrs = ent_prrs(:, 1:end-1);
    
    oracle_tput = cal_actual_tput(this_actual_bers, TPUT_TABLE);
    eff_tput = cal_tput(this_eff_prrs, this_actual_bers, TPUT_TABLE);
    ent_tput = cal_tput(this_ent_prrs, this_actual_bers, TPUT_TABLE);
    fprintf('oracle=%d, EffSNR=%d(%f), EntireSNR=%d(%f)\n', oracle_tput, eff_tput, eff_tput/oracle_tput, ent_tput, ent_tput/oracle_tput);


    %% ----------------------------------
    %  iii) without prediction, evm 
    oracle_tput = cal_actual_tput(actual_bers, TPUT_TABLE);
    eff_tput = cal_tput(eff_prrs, actual_bers, TPUT_TABLE);
    ent_tput = cal_tput(ent_evmprrs, actual_bers, TPUT_TABLE);
    fprintf('oracle=%d, EffSNR=%d(%f), EntireSNR=%d(%f)\n', oracle_tput, eff_tput, eff_tput/oracle_tput, ent_tput, ent_tput/oracle_tput);

    %% ----------------------------------
    %  iv) prediction 1: previous pkt, evm
    this_actual_bers = actual_bers(:, 2:end);
    this_eff_prrs = eff_prrs(:, 1:end-1);
    this_ent_prrs = ent_evmprrs(:, 1:end-1);
    
    oracle_tput = cal_actual_tput(this_actual_bers, TPUT_TABLE);
    eff_tput = cal_tput(this_eff_prrs, this_actual_bers, TPUT_TABLE);
    ent_tput = cal_tput(this_ent_prrs, this_actual_bers, TPUT_TABLE);
    fprintf('oracle=%d, EffSNR=%d(%f), EntireSNR=%d(%f)\n\n', oracle_tput, eff_tput, eff_tput/oracle_tput, ent_tput, ent_tput/oracle_tput);
    
end




%% cal_actual_tput: function description
function [tput] = cal_actual_tput(bers, tput_table)
    % ber: modulation * pkt
    tput = 0; 
    for pkt_i = 1:size(bers, 2)
        min_ber = min(bers(:, pkt_i));
        min_ber_ind = find(bers(:, pkt_i) == min_ber);
        selected_ind = min_ber_ind(end);

        if min_ber == 0
            tput = tput + tput_table(selected_ind);
        end
        % fprintf('pkt %d: bers=(%f, %f, %f, %f), oracle=(%d, %f), tput=%d\n', pkt_i, bers(:, pkt_i), selected_ind, min_ber, tput);
    end
end


%% cal_tput: function description
function [tput] = cal_tput(prrs, bers, tput_table)
    tput = 0; 
    for pkt_i = 1:size(bers, 2)
        if prrs(4, pkt_i) >= 0.9
            selected_ind = 4;
        elseif prrs(3, pkt_i) >= 0.9
            selected_ind = 3;
        elseif prrs(2, pkt_i) >= 0.9
            selected_ind = 2;
        else
            selected_ind = 1;
        end

        if bers(selected_ind, pkt_i) == 0
            tput = tput + tput_table(selected_ind);
        end

        % fprintf('pkt %d: bers=(%f, %f, %f, %f), prrs=(%f, %f, %f, %f),  scheme=(%d, %f), tput=%d\n', pkt_i, bers(:, pkt_i), prrs(:, pkt_i), selected_ind, prrs(selected_ind, pkt_i), tput);
    end
end



%% cal_eff_ber:bers, num_eff_sym description
function [eff_bers] = cal_eff_ber(bers, sym_per_pkt, num_eff_sym)
    [num_sc, num_sym] = size(bers);
    num_pkts = num_sym / sym_per_pkt;

    eff_bers = zeros(1, num_pkts);
    for pkt_i = 1:num_pkts
        std_ind = (pkt_i - 1) * sym_per_pkt + 1;
        end_ind = std_ind + num_eff_sym - 1;
        eff_bers(1, pkt_i) = mean(mean(bers(:, std_ind:end_ind)));
    end
end



%% cal_prr: function description
% function [prr] = cal_prr(bers, num_bits, modulation)
%     num_pkts = length(bers);
%     prr = zeros(1, num_pkts);

%     for pkt_i = 1:num_pkts
%         ber = bers(pkt_i);

%         if strcmp(modulation, 'BPSK')
%             prr(1, pkt_i) = power(1-ber, num_bits);
%         elseif strcmp(modulation, 'QPSK')
%             if ber > 0.001
%                 prr(1, pkt_i) = 0;
%             else
%                 prr(1, pkt_i) = 1;
%             end
%         elseif strcmp(modulation, '16QAM')
%             if ber > 0.000207
%                 prr(1, pkt_i) = 0;
%             else
%                 prr(1, pkt_i) = 1;
%             end
%         elseif strcmp(modulation, '64QAM')
%             % if ber > 0.0026631
%             if ber > 0.0025
%                 prr(1, pkt_i) = 0;
%             else
%                 prr(1, pkt_i) = 1;
%             end
%         end
%     end
% end


function [prr] = cal_prr(bers, num_bits, modulation, threshold)
    num_pkts = length(bers);
    prr = zeros(1, num_pkts);

    for pkt_i = 1:num_pkts
        ber = bers(pkt_i);

        if strcmp(modulation, 'BPSK')
            prr(1, pkt_i) = power(1-ber, num_bits);
        elseif strcmp(modulation, 'QPSK')
            if ber > threshold(2)
                prr(1, pkt_i) = 0;
            else
                prr(1, pkt_i) = 1;
            end
        elseif strcmp(modulation, '16QAM')
            if ber > threshold(3)
                prr(1, pkt_i) = 0;
            else
                prr(1, pkt_i) = 1;
            end
        elseif strcmp(modulation, '64QAM')
            % if ber > 0.0026631
            if ber > threshold(4)
                prr(1, pkt_i) = 0;
            else
                prr(1, pkt_i) = 1;
            end
        end
    end
end


