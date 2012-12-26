%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%
% Output:
%
%
% Example:
%   pred_fer('rx_actual_bers_run1.dat', 'rx_bers_bpsk_run1.dat', 'rx_bers_qpsk_run1.dat', 'rx_bers_16qam_run1.dat', 'rx_bers_64qam_run1.dat')
%

function pred_fer(file_ind)
    
    %% ----------------------------------
    % constants
    NUM_MOD = 4;
    NUM_EFF_SYM = 1;
    TPUT_TABLE = [1 2 4 6];

    granularity = 0.000001;
    max_ber = 0.1;


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
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_pred_per/';
    figure_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/figures_pred_per/';
    actual_ber_file = ['rx_actual_bers_run' int2str(file_ind) '.dat'];
    
    ber_bpsk_file   = ['rx_bers_bpsk_run' int2str(file_ind) '.dat'];
    ber_qpsk_file   = ['rx_bers_qpsk_run' int2str(file_ind) '.dat'];
    ber_16qam_file  = ['rx_bers_16qam_run' int2str(file_ind) '.dat'];
    ber_64qam_file  = ['rx_bers_64qam_run' int2str(file_ind) '.dat'];

    evmber_bpsk_file   = ['rx_evmbers_bpsk_run' int2str(file_ind) '.dat'];
    evmber_qpsk_file   = ['rx_evmbers_qpsk_run' int2str(file_ind) '.dat'];
    evmber_16qam_file  = ['rx_evmbers_16qam_run' int2str(file_ind) '.dat'];
    evmber_64qam_file  = ['rx_evmbers_64qam_run' int2str(file_ind) '.dat'];

    effber_bad = zeros(NUM_MOD, floor(max_ber / granularity) + 1);
    effber_good = zeros(NUM_MOD, floor(max_ber / granularity) + 1);
    entber_bad = zeros(NUM_MOD, floor(max_ber / granularity) + 1);
    entber_good = zeros(NUM_MOD, floor(max_ber / granularity) + 1);
    evmentber_bad = zeros(NUM_MOD, floor(max_ber / granularity) + 1);
    evmentber_good = zeros(NUM_MOD, floor(max_ber / granularity) + 1);


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
    eff_bers = [eff_ber_bpsk; eff_ber_qpsk; eff_ber_16qam; eff_ber_64qam];
    [effber_bad, effber_good] = cal_distribution(eff_bers, actual_bers, max_ber, granularity);

    for mod_i = 1:NUM_MOD
        range = 1:floor(max_ber/granularity)+1;
        if mod_i == 1
            range = 1:5;
        elseif mod_i == 2
            range = 1:5;
        elseif mod_i == 3
            range = 1:250;
        elseif mod_i == 4
            range = 1:40000;
        end
                
        f1 = figure;
        plot(range * granularity, effber_bad(mod_i, range), range * granularity, effber_good(mod_i, range));
        legend('bad', 'good')
        print(f1, '-dpsc', [figure_dir 'rx_eff_mod' int2str(mod_i) '_run' int2str(file_ind) '.eps']);

        % fprintf('%f: %f\n', [range; effber_bad(mod_i, range)]);
    end
    

    %% -------------------
    ent_ber_bpsk = cal_eff_ber(ber_bpsk_grid, num_ofdm_symbol_BPSK, num_ofdm_symbol_BPSK);
    ent_ber_qpsk = cal_eff_ber(ber_qpsk_grid, num_ofdm_symbol_QPSK, num_ofdm_symbol_QPSK);
    ent_ber_16qam = cal_eff_ber(ber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
    ent_ber_64qam = cal_eff_ber(ber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    ent_bers = [ent_ber_bpsk; ent_ber_qpsk; ent_ber_16qam; ent_ber_64qam];
    [entber_bad, entber_good] = cal_distribution(ent_bers, actual_bers, max_ber, granularity);

    for mod_i = 1:NUM_MOD
        range = 1:floor(max_ber/granularity)+1;
        if mod_i == 1
            range = 1:5;
        elseif mod_i == 2
            range = 1:10;
        elseif mod_i == 3
            range = 1:30;
        % elseif mod_i == 4
        %    range = 1:300;
        end

        f1 = figure;
        plot(range * granularity, entber_bad(mod_i, range), range * granularity, entber_good(mod_i, range));
        legend('bad', 'good')
        print(f1, '-dpsc', [figure_dir 'rx_ent_mod' int2str(mod_i) '_run' int2str(file_ind) '.eps']);
    end

    %% -------------------
    ent_evmber_bpsk  = cal_eff_ber(evmber_bpsk_grid,  num_ofdm_symbol_BPSK,  num_ofdm_symbol_BPSK);
    ent_evmber_qpsk  = cal_eff_ber(evmber_qpsk_grid,  num_ofdm_symbol_QPSK,  num_ofdm_symbol_QPSK);
    ent_evmber_16qam = cal_eff_ber(evmber_16qam_grid, num_ofdm_symbol_16QAM, num_ofdm_symbol_16QAM);
    ent_evmber_64qam = cal_eff_ber(evmber_64qam_grid, num_ofdm_symbol_64QAM, num_ofdm_symbol_64QAM);
    ent_evmbers = [ent_evmber_bpsk; ent_evmber_qpsk; ent_evmber_16qam; ent_evmber_64qam];
    [entevmber_bad, entevmber_good] = cal_distribution(ent_evmbers, actual_bers, max_ber, granularity);

    for mod_i = 1:NUM_MOD
        range = 1:floor(max_ber/granularity)+1;
        % if mod_i == 1
        %     range = 1:5;
        % elseif mod_i == 2
        %     range = 1:5;
        % elseif mod_i == 3
        %     range = 1:5000;
        % elseif mod_i == 4
        %    range = 1:10000;
        % end

        f1 = figure;
        plot(range * granularity, entevmber_bad(mod_i, range), range * granularity, entevmber_good(mod_i, range));
        legend('bad', 'good')
        print(f1, '-dpsc', [figure_dir 'rx_evment_mod' int2str(mod_i) '_run' int2str(file_ind) '.eps']);
    end
    return;

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
    fprintf('oracle=%d, EffSNR=%d(%f), EntireSNR=%d(%f)\n', oracle_tput, eff_tput, eff_tput/oracle_tput, ent_tput, ent_tput/oracle_tput);
    
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



%% cal_distribution: function description
function [bad, good] = cal_distribution(bers, actual_bers, max_ber, granularity)
    [num_mod, num_pkts] = size(bers);
    if num_mod ~= 4
        error('wrong input');
    end

    bad = zeros(num_mod, floor(max_ber / granularity) + 1);
    good = zeros(num_mod, floor(max_ber / granularity) + 1);

    for pkt_i = 1:num_pkts
        for mod_i = 1:num_mod
            ber = bers(mod_i, pkt_i);
            if ber > max_ber
                fprintf('ber > max_ber\n');
                continue;
            end

            ind = floor(ber / granularity) + 1;

            actual_ber = actual_bers(mod_i, pkt_i);
            if actual_ber > 0
                %% bad
                bad(mod_i, ind) = bad(mod_i, ind) + 1;
            else
                %% good
                good(mod_i, ind) = good(mod_i, ind) + 1;
            end
        end
    end
end
