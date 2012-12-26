%% -------------------------------
% Yi-Chao Chen @ UT Austin CS
% 
% Input:
%   input_file
%     format:
%       <pkt> <actual_BER                : BPSK, QPSK, 16QAM, 64QAM, mod, pred_mod> 
%             <effSNR BER (formula)      : BPSK, QPSK, 16QAM, 64QAM>
%             <effSNR PRR                : BPSK, QPSK, 16QAM, 64QAM, mod, pred_mod>
%             <effSNR BER (threshold)    : BPSK, QPSK, 16QAM, 64QAM>
%             <effSNR PRR                : BPSK, QPSK, 16QAM, 64QAM, mod, pred_mod>
%             <entire SNR BER (formula)  : BPSK, QPSK, 16QAM, 64QAM>
%             <entire SNR PRR            : BPSK, QPSK, 16QAM, 64QAM, mod, pred_mod>
%             <entire SNR BER (threshold): BPSK, QPSK, 16QAM, 64QAM>
%             <entire SNR PRR            : BPSK, QPSK, 16QAM, 64QAM, mod, pred_mod>
%
% Output:
%
%
% Example:
%   pred_ber('RXDATA.s96-2.pdat')
%

function pred_ber(input_file)
    
    %% ----------------------------------
    % constants
    % CONST_BER_SHIFT_EFFSNR = [0, 0.0016, 0.00078531, 0.003050834];
    % CONST_BER_SHIFT_ENTSNR = [0, 0.0016, 0.000131142, 0.00407008];
    CONST_BER_SHIFT_EFFSNR = [0, 0.0016, 0.000322263, 0.002137679];
    CONST_BER_SHIFT_ENTSNR = [0, 0.0016, -5.60631E-05, 0.003986208];

    DATA_COL = 47;
    NUM_MOD = 4;
    ACL_BER_IND = [2:5];
    EFFSNR_BER_IND = [8:11];
    ENTIRESNR_BER_IND = [28:31];


    %% ----------------------------------
    % global variables
    input_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';
    output_dir = '/v/filer4b/v27q002/ut-wireless/yichao/mobile_streaming/effective_snr/OUTPUT_sim/';

    effSNR_pred_errs = [-1];
    entSNR_pred_errs = [-1];
    effSNR_pred_errs2 = [-1];
    entSNR_pred_errs2 = [-1];


    %% ----------------------------------
    % initinalization
    data = load([input_dir input_file]);
    num_pkt = size(data, 1);
    if DATA_COL ~= size(data, 2)
        error('wrong input size');
    end
    

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
    
    actual_BER = data(2:end, ACL_BER_IND);
    actual_FER = [length(find(actual_BER(:, 1) > 0)) length(find(actual_BER(:, 2) > 0)) length(find(actual_BER(:, 3) > 0)) length(find(actual_BER(:, 4) > 0))] / (num_pkt - 1);
    
    %% ----------------------------------
    %  i) without prediction
    %     a) without BER shift
    effSNR_BER = data(2:end, EFFSNR_BER_IND);
    entSNR_BER = data(2:end, ENTIRESNR_BER_IND);
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     b) with BER shift 1
    effSNR_BER = data(2:end, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt-1, 1);
    entSNR_BER = data(2:end, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt-1, 1);
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     c) with BER shift 2
    effSNR_BER = data(2:end, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt-1, 1);
    negatives_ind = find(effSNR_BER < 0);
    effSNR_BER(negatives_ind) = 0;
    entSNR_BER = data(2:end, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt-1, 1);
    negatives_ind = find(entSNR_BER < 0);
    entSNR_BER(negatives_ind) = 0;
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);



    %% ----------------------------------
    %  ii) prediction 1: previous pkt
    %     a) without BER shift
    effSNR_BER = data(1:end-1, EFFSNR_BER_IND);
    entSNR_BER = data(1:end-1, ENTIRESNR_BER_IND);
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     b) with BER shift 1
    effSNR_BER = data(1:end-1, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt-1, 1);
    entSNR_BER = data(1:end-1, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt-1, 1);
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     c) with BER shift 2
    effSNR_BER = data(1:end-1, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt-1, 1);
    negatives_ind = find(effSNR_BER < 0);
    effSNR_BER(negatives_ind) = 0;
    entSNR_BER = data(1:end-1, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt-1, 1);
    negatives_ind = find(entSNR_BER < 0);
    entSNR_BER(negatives_ind) = 0;
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);



    %% ----------------------------------
    %  iii) prediction 2: EWMA
    granularity = 0.1;
    %     a) without BER shift
    [effSNR_alpha, effSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', data(:, EFFSNR_BER_IND)', granularity);
    [entSNR_alpha, entSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', data(:, ENTIRESNR_BER_IND)', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     b) with BER shift 1
    effSNR_BER = data(:, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt, 1);
    [effSNR_alpha, effSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', effSNR_BER', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = data(:, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt, 1);
    [entSNR_alpha, entSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', entSNR_BER', granularity);
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     c) with BER shift 2
    effSNR_BER = data(:, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt, 1);
    negatives_ind = find(effSNR_BER < 0);
    effSNR_BER(negatives_ind) = 0;
    [effSNR_alpha, effSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', effSNR_BER', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = data(:, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt, 1);
    negatives_ind = find(entSNR_BER < 0);
    entSNR_BER(negatives_ind) = 0;
    [entSNR_alpha, entSNR_BER] = ewma_trial(data(:, ACL_BER_IND)', entSNR_BER', granularity);
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    

    %% ----------------------------------
    %  iv) prediction 3: Holt-Winters
    granularity = 0.1;
    %     a) without BER shift
    [effSNR_alpha, effSNR_beta, effSNR_BER] = hw_trial(data(:, ACL_BER_IND)', data(:, EFFSNR_BER_IND)', granularity);
    [entSNR_alpha, entSNR_beta, entSNR_BER] = hw_trial(data(:, ACL_BER_IND)', data(:, ENTIRESNR_BER_IND)', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     b) with BER shift 1
    effSNR_BER = data(:, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt, 1);
    [effSNR_alpha, effSNR_beta, effSNR_BER] = hw_trial(data(:, ACL_BER_IND)', effSNR_BER', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = data(:, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt, 1);
    [entSNR_alpha, entSNR_beta, entSNR_BER] = hw_trial(data(:, ACL_BER_IND)', entSNR_BER', granularity);
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);
    %     c) with BER shift 2
    effSNR_BER = data(:, EFFSNR_BER_IND) - repmat(CONST_BER_SHIFT_EFFSNR, num_pkt, 1);
    negatives_ind = find(effSNR_BER < 0);
    effSNR_BER(negatives_ind) = 0;
    [effSNR_alpha, effSNR_beta, effSNR_BER] = hw_trial(data(:, ACL_BER_IND)', effSNR_BER', granularity);
    effSNR_BER = effSNR_BER(:, 2:end-1)';
    entSNR_BER = data(:, ENTIRESNR_BER_IND) - repmat(CONST_BER_SHIFT_ENTSNR, num_pkt, 1);
    negatives_ind = find(entSNR_BER < 0);
    entSNR_BER(negatives_ind) = 0;
    [entSNR_alpha, entSNR_beta, entSNR_BER] = hw_trial(data(:, ACL_BER_IND)', entSNR_BER', granularity);
    entSNR_BER = entSNR_BER(:, 2:end-1)';
    effSNR_pred_err = mean(abs(actual_BER - effSNR_BER));
    entSNR_pred_err = mean(abs(actual_BER - entSNR_BER));
    effSNR_pred_errs = store_pred_err(effSNR_pred_errs, effSNR_pred_err);
    entSNR_pred_errs = store_pred_err(entSNR_pred_errs, entSNR_pred_err);
    effSNR_pred_err = abs(mean(actual_BER) - mean(effSNR_BER));
    entSNR_pred_err = abs(mean(actual_BER) - mean(entSNR_BER));
    effSNR_pred_errs2 = store_pred_err(effSNR_pred_errs2, effSNR_pred_err);
    entSNR_pred_errs2 = store_pred_err(entSNR_pred_errs2, entSNR_pred_err);


    %% ----------------
    % output
    %
    mean_actual_BER = mean(actual_BER);
    fprintf('%.10f, %.10f, %.10f, %.10f\n', mean_actual_BER);
    % effSNR_pred_errs([2, 5, 8, 11], 3:4)
    % entSNR_pred_errs([2, 5, 8, 11], 3:4)

    % effSNR_pred_errs2([2, 5, 8, 11], 3:4) ./ repmat(mean_actual_BER(1, 3:4), 4, 1)
    % entSNR_pred_errs2([2, 5, 8, 11], 3:4) ./ repmat(mean_actual_BER(1, 3:4), 4, 1)

    % dlmwrite([output_dir input_file '.effSNR.pred2.dat'], effSNR_pred_errs');
    % dlmwrite([output_dir input_file '.entSNR.pred2.dat'], entSNR_pred_errs');
    out_array = [effSNR_pred_errs([3, 6, 9, 12], 3) entSNR_pred_errs([3, 6, 9, 12], 3) effSNR_pred_errs([3, 6, 9, 12], 4) entSNR_pred_errs([3, 6, 9, 12], 4)]
    dlmwrite([output_dir input_file '.pred.out.dat'], out_array);
    
end



%% store_pred_err: function description
function [store_errs] = store_pred_err(store_errs, new_err)
    if length(store_errs) == 1
        store_errs = new_err;
    else
        store_errs = [store_errs; new_err];
    end
end


