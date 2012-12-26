% EVM2SNR: EVM = sqrt(1/SNR)
function [snr] = EVM2SNR(evm)
    snr = 1 ./ power(evm, 2);

    snr = 10 * log(snr) / log(10);  % power to dB

    % good_ind = find(evm == 0 | snr > 50);
    % snr(good_ind) = 50;
end
