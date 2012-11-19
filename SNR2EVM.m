% SNR2EVM: EVM = sqrt(1/SNR)
function [evm] = SNR2EVM(snr)
    snr = power(10, snr / 10);    %% dB to power

    evm = abs(sqrt(1 ./ snr));
end
