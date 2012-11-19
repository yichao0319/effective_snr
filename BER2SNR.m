function [snr] = BER2SNR(mod_type, ber)
    % BPSK
    if strcmp(mod_type,'BPSK')
        snr = power(erfcinv(2*ber), 2);
    % QPSK modulation
    elseif strcmp(mod_type,'QPSK')
        snr = 2 * power(erfcinv(2*ber), 2);
    % 16-QAM modulation
    elseif strcmp(mod_type,'16QAM')
        snr = 7 / 3 * power(erfcinv(ber/4), 2);
    % 64-QAM modulation
    elseif strcmp(mod_type,'64QAM')
        snr = power(erfcinv(8/3*ber), 2) * 10;
    else
        error('Unimplemented modulation');
    end
end
