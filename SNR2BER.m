% function [ber] = SNR2BER(mod_type, snr)
%     snr = power(10, snr / 10);    %% dB to power
%     % snr = 80 / 52 * snr;

%     % BPSK
%     if strcmp(mod_type,'BPSK')
%         ber = erfc(power(snr, 0.5) ) / 2;
%     % QPSK modulation
%     elseif strcmp(mod_type,'QPSK')
%         ber = erfc(power(snr / 2, 0.5) ) / 2;
%     % 16-QAM modulation
%     elseif strcmp(mod_type,'16QAM')
%         ber = erfc(power(snr / 10, 0.5) ) * 3 / 8;
%     % 64-QAM modulation
%     elseif strcmp(mod_type,'64QAM')
%         ber = erfc(power(snr / 42, 0.5) ) * 7 / 24;
%     else
%         error('Unimplemented modulation');
%     end
% end


function [ber] = SNR2BER(mod_type, snr)
    snr = power(10, snr / 10);    %% dB to power
    % snr = 80 / 52 * snr;

    % BPSK
    if strcmp(mod_type,'BPSK')
        ber = qfunc(sqrt(2*snr));
    % QPSK modulation
    elseif strcmp(mod_type,'QPSK')
        ber = qfunc(sqrt(snr));
    % 16-QAM modulation
    elseif strcmp(mod_type,'16QAM')
        ber = qfunc(sqrt(snr/5))*3/4;
    % 64-QAM modulation
    elseif strcmp(mod_type,'64QAM')
        ber = qfunc(sqrt(snr/21))*7/12;
    else
        error('Unimplemented modulation');
    end
end


%% qfunc: function description
function q = qfunc(x)
    q = erfc(x/sqrt(2)) / 2;
end

