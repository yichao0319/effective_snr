function plot_relationship()

    
    %% ----------------------------------
    % global variables
    figure_dir = './figures/';


    %% ----------------------------------
    % SNR to BER
    f1 = figure;
    x_values = -10:25;
    semilogy(x_values, SNR2BER('BPSK', x_values), ...
         x_values, SNR2BER('QPSK', x_values), ...
         x_values, SNR2BER('16QAM', x_values), ...
         x_values, SNR2BER('64QAM', x_values) );
    xlabel('SNR(dB)');
    ylabel('BER');
    axis([-10 25 0.0001 1]);
    legend('BPSK', 'QPSK', '16QAM', '64QAM');
    print(f1, '-dpsc', [figure_dir 'SNR2BER.ps']);



    %% ----------------------------------
    % EVM to SNR
    f2 = figure;
    x_values = 0:0.01:1.8;
    plot(x_values, EVM2SNR(x_values));
    xlabel('EVM');
    ylabel('SNR');
    print(f2, '-dpsc', [figure_dir 'EVM2SNR.ps']);


    %% ----------------------------------
    % SNR to EVM
    f3 = figure;
    % y_values = 0.2:0.01:1.8;
    % x_values = EVM2SNR(y_values);
    x_values = -5:15;
    y_values = SNR2EVM(x_values);
    plot(x_values, y_values);
    xlabel('SNR (dB)');
    ylabel('EVM');
    print(f3, '-dpsc', [figure_dir 'SNR2EVM.ps']);


    %% ----------------------------------
    % EVM to BER
    f2 = figure;
    x_values = 0:0.001:1.8;
    y_values = SNR2BER('QPSK', EVM2SNR(x_values));
    plot(x_values, y_values);
    xlabel('EVM');
    ylabel('BER');
    print(f2, '-dpsc', [figure_dir 'EVM2BER.ps']);



end