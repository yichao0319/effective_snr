% batch file to run analyze_csi
% e.g.:
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
%   method_FEC
%     [NO_FEC | FEC_RS | FEC_DIVERSITY_GAIN]
%   method_estimate_frame_err: the way to estimate frame error rate
%     [FROM_BER | GILBERT]
%   num_ofdm_symbol_per_pkt:
%   modulation:
%     [BPSK | QPSK | 16QAM | 64QAM]
%   code_rate_k: message length
%   code_rate_n: block length
%   alpha: the parameter for EWMA


function batch_analyze_csi()

    % input_sym_files = {'RXDATA.static.d0.pdat'};
    input_sym_files = {'XXXXX'};
    % input_sym_files = {'RXDATA.static.d0.pdat', 'RXDATA.static.d6.pdat', 'RXDATA.mobile.2.pdat', 'RXDATA.mobile.b1.pdat', 'RXDATA.mobile.b2.pdat', 'RXDATA.s96.pdat', 'RXDATA.s96-2.pdat'};
    % input_sym_files = {'RXDATA.static.highSNR.pdat', 'RXDATA.static.highSNR.2.pdat', 'RXDATA.static.highSNR.d1.pdat'};
    % input_sym_files = {'RXDATA.static.highSNR.d1.2.pdat', 'RXDATA.static.highSNR.d2.pdat'};
    % input_sym_files = {'RXDATA.static.highSNR.d6.pdat', 'RXDATA.static.highSNR.d6.2.pdat', 'RXDATA.mobile.highSNR.pdat'};
    % input_sym_files = {'RXDATA.mobile.highSNR.1.pdat', 'RXDATA.mobile.highSNR.b1.pdat', 'RXDATA.mobile.highSNR.b2.pdat'};
    
    methods_prediction = {'ORACLE', 'EWMA'};
    methods_get_ber = {'ACTUAL_SNR', 'GUESSED_EVM', 'ACTUAL_EVM', 'AVERAGE_SNR'};
    methods_pkt_ber = {'PREAMBLE', 'ENTIRE'};
    methods_snr2ber = {'FORMULA', 'THRESHOLD'};
    methods_fec = {'NO_FEC', 'RS', 'DIVERSITY_GAIN'}; % 
    methods_estimate_frame_err = {'FROM_BER'};
    num_codes = 1;
    codes_k = [1];
    codes_n = [2];


    avg_actual_bers = zeros(length(input_sym_files), length(methods_snr2ber), length(methods_prediction), length(methods_fec), length(methods_estimate_frame_err), num_codes, length(methods_get_ber), length(methods_pkt_ber) );
    avg_predit_bers = zeros(length(input_sym_files), length(methods_snr2ber), length(methods_prediction), length(methods_fec), length(methods_estimate_frame_err), num_codes, length(methods_get_ber), length(methods_pkt_ber) );

    actual_frame_err_rates = zeros(length(input_sym_files), length(methods_snr2ber), length(methods_prediction), length(methods_fec), length(methods_estimate_frame_err), num_codes, length(methods_get_ber), length(methods_pkt_ber) );
    estimated_frame_err_rates = zeros(length(input_sym_files), length(methods_snr2ber), length(methods_prediction), length(methods_fec), length(methods_estimate_frame_err), num_codes, length(methods_get_ber), length(methods_pkt_ber) );


    output_dir = './OUTPUT/';
    fid = fopen([output_dir 'batch_out.txt'], 'w');

    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));

        for snr2ber_i = 1:length(methods_snr2ber)
            method_snr2ber = char(methods_snr2ber(snr2ber_i) );

            for pdct_i = 1:length(methods_prediction)
                method_prediction = char(methods_prediction(pdct_i));

                for fec_i = 1:length(methods_fec)
                    method_fec = char(methods_fec(fec_i));

                    for frame_err_i = 1:length(methods_estimate_frame_err)
                        method_estimate_frame_err = char(methods_estimate_frame_err(frame_err_i));

                        for code_i = 1:num_codes

                            % for getber_i = 1:length(methods_get_ber)
                            %     method_get_ber = char(methods_get_ber(getber_i));

                            %     for pktber_i = 1:length(methods_pkt_ber)
                            %         method_pkt_ber = char(methods_pkt_ber(pktber_i));

                            %         [avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate] = analyze_csi(input_sym_file, method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, method_fec, method_estimate_frame_err, 50, 'QPSK', codes_k(code_i), codes_n(code_i), 1);

                            %         fprintf('%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);
                            %         fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);

                            %         avg_actual_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_actual_ber;
                            %         avg_predit_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_predit_ber;
                            %         actual_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = actual_frame_err_rate;
                            %         estimated_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = estimated_frame_err_rate;

                            %     end
                            % end


                            %% ---------------------
                            % 2. ACTUAL_SNR + PREAMBLE
                            getber_i = 1;
                            method_get_ber = 'ACTUAL_SNR';
                            pktber_i = 1;
                            method_pkt_ber = 'PREAMBLE';

                            [avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate] = analyze_csi(input_sym_file, method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, method_fec, method_estimate_frame_err, 50, 'QPSK', codes_k(code_i), codes_n(code_i), 1);
                            fprintf('%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);
                            fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);

                            avg_actual_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_actual_ber;
                            avg_predit_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_predit_ber;
                            actual_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = actual_frame_err_rate;
                            estimated_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = estimated_frame_err_rate;




                            %% ---------------------
                            % 3. ACTUAL_SNR + ENTIRE
                            getber_i = 1;
                            method_get_ber = 'ACTUAL_SNR';
                            pktber_i = 2;
                            method_pkt_ber = 'ENTIRE';

                            [avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate] = analyze_csi(input_sym_file, method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, method_fec, method_estimate_frame_err, 50, 'QPSK', codes_k(code_i), codes_n(code_i), 1);
                            fprintf('%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);
                            fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);

                            avg_actual_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_actual_ber;
                            avg_predit_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_predit_ber;
                            actual_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = actual_frame_err_rate;
                            estimated_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = estimated_frame_err_rate;


                            %% ---------------------
                            % 4. GUESSED_EVM + ENTIRE
                            getber_i = 2;
                            method_get_ber = 'GUESSED_EVM';
                            pktber_i = 2;
                            method_pkt_ber = 'ENTIRE';

                            [avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate] = analyze_csi(input_sym_file, method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, method_fec, method_estimate_frame_err, 50, 'QPSK', codes_k(code_i), codes_n(code_i), 1);
                            fprintf('%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);
                            fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);

                            avg_actual_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_actual_ber;
                            avg_predit_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_predit_ber;
                            actual_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = actual_frame_err_rate;
                            estimated_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = estimated_frame_err_rate;


                            %% ---------------------
                            % 1. AVERAGE_SNR + PREAMBLE
                            getber_i = 4;
                            method_get_ber = 'AVERAGE_SNR';
                            pktber_i = 1;
                            method_pkt_ber = 'PREAMBLE';

                            [avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate] = analyze_csi(input_sym_file, method_prediction, method_get_ber, method_pkt_ber, method_snr2ber, method_fec, method_estimate_frame_err, 50, 'QPSK', codes_k(code_i), codes_n(code_i), 1);
                            fprintf('%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);
                            fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %d/%d: %f, %f, %f, %f\n', input_sym_file, method_snr2ber, method_prediction, method_fec, method_estimate_frame_err, method_get_ber, method_pkt_ber, codes_k(code_i), codes_n(code_i), avg_actual_ber, avg_predit_ber, actual_frame_err_rate, estimated_frame_err_rate);

                            avg_actual_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_actual_ber;
                            avg_predit_bers(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = avg_predit_ber;
                            actual_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = actual_frame_err_rate;
                            estimated_frame_err_rates(file_i, snr2ber_i, pdct_i, fec_i, frame_err_i, code_i, getber_i, pktber_i) = estimated_frame_err_rate;

                        end
                    end
                end
            end
        end
    end

    fclose(fid);



    %% -------------------------------------
    % input format:
    %    <filename>, <method_snr2ber             [FORMULA | THRESHOLD]>, 
    %                <method_prediction          [ORACLE | EWMA]>, 
    %                <method_fec                 [NO_FEC | FEC_RS | FEC_DIVERSITY_GAIN]>
    %                <method_estimate_frame_err  [FROM_BER | GILBERT]>
    %                <fec_codes>
    %                <method_get_ber             [ACTUAL_SNR | GUESSED_EVM | ACTUAL_EVM | AVERAGE_SNR]>
    %                <method_pkt_ber             [PREAMBLE | ENTIRE]>
    fid2 = fopen([output_dir 'batch_err_ratio.txt'], 'w');
    
    %% -------------------------------------
    %  BER
    ber_error = abs(avg_actual_bers - avg_predit_bers);
    ber_error_ratio = ber_error ./ avg_actual_bers;
    % output format:
    %    <filename>, <avg pkt snr:         [AVERAGE_SNR + PREAMBLE]>, 
    %                <preamble-SNR:        [ACTUAL_SNR + PREAMBLE]>, 
    %                <actual entire-SNR:   [ACTUAL_SNR + ENTIRE]>, 
    %                <entire-SNR from EVM: [GUESSED_EVM + ENTIRE]>

    % SNR2BER formula, oracle prediction:  [FORMULA + ORACLE]
    fprintf(fid2, 'SNR2BER formula, oracle prediction\n');
    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));
        fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ber_error_ratio(file_i, 1, 1, 1, 1, 1, 4, 1), ...
                                                              ber_error_ratio(file_i, 1, 1, 1, 1, 1, 1, 1), ...
                                                              ber_error_ratio(file_i, 1, 1, 1, 1, 1, 1, 2), ...
                                                              ber_error_ratio(file_i, 1, 1, 1, 1, 1, 2, 2) );
    end

    % SNR2BER formula, EWMA prediction:  [FORMULA + EWMA]
    fprintf(fid2, '\nSNR2BER formula, EWMA prediction\n');
    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));
        fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ber_error_ratio(file_i, 1, 2, 1, 1, 1, 4, 1), ...
                                                              ber_error_ratio(file_i, 1, 2, 1, 1, 1, 1, 1), ...
                                                              ber_error_ratio(file_i, 1, 2, 1, 1, 1, 1, 2), ...
                                                              ber_error_ratio(file_i, 1, 2, 1, 1, 1, 2, 2) );
    end

    % SNR2BER threshold, oracle prediction: [THRESHOLD + ORACLE]
    fprintf(fid2, '\nSNR2BER threshold, oracle prediction\n');
    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));
        fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ber_error_ratio(file_i, 2, 1, 1, 1, 1, 4, 1), ...
                                                              ber_error_ratio(file_i, 2, 1, 1, 1, 1, 1, 1), ...
                                                              ber_error_ratio(file_i, 2, 1, 1, 1, 1, 1, 2), ...
                                                              ber_error_ratio(file_i, 2, 1, 1, 1, 1, 2, 2) );
    end

    % SNR2BER threshold, EWMA prediction:  [FORMULA + EWMA]
    fprintf(fid2, '\nSNR2BER threshold, EWMA prediction\n');
    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));
        fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ber_error_ratio(file_i, 2, 2, 1, 1, 1, 4, 1), ...
                                                              ber_error_ratio(file_i, 2, 2, 1, 1, 1, 1, 1), ...
                                                              ber_error_ratio(file_i, 2, 2, 1, 1, 1, 1, 2), ...
                                                              ber_error_ratio(file_i, 2, 2, 1, 1, 1, 2, 2) );
    end


    %% -------------------------------------
    %  frame error rate
    frame_error = abs(actual_frame_err_rates - estimated_frame_err_rates);
    fprintf(fid2, '\n\n\nframe error rate\n');
    % output format:
    %    <filename>, <avg pkt snr:         [AVERAGE_SNR + PREAMBLE]>, 
    %                <preamble-SNR:        [ACTUAL_SNR + PREAMBLE]>, 
    %                <actual entire-SNR:   [ACTUAL_SNR + ENTIRE]>, 
    %                <entire-SNR from EVM: [GUESSED_EVM + ENTIRE]>

    for fec_i = 1:length(methods_fec)
        method_fec = char(methods_fec(fec_i));

        for frame_err_i = 1:length(methods_estimate_frame_err)
            method_estimate_frame_err = char(methods_estimate_frame_err(frame_err_i));

            for code_i = 1:num_codes
                fprintf(fid2, '\n%s, %s, code rate %d:\n', method_fec, method_estimate_frame_err, code_i);
                % SNR2BER formula, oracle prediction:  [FORMULA + ORACLE]
                fprintf(fid2, 'SNR2BER formula, oracle prediction\n');
                for file_i = 1:length(input_sym_files)
                    input_sym_file = char(input_sym_files(file_i));
                    fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ...
                                    frame_error(file_i, 1, 1, fec_i, frame_err_i, code_i, 4, 1), ...
                                    frame_error(file_i, 1, 1, fec_i, frame_err_i, code_i, 1, 1), ...
                                    frame_error(file_i, 1, 1, fec_i, frame_err_i, code_i, 1, 2), ...
                                    frame_error(file_i, 1, 1, fec_i, frame_err_i, code_i, 2, 2) );
                end

                % SNR2BER formula, EWMA prediction:  [FORMULA + EWMA]
                fprintf(fid2, '\nSNR2BER formula, EWMA prediction\n');
                for file_i = 1:length(input_sym_files)
                    input_sym_file = char(input_sym_files(file_i));
                    fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ...
                                    frame_error(file_i, 1, 2, fec_i, frame_err_i, code_i, 4, 1), ...
                                    frame_error(file_i, 1, 2, fec_i, frame_err_i, code_i, 1, 1), ...
                                    frame_error(file_i, 1, 2, fec_i, frame_err_i, code_i, 1, 2), ...
                                    frame_error(file_i, 1, 2, fec_i, frame_err_i, code_i, 2, 2) );
                end

                % SNR2BER threshold, oracle prediction: [THRESHOLD + ORACLE]
                fprintf(fid2, '\nSNR2BER threshold, oracle prediction\n');
                for file_i = 1:length(input_sym_files)
                    input_sym_file = char(input_sym_files(file_i));
                    fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ...
                                    frame_error(file_i, 2, 1, fec_i, frame_err_i, code_i, 4, 1), ...
                                    frame_error(file_i, 2, 1, fec_i, frame_err_i, code_i, 1, 1), ...
                                    frame_error(file_i, 2, 1, fec_i, frame_err_i, code_i, 1, 2), ...
                                    frame_error(file_i, 2, 1, fec_i, frame_err_i, code_i, 2, 2) );
                end

                % SNR2BER threshold, EWMA prediction:  [FORMULA + EWMA]
                fprintf(fid2, '\nSNR2BER threshold, EWMA prediction\n');
                for file_i = 1:length(input_sym_files)
                    input_sym_file = char(input_sym_files(file_i));
                    fprintf(fid2, '%s, %f, %f, %f, %f\n', input_sym_file, ...
                                    frame_error(file_i, 2, 2, fec_i, frame_err_i, code_i, 4, 1), ...
                                    frame_error(file_i, 2, 2, fec_i, frame_err_i, code_i, 1, 1), ...
                                    frame_error(file_i, 2, 2, fec_i, frame_err_i, code_i, 1, 2), ...
                                    frame_error(file_i, 2, 2, fec_i, frame_err_i, code_i, 2, 2) );
                end
            end
        end
    end

    fclose(fid2);



end