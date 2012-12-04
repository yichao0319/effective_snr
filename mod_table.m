%% mod_table: function description
function [table, m, table2] = mod_table(mod_type)
    % BPSK
    if strcmp(mod_type,'BPSK')
        m = 1;
        table=exp(j*[0 -pi]);  % generates BPSK symbols
        table2 = [1 0];
    % QPSK modulation
    elseif strcmp(mod_type,'QPSK')
        m = 2;
        table=exp(j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
        table2 = [0 1 3 2];
    % 16-QAM modulation
    elseif strcmp(mod_type,'16QAM')
        % generates 16QAM symbols
        m=1;
        for k=-3:2:3
            for l=-3:2:3
                table(m) = (k+j*l)/sqrt(10); % power normalization
                m=m+1;
            end;
        end;
        m = 4;
        table2 = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
    % 64-QAM modulation
    elseif strcmp(mod_type,'64QAM')
        % generates 64QAM symbols
        m=1;
        for k=-7:2:7
            for l=-7:2:7
                table(m) = (k+j*l)/sqrt(42); % power normalization
                m=m+1;
            end;
        end;
        m = 6;
        table2 = [0 1 3 2 6 7 5 4 8 9 11 10 14 15 13 12 24 25 27 26 30 31 29 28 16 17 19 18 22 23 21 20 48 49 51 50 54 55 53 52 56 57 59 58 62 63 61 60 40 41 43 42 46 47 45 44 32 33 35 34 38 39 37 36];

    else
        error('Unimplemented modulation');
    end
end