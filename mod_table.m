%% mod_table: function description
function [table, m] = mod_table(mod_type)
    % BPSK
    if strcmp(mod_type,'BPSK')
        m = 1;
        table=exp(j*[0 -pi]);  % generates BPSK symbols
    % QPSK modulation
    elseif strcmp(mod_type,'QPSK')
        m = 2;
        table=exp(j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
    % 16-QAM modulation
    elseif strcmp(mod_type,'16QAM')
        % generates 16QAM symbols
        m=4;
        for k=-3:2:3
            for l=-3:2:3
                table(m) = (k+j*l)/sqrt(10); % power normalization
                m=m+1;
            end;
        end;
    % 64-QAM modulation
    elseif strcmp(mod_type,'64QAM')
        % generates 64QAM symbols
        m=6;
        for k=-7:2:7
            for l=-7:2:7
                table(m) = (k+j*l)/sqrt(42); % power normalization
                m=m+1;
            end;
        end;
    else
        error('Unimplemented modulation');
    end
end