%% symbol2bit: function description
function [bit] = symbol2bit(mod_type, sym)
    [row col] = size(sym);
    table = mod_table(mod_type);

    % BPSK
    if strcmp(mod_type,'BPSK')
        bit = zeros(row, col * 2);
    % QPSK modulation
    elseif strcmp(mod_type,'QPSK')
        m = 2;
        bit = zeros(row, col * m);
        for row_i = 1:row
            for col_i = 1:col
                real_part = real(sym(row_i, col_i) );
                imag_part = imag(sym(row_i, col_i) );
                if(real_part > 0)
                    bit(row_i, 2*(col_i-1)+1) = 1;
                else
                    bit(row_i, 2*(col_i-1)+1) = 0;
                end

                if(imag_part > 0)
                    bit(row_i, 2*col_i) = 1;
                else
                    bit(row_i, 2*col_i) = 0;
                end
            end
        end

    % 16-QAM modulation
    elseif strcmp(mod_type,'16QAM')
        bit = zeros(row, col * 2);
    % 64-QAM modulation
    elseif strcmp(mod_type,'64QAM')
        bit = zeros(row, col * 2);
    else
        error('Unimplemented modulation');
    end
end