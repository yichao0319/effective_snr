%% symbol2bit: function description
function [bit] = symbol2bit(mod_type, sym)
    [row col] = size(sym);
    [table, m, table2] = mod_table(mod_type);
    

    bit = zeros(row, col * m);
    for row_i = 1:row
        for col_i = 1:col
            ind = find(table == sym(row_i, col_i));
            str = dec2bin(table2(ind), m);
            for bi = 1:m
                bit(row_i, m*(col_i-1)+bi) = str2num(str(bi));
            end
        end
    end


    % % BPSK
    % if strcmp(mod_type,'BPSK')
    %     bit = zeros(row, col * m);
    %     for row_i = 1:row
    %         for col_i = 1:col
    %             ind = find(table == sym(row_i, col_i));
    %             bit(row_i, m*(col_i-1)+1:m*col_i) = dec2bin(table2(ind));
    %         end
    %     end

    % % QPSK modulation
    % elseif strcmp(mod_type,'QPSK')
    %     m = 2;
    %     bit = zeros(row, col * m);
    %     for row_i = 1:row
    %         for col_i = 1:col
    %             real_part = real(sym(row_i, col_i) );
    %             imag_part = imag(sym(row_i, col_i) );
    %             if(real_part > 0)
    %                 bit(row_i, m*(col_i-1)+1) = 1;
    %             else
    %                 bit(row_i, m*(col_i-1)+1) = 0;
    %             end

    %             if(imag_part > 0)
    %                 bit(row_i, m*col_i) = 1;
    %             else
    %                 bit(row_i, m*col_i) = 0;
    %             end
    %         end
    %     end

    % % 16-QAM modulation
    % elseif strcmp(mod_type,'16QAM')
    %     bit = zeros(row, col * 2);
    % % 64-QAM modulation
    % elseif strcmp(mod_type,'64QAM')
    %     bit = zeros(row, col * 2);
    % else
    %     error('Unimplemented modulation');
    % end
end