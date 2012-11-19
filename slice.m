%% ----------------------------------
% slice: 
%    EVM denote the euclidean distance to the closest constellation point in percentage
function [out_grid, evm] = slice(mod_type, rx_grid)
    table = mod_table(mod_type);
    p0 = mean(abs(table) );   %% should be 1 for qpsk though...
    evm = zeros(size(rx_grid));
    out_grid = zeros(size(rx_grid));

    [row, col] = size(rx_grid);
    for row_i = 1:row
        for col_i = 1:col
            dist = zeros(length(table), 1);
            for dist_i = 1:length(table)
                dist(dist_i) = abs(rx_grid(row_i, col_i) - table(dist_i));
            end
            [m, index] = min(dist);
            evm(row_i, col_i) = m / p0;
            out_grid(row_i, col_i) = table(index);
        end
    end
end