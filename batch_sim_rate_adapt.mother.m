
function batch_analyze_csi()

    input_sym_files = {'XXXXX'};
    
    
    % output_dir = './OUTPUT_sim/';
    % fid = fopen([output_dir 'batch_sim.txt'], 'w');

    for file_i = 1:length(input_sym_files)
        input_sym_file = char(input_sym_files(file_i));
        sim_rate_adapt(input_sym_file, 'THRESHOLD');
    end

    % fclose(fid);

end