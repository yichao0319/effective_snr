
function generate_tx_data()

  s = RandStream('mcg16807','Seed',0);    % fix the seed 
  RandStream.setGlobalStream(s);          % so it generates the same random bit stream each time

  num_sc = 48;
  num_ofdmsym = 24;

  output_dir = 'rawofdm.4modulation/';
  
  fid1 = fopen([output_dir 'tx_bits.dat'], 'w');
  fid2 = fopen([output_dir 'tx_syms.dat'], 'w');
  fid3 = fopen([output_dir 'tx_syms_plane.dat'], 'w');


  %% -----------------------------------------
  % BPSK
  modulation = 'BPSK';
  [table, m, table2] = mod_table(modulation);
  k = power(2, m);
  x = randi([1 k], num_sc, num_ofdmsym);
  syms = table(x);
  syms_line = reshape(syms, [], 1);
  syms_real = real(syms_line)';
  syms_imag = imag(syms_line)';
  syms_combine = [syms_real; syms_imag];
  syms_save = reshape(syms_combine, [], 1);

  bits = symbol2bit(modulation, syms_line);
  % size(bits)
  bits_line = reshape(bits', [], 1);
  % size(bits_line)
  fwrite(fid1, bits_line, 'uint8');
  fwrite(fid2, syms_save, 'float');
  fprintf(fid3, '%.4f, %.4f\n', [real(syms_line)'; imag(syms_line)']);


  %% -----------------------------------------
  % QPSK
  modulation = 'QPSK';
  [table, m, table2] = mod_table(modulation);
  k = power(2, m);
  x = randi([1 k], num_sc, num_ofdmsym);
  syms = table(x);
  syms_line = reshape(syms, [], 1);
  syms_real = real(syms_line)';
  syms_imag = imag(syms_line)';
  syms_combine = [syms_real; syms_imag];
  syms_save = reshape(syms_combine, [], 1);

  bits = symbol2bit(modulation, syms_line);
  % size(bits)
  bits_line = reshape(bits', [], 1);
  % size(bits_line)
  fwrite(fid1, bits_line, 'uint8');
  fwrite(fid2, syms_save, 'float');
  fprintf(fid3, '%.4f, %.4f\n', [real(syms_line)'; imag(syms_line)']);

  % tmp = syms_line(1:10)
  % fprintf('%.4f, %.4f\n', [real(tmp)'; imag(tmp)']);


  %% -----------------------------------------
  % 16QAM
  modulation = '16QAM';
  [table, m, table2] = mod_table(modulation);
  k = power(2, m);
  x = randi([1 k], num_sc, num_ofdmsym);
  syms = table(x);
  syms_line = reshape(syms, [], 1);
  syms_real = real(syms_line)';
  syms_imag = imag(syms_line)';
  syms_combine = [syms_real; syms_imag];
  syms_save = reshape(syms_combine, [], 1);

  bits = symbol2bit(modulation, syms_line);
  % size(bits)
  bits_line = reshape(bits', [], 1);
  % size(bits_line)
  fwrite(fid1, bits_line, 'uint8');
  fwrite(fid2, syms_save, 'float');
  fprintf(fid3, '%.4f, %.4f\n', [real(syms_line)'; imag(syms_line)']);


  %% -----------------------------------------
  % 64QAM
  modulation = '64QAM';
  [table, m, table2] = mod_table(modulation);
  k = power(2, m);
  x = randi([1 k], num_sc, num_ofdmsym);
  syms = table(x);
  syms_line = reshape(syms, [], 1);
  syms_real = real(syms_line)';
  syms_imag = imag(syms_line)';
  syms_combine = [syms_real; syms_imag];
  syms_save = reshape(syms_combine, [], 1);

  bits = symbol2bit(modulation, syms_line);
  % size(bits)
  bits_line = reshape(bits', [], 1);
  % size(bits_line)
  fwrite(fid1, bits_line, 'uint8');
  fwrite(fid2, syms_save, 'float');
  fprintf(fid3, '%.4f, %.4f\n', [real(syms_line)'; imag(syms_line)']);



  fclose(fid1);
  fclose(fid2);
  fclose(fid3);



end
