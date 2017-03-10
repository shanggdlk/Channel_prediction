function channel = get_channel(ofdm_data, freq_symbol, offset, n_fft, n_occupied)
    channel = fft(ofdm_data(offset:offset+n_fft-1))./freq_symbol; 
    channel(convert_bin_index_normal_to_fft([-n_fft/2:(-n_occupied/2-1) 0 (n_occupied/2+1):n_fft/2-1], n_fft)) = NaN;
end