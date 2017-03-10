function channels = ajacent_channels(ofdm_data, freq_symbol, start_index, end_index, n_fft, n_occupied)
    channels= zeros(64,1);
    for i=start_index:n_fft:(end_index - n_fft)
        c = get_channel(ofdm_data, freq_symbol, i, n_fft, n_occupied);
        channels = [channels c];
    end
    channels = channels(1:end, 2:end);
end