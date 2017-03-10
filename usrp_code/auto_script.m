cd ~/Desktop/matlab/;
load ~/data.mat;
n_fft = 64;
n_occupied = 52;
[time_symbol, freq_symbol] = create_symbol(data, n_fft, n_occupied);
ofdm_data = read_complex_binary('/home/zijian/Desktop/ofdm_recv.dat');
%ofdm_data = read_complex_binary('/home/zijian/Desktop/doppler_recv6.dat');
%ofdm_data = read_complex_binary('/home/zijian/Desktop/uhd/host/build/examples/ofdm_recv_0.dat');
%ofdm_data = read_complex_binary('/home/zijian/Desktop/remove_system_offset.dat');
[row, line] = size(ofdm_data);
start_index = 1e5;
all_channels = ajacent_channels(ofdm_data, freq_symbol, start_index, row, n_fft, n_occupied);
%origin_channels = ajacent_channels(origin_data, freq_symbol, start_index, row, n_fft, n_occupied);

%all_channels = all_channels ./ origin_channels;
% get the offsets
temp = angle(all_channels(convert_bin_index_normal_to_fft([-1], 64), 1:end));

%%%%%%%%%%%
%sample_end = 1e5;
%sample_angle = temp(1:sample_end);
%estimate_result = estimate_jump(sample_angle, 30, 2);
%tmp1 = abs(estimate_result);
%pre_process = tmp1 - 10 * std(tmp1);
%pre_process(pre_process<0,:) = 0;
%target_index = find(pre_process);
%index = select_neighbor(target_index, 30);
%start_jmp_index = index(1);
%period = median(diff(index));
%values = get_period_median(temp, start_jmp_index, period);
%jumps = diff(values);
%jump = median(jumps);
%%%%%%%%%%%%%%%
mean(unwrap(temp))
figure;plot(unwrap(temp));

jump = 0.1;
period = 15625;
%start_jmp_index = 9931;
start_jmp_index = 10451;
result = remove_shift(temp, start_jmp_index, period, jump);
figure;plot(temp);
hold on;plot(unwrap(result));
offsets = repelem(0, 7);
for i=-26:-1
    %hold on;plot(unwrap(angle(all_channels(convert_bin_index_normal_to_fft([i], 64), 2000:8000))));
    x = 1:8000;
    y = unwrap(angle(all_channels(convert_bin_index_normal_to_fft([i], 64), x)));
    p = polyfit(x, y, 1);
    offsets = [offsets p(1)];
end
offsets = [offsets 0];
for i = 1:26
    x = 1:8000;
    y = unwrap(angle(all_channels(convert_bin_index_normal_to_fft([i], 64), x)));
    p = polyfit(x, y, 1);
    offsets = [offsets p(1)];
end
for i = 27:31
    offsets = [offsets 0];
end
%offsets = offsets(2:end) ./ (2*pi);
offsets = offsets(2:end) ./ 64;
new_data = ofdm_data;
for i = start_index:n_fft:row-n_fft
    new_data(i:i+n_fft-1) = ofdm_data(i:i+n_fft-1) .* (exp((offsets .* (i:i+n_fft-1)) * 1i) )';
end
new_channels = ajacent_channels(new_data, freq_symbol, start_index, row, n_fft, n_occupied);

for i = -26:-26
    hold on; plot(unwrap(angle(new_channels(convert_bin_index_normal_to_fft([i], 64), x))));
    hold on; plot(unwrap(angle(all_channels(convert_bin_index_normal_to_fft([i], 64), x))));
end