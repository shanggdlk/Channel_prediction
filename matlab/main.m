clc,clear
% csi_trace = read_bf_file('sample_data/testmove1.dat');
ath_csi_trace = read_log_file('sample_data/move2.dat');
% num_csi = length(csi_trace);
num_ath_csi = length(ath_csi_trace);
csi_trace = cell(0);
id = 0;
for i = 1:num_ath_csi
    if ~ath_csi_trace{i}.csi_len
        continue
    end
    csi_trace = [csi_trace; cell(1)];
    id = id + 1;
    csi_trace{id}.csi = permute(ath_csi_trace{i}.csi,[2,1,3]);
    csi_trace{id}.rssi1 = double(ath_csi_trace{i}.rssi1);
    csi_trace{id}.rssi2 = double(ath_csi_trace{i}.rssi2);
    csi_trace{id}.rssi3 = double(ath_csi_trace{i}.rssi3);
    csi_trace{id}.nr = double(ath_csi_trace{i}.nr);
    csi_trace{id}.nc = double(ath_csi_trace{i}.nc);
end
num_csi = length(csi_trace);
esnr = zeros(1,num_csi);
rssi = zeros(1,num_csi);
unscaled_csi = zeros(1, num_csi);
for i = 1:num_csi
    csi_entry = csi_trace{i};
    rssi(i) = csi_entry.rssi1;
    tmp_unscaled_csi = csi_entry.csi;
    tmp_unscaled_esnr = db(get_eff_SNRs(tmp_unscaled_csi), 'pow');
    unscaled_csi(i) = tmp_unscaled_esnr(1,4);
    csi = get_scaled_csi(csi_entry);
    tmp_esnr = db(get_eff_SNRs(csi), 'pow');
    esnr(i) = tmp_esnr(1, 4);
    if isinf(esnr(i))
        esnr(i) = 10;
    end
    
end
figure
plot(esnr);
%ylim([0,50]);
title('esnr')

figure
plot(rssi);
title('rssi');
%ylim([0,50]);

figure
plot(unscaled_csi);
title('unscaled csi');
%ylim([0,50]);
 
 
 