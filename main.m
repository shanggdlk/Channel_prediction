clc,clear
globals_init();
% --------------------------------------

my_dir = 'csi5g1/';
prefix_T1 = 'csi_';
D = dir([my_dir,'*.log']);
fnum =size(D,1);                       
%Read files
% different frequency
for i=1:1
    filename = [my_dir,D(i).name];
    csi_entries = read_log_file(filename);
    n_entry = length(csi_entries);
    csis = cell(1, n_entry);
    for j=1:length(csi_entries)
        csi = get_scaled_csi(csi_entries{j});
        csi_5300format = permute(csi, [2 1 3]);
        csis{j} = csi_5300format;
    end
end
    
csi_amp = zeros(1, n_entry);
csi_angle = zeros(1, n_entry);
for i = 1:n_entry
    csi_amp(i) = abs(csis{i}(1,1,1));
    csi_angle(i) = angle(csis{i}(1,1,1));
end
plot(1:n_entry, csi_amp);
title('csi amplitude');
figure
plot(1:n_entry, csi_angle);
title('csi phase');
