function [ esnr ] = esnr_extract(csi_entry)
%   compute esnr from normalizd csi and rssi value
%   created by Zhenyu Song, sunnyszy@gmail.com

    csi = get_scaled_csi(csi_entry);
    csi_5300format = permute(csi, [2 1 3]);
%     size(csi)
%     esnr1 = db(get_eff_SNRs(csi), 'pow');
    esnr = db(get_eff_SNRs(csi_5300format), 'pow');
    
end

