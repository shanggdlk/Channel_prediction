my_dir = '/Users/zhenyus/github/Channel_prediction/csi5g1/';
prefix_T1 = 'csi_';
D = dir([my_dir,'*.log']);
fnum =size(D,1);                       
%Read files
c_quality = cell(1,fnum);
ESNR = [];
for i=1:5
    filename = [my_dir,D(i).name];
    csi_entries = read_log_file(filename);
    csi_entry = csi_entries{1};
    %esnr = esnr_extract(csi_entry);
    dd = csi_entry.csi;
    c_quality{i} = dd;
    
esnr = esnr_extract(csi_entry);
g = [esnr(1,4);esnr(2,4)];
ESNR = [ESNR,g];
end


ESNR = [];
for i = 1:6
%csi = c_quality{i}(:,1,:);
%csi = reshape(csi,3,56);
%plot(db(abs(csi.')))
esnr = esnr_extract(csi_entry);
ESNR = [ESNR,esnr];
end

%csi_entries = read_log_file(csi_file_name);
%csi_entry = csi_entries{1};
%g = csi_entry.csi;
%esnr = esnr_extract(csi_entry);