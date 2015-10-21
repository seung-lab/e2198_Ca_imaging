function [roi_struct] = f(full_filename,typeclrs)

%
% load in the ROI data from the file, w/ error checking
%

% open the file
% full_filename=strcat(pathname,filename);
fid=fopen(full_filename,'r','ieee-be');
if (fid == -1)
  errordlg(sprintf('Unable to open file %s',filename),...
           'File Error');
  return;
end

% read the number of rois
[n_rois,count]=fread(fid,1,'uint32');
%n_rois
if (count ~= 1)
  errordlg(sprintf('Error loading ROIs from file %s',filename),...
           'File Error');
  fclose(fid);
  return;
end

% dimension cell arrays to hold the ROI labels and vertex lists
labels=cell(n_rois,1);
borders=cell(n_rois,1);
types=zeros(n_rois,1);

% for each ROI, read the label and the vertex list
for j=1:n_rois
  % the label
  [n_chars,count]=fread(fid,1,'uint32');
  if (count ~= 1)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  [temp,count]=fread(fid,[1 n_chars],'uchar');
  if (count ~= n_chars)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  labels{j}=char(temp);
  [temp,count]=fread(fid,1,'uint8');
  types(j) = temp;
  % the vertex list
  [n_vertices,count]=fread(fid,1,'uint32');
  if (count ~= 1)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  this_border=zeros(2,n_vertices);
  [this_border,count]=fread(fid,[2 n_vertices],'float32');
  %this_border
  if (count ~= 2*n_vertices)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  borders{j}=this_border;
end

% close the file
fclose(fid);

roi_ids=(1:n_rois)';
label_h=zeros(n_rois,1);
border_h=zeros(n_rois,1);

for j=1:n_rois
  this_border=borders{j};
  com=border_com(this_border);
  label_h(j)=...
    text('Parent',gca,...
         'Position',[com(1) com(2) 1],...
         'String',labels{j},...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Color',typeclrs(types(j),:),...
         'Tag','label_h',...
         'Clipping','on');
  border_h(j)=...
    line('Parent',gca,...
         'Color',typeclrs(types(j),:),...
         'Tag','border_h',...
         'XData',this_border(1,:),...
         'YData',this_border(2,:),...
         'ZData',repmat(1,[1 size(this_border,2)]));
end

roi_struct = struct('roi_ids',roi_ids,....
                    'border_h',border_h,...
                    'label_h',label_h,...
                    'type',types);
   