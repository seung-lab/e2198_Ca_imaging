function [nodes, edges, root_node] = TEASAR(im_points, input_resolution, dsmp_resolution, scale, const)
    %im_points - n x 3 array of point coordinates
    %input_resolution - self explanatory, difference in real length between
    %                   point [1 1 1] and [0 0 0]
    %dsmp_resolution - resolution of the coordinates after downsampling
    %scale and const - "scale" and "const" parameter from TEASAR paper,
    %                  larger values mean fewer hairs
    %


    NDIMS = 3;
    quick_ind = @(S, I) (sub2ind(S, I(:,1), I(:,2), I(:,3)));
    
    im_min = double(min(im_points));
    
    for d = 1:NDIMS
        im_points(:,d) = 2+floor((im_points(:,d)-im_min(d))*input_resolution(d)/dsmp_resolution(d));
    end
    
    im_max = double(max(im_points))+1;
    
    bin_im = false(im_max);
    
    bin_im(quick_ind(im_max, im_points)) = true;
%     bin_im(sub2ind(im_max, im_points(:,1), im_points(:,2), im_points(:,3))) = true;
    
    clear im_points;

    node2ind = find(bin_im(:));
    N = length(node2ind);
    ind2node = sparse(node2ind,ones(N,1),(1:N)',prod(im_max),1);
    
    
    
    
    DBF = bwdist(~bin_im);
    
    
    %below is a change I added to avoid problems with flat portions of the
    %neurite. Change DBF_adjustment_constant to 0 or delete to remove.
    
    DBF_adjustment_constant = (sqrt(2)-1)/100; %pretty small number, this is just to break ties
    for d = 1:NDIMS
        not_d = 1:NDIMS;
        not_d(d) = [];
        
        temp_bin_im = permute(bin_im, [not_d d]);
        DBF_inc = zeros(size(temp_bin_im));
        for n = 1:size(bin_im,d)
            DBF_inc(:,:,n) = bwdist(~temp_bin_im(:,:,n));
        end
        
        DBF = DBF + DBF_adjustment_constant*permute(DBF_inc, [1:d-1 NDIMS d:NDIMS-1]);
    end
            
            
    
    
    M = max(DBF(:))*1.01;
    
    p_v = 5000*(1-DBF/M);
    
    
    %create graph, use euclidean distance for edge weights
    %26-connectivity
    
    [nhood(:,1) nhood(:,2) nhood(:,3)] = ind2sub([3 3 3], 1:27);
    nhood = nhood - 2;
    nhood(all(nhood==0,2),:) = [];
    hood_weight = sqrt(sum(nhood.^2,2));
    hood_length = size(nhood,1);
    
    dest_node = zeros(N,hood_length);
    my_node = zeros(N,hood_length);
    edge_weight = zeros(N,hood_length);
    
    for k = 1:N
        my_node(k,:) = k;
        my_ind = node2ind(k);        
        [my_loc(1) my_loc(2) my_loc(3)] = ind2sub(im_max, my_ind);            
        
        for n = 1:hood_length
            dest_loc = my_loc + nhood(n,:);
            if all(dest_loc > 0) && all(dest_loc <= im_max)
                dest_node(k,n) = ind2node(quick_ind(im_max,dest_loc));
                if dest_node(k,n) ~= 0
                    edge_weight(k,n) = hood_weight(n) * p_v(node2ind(dest_node(k,n)));
                    
                    %edge_weight(k,n) = hood_weight(n) + p_v(node2ind(dest_node(k,n)));                    
                    % ^ is as written in the paper. This makes little
                    % sense, imo, since it just makes the actual length of
                    % the path irrelevant
                    
                end
            end
        end
    end
    
    has_edge = dest_node(:) ~= 0;
    G = sparse(my_node(has_edge), dest_node(has_edge), edge_weight(has_edge), N, N);
    
    
    %the commented code below is intended to find the soma and choose it as
    %the root node. For general purpose, we'll instead choose an arbitrary
    %point that's probably a leaf or very near one.
        
    
%     strat = zeros(im_max(1),1);
%     for c = 1:im_max(1)
%         strat(c) = sum(bin_im(c,:));
%     end
% 
%     cumstrat = cumsum(strat);
%     cumstrat = cumstrat/cumstrat(end);
%     
%     if ~exist('root_position','var') || isempty(root_position)    
%         %guess soma loc from stratification
%         
%         strat_med = find(cumstrat>.5,1,'first');
%         if strat_med > im_max(1)/2
%             root_position = 'GCL';
%         else
%             root_position = 'INL';
%         end
%         
%     end
%     
%     if strcmp(root_position, 'GCL');
%         slice_pos = find(cumstrat>=.05,1,'first');
%     else
%         slice_pos = find(cumstrat>=.95,1,'first');        
%     end
%         
%     im_slice = squeeze(bin_im(slice_pos,:,:));
%     [x_coord y_coord] = find(im_slice);
%     x_mean = mean(x_coord(:));
%     y_mean = mean(y_coord(:));
%     
%     dist_from_mean = (x_coord-x_mean).^2 + (y_coord-y_mean).^2;
%     [~, min_loc] = min(dist_from_mean);
%     
%     root_location = [slice_pos, x_coord(min_loc), y_coord(min_loc)];
%     root_node = ind2node(quick_ind(im_max, root_location));
    
    root_ind = find(bin_im(:),1,'first');
    root_node = ind2node(root_ind);
    

    
    PDRF_G = graphshortestpath(G,root_node);
    PDRF = zeros(im_max);
    PDRF(node2ind) = PDRF_G;
    
    
    
    paths = [];
    k = 1;
    
    mask_im = bin_im & ~isinf(PDRF);
    
    
    while any(mask_im(:))
        
        mask_inds = find(mask_im(:));
        
        [~, farthest_in_ind] = max(PDRF(mask_inds));
        farthest_node = ind2node(mask_inds(farthest_in_ind));
        
        [~, paths{k}] = graphshortestpath(G,root_node,farthest_node);
        
        %prob a faster way to do this, but lazy
        for c = 1:length(paths{k})
            my_ind = node2ind(paths{k}(c));
            [my_loc(1) my_loc(2) my_loc(3)] = ind2sub(im_max, my_ind);
            
            r = DBF(my_ind)*scale + const;
            
            r_begin = ceil(my_loc - r);
            r_begin(r_begin<1) = 1;
            
            r_end = floor(my_loc + r);
            r_end = min([r_end; im_max]);
            
            mask_im(r_begin(1):r_end(1), r_begin(2):r_end(2), r_begin(3):r_end(3)) = false;
        end
        
        k = k + 1;
        
    end
            
        
        
    %consolidate paths
    
    tree_nodes = [];
    for k = 1:length(paths);
        tree_nodes = [tree_nodes; paths{k}'];
    end
    tree_nodes = unique(tree_nodes);
    tree_size = length(tree_nodes);
    
    inv_tree_nodes = sparse(tree_nodes, ones(tree_size,1), (1:tree_size)', N,1);
    root_node = inv_tree_nodes(paths{1}(1));
    
    edges = []; %probably losing a lot of speed here by not preallocating, but lazy
    for k = 1:length(paths);
        edges = [edges; full([inv_tree_nodes(paths{k}(1:end-1)') inv_tree_nodes(paths{k}(2:end)')])];
    end
    edges = unique(edges, 'rows');
    
    
    [nodes(:,1) nodes(:,2) nodes(:,3)] = ind2sub(im_max, node2ind(tree_nodes));
    
    for d = 1:NDIMS
        nodes(:,d) = (nodes(:,d)-2)*dsmp_resolution(d)/input_resolution(d) + im_min(d);
    end
%     im_points(:,d) = 2+floor((im_points(:,d)-im_min(d))*input_resolution(d)/dsmp_resolution(d));
    
    
    
    
    
end
        
    
    