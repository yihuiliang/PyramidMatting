function [FB_pairs] = MatchFBPairs(pyramid,n)
%MATCHFBPAIRS 此处显示有关此函数的摘要
%   此处显示详细说明
        img = pyramid{n,1};
        trimap = pyramid{n,2};
        F_ind = find(trimap == 255);
        B_ind = find(trimap == 0);
        U_ind = find(trimap == 128);
        FB_pairs = zeros(length(U_ind),2);
        %% 不同层之间像素对匹配
        trimap_s = pyramid{n+1,2};
        ind_s = reshape(1:numel(pyramid{n+1,2}),size(trimap_s));
        map_ind_l2s = repelem(ind_s,2,2); 
        map_ind_l2s = map_ind_l2s(1:size(trimap,1),1:size(trimap,2));
        s = regionprops(map_ind_l2s,'PixelIdxList');
        map_ind_s2l = struct2cell(s);
        %建立全局坐标到局部坐标的映射
        g2l_Find = zeros(size(trimap)); g2l_Find(F_ind) = 1:length(F_ind);
        g2l_Bind = zeros(size(trimap)); g2l_Bind(B_ind) = 1:length(B_ind);
        g2l_Uind = zeros(size(trimap)); g2l_Uind(U_ind) = 1:length(U_ind);
        
        g2l_Uind_s = zeros(size(trimap_s)); g2l_Uind_s(trimap_s==128) = 1:nnz(trimap_s==128);
        
        FB_pairs_s = pyramid{n+1,4};
        for i = 1:length(U_ind)
            % 未知像素之间的映射
            U_ind_i = U_ind(i);
            U_ind_s_i = map_ind_l2s(U_ind_i);
            share_ind = map_ind_s2l{U_ind_s_i};
            share_ind = g2l_Uind(share_ind);
            share_ind(share_ind==0) = [];
            
            % 已知像素的映射
            local_U_ind_s = g2l_Uind_s(U_ind_s_i);
            if local_U_ind_s ==0
                continue;
            end
            FB_s_i = FB_pairs_s(local_U_ind_s,:);

            F_i_cand_ind = g2l_Find(map_ind_s2l{FB_s_i(1)});
            F_i_cand_ind(F_i_cand_ind == 0) = [];
            F_i_cand_ind = F_i_cand_ind(randperm(length(F_i_cand_ind)));
            if numel(F_i_cand_ind)<numel(share_ind)
                F_i_cand_ind = repmat(F_i_cand_ind,numel(share_ind),1);
            end
            B_i_cand_ind = g2l_Bind(map_ind_s2l{FB_s_i(2)});
            B_i_cand_ind(B_i_cand_ind == 0) = [];
            B_i_cand_ind = B_i_cand_ind(randperm(length(B_i_cand_ind)));            
            if numel(B_i_cand_ind)<numel(share_ind)
                B_i_cand_ind = repmat(B_i_cand_ind,numel(share_ind),1);
            end          
            
            FB_pairs(share_ind,1) = F_i_cand_ind(1:length(share_ind));
            FB_pairs(share_ind,2) = B_i_cand_ind(1:length(share_ind));
        end
        bw = false(size(trimap));
        bw(U_ind(FB_pairs(:,1)~=0)) = true;
        [~,U_mindist_idx] = bwdist(bw);
        for i = 1:length(U_ind)
            if FB_pairs(i,1) ==0
                nearest_U_ind = g2l_Uind(U_mindist_idx(U_ind(i)));
                FB_pairs(i,:) = FB_pairs(nearest_U_ind,:);
            end
        end
end

