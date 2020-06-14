function [F_ind,B_ind,U_ind,F_rgb,B_rgb,U_rgb,F_s,B_s,U_s,F_mindist,B_mindist] = GetMattingInfo(img,trimap)
%GETMATTINGINFO 此处显示有关此函数的摘要
%   此处显示详细说明
    F_ind = find(trimap == 255);
    B_ind = find(trimap == 0);
    U_ind = find(trimap == 128);
    img_rgb = single(reshape(img,[numel(trimap),3]));
    %% 样本属性计算
    % 颜色
    F_rgb = img_rgb(F_ind,:);
    B_rgb = img_rgb(B_ind,:);
    U_rgb = img_rgb(U_ind,:);
    % 距离
    [F_y,F_x] = ind2sub(size(trimap),F_ind); F_yx = single([F_y,F_x]);
    [B_y,B_x] = ind2sub(size(trimap),B_ind); B_yx = single([B_y,B_x]);
    [U_y,U_x] = ind2sub(size(trimap),U_ind); U_yx = single([U_y,U_x]);
    F_s  = [F_y,F_x];
    B_s  = [B_y,B_x];
    U_s  = [U_y,U_x];
    
    % 最小距离
    F_mindist = bwdist(trimap == 255);F_mindist = F_mindist(U_ind);
    B_mindist = bwdist(trimap == 0);B_mindist = B_mindist(U_ind);
end

