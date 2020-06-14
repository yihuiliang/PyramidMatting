function [ alpha, fitness ] = FB2alpha( FBpairs, img, trimap, is_sort )
%FB2ALPHA 将FB输出成alpha
%   Detailed explanation goes here
    F_ind = find(trimap == 255);
    B_ind = find(trimap == 0);
    U_ind = find(trimap == 128);
    img_rgb = single(reshape(img,[numel(trimap),3]));
    F_rgb = img_rgb(F_ind,:); 
    B_rgb = img_rgb(B_ind,:); 
    U_rgb = img_rgb(U_ind,:); 
    [F_y,F_x] = ind2sub(size(trimap),F_ind); F_yx = [F_y,F_x];
    [B_y,B_x] = ind2sub(size(trimap),B_ind); B_yx = [B_y,B_x];
    [U_y,U_x] = ind2sub(size(trimap),U_ind); U_yx = [U_y,U_x];
    %% 前景背景点按照HSV排序
    if is_sort
        if is_sort ==1
            img_hsv = single(reshape(rgb2hsv(img),[numel(trimap),3]));
        else
            img_hsv = single(reshape(rgb2gray(img),[numel(trimap),1]));
        end
        F_hsv = img_hsv(F_ind,:);
        B_hsv = img_hsv(B_ind,:);
%         U_hsv = img_hsv(U_ind,:);
        [~,F_hsv_ind] = sortrows(F_hsv); 
        [~,B_hsv_ind] = sortrows(B_hsv); 
        F_rgb = F_rgb(F_hsv_ind,:);      B_rgb = B_rgb(B_hsv_ind,:);     
        F_yx = F_yx(F_hsv_ind,:);        B_yx = B_yx(B_hsv_ind,:);
    end  
    
    
    Fx_rgb = F_rgb(FBpairs(:,1),:);
    Bx_rgb = B_rgb(FBpairs(:,2),:);
    Fx_Bx_rgb = Fx_rgb-Bx_rgb;
    alpha_U = sum((U_rgb-Bx_rgb).*(Fx_Bx_rgb),2)./(sum(Fx_Bx_rgb.*Fx_Bx_rgb,2)+1);
    
    alpha = trimap;
    alpha(U_ind) = 255*alpha_U;
    
    %% fitness    
    alpha_U(alpha_U>1) = 1;
    alpha_U(alpha_U<0) = 0;
    alphaX3 = [alpha_U,alpha_U,alpha_U];
    Obj1 = sum((U_rgb - (alphaX3.*Fx_rgb+(1-alphaX3).*Bx_rgb)).^2,2).^0.5;
    
    %% 最小距离计算
    F_mindist = bwdist(trimap == 255);F_mindist = F_mindist(U_ind);
    B_mindist = bwdist(trimap == 0);B_mindist = B_mindist(U_ind);
    Fx_yx = F_yx(FBpairs(:,1),:);
    Bx_yx = B_yx(FBpairs(:,2),:);
    Obj2 = sum((Fx_yx-U_yx).^2,2).^0.5./F_mindist;
    Obj3 = sum((Bx_yx-U_yx).^2,2).^0.5./B_mindist;
    fitness = mean(sum([Obj1,Obj2,Obj3]'));
end

