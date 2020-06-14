clear
max_fitness_evaluation_per_unknow_pixel = 5e3;


img_g = imread('GT11.png');
trimap_g = imread('trimap.png');
MFE_p = max_fitness_evaluation_per_unknow_pixel*nnz(trimap_g==128);
FE_counter = 0;
FE_per_pixel = 5e3;
%% constract image pyramid
pyramid{1,1} = img_g;
pyramid{1,2} = trimap_g;
for n = 2:1e10
    pyramid{n,1} = impyramid(pyramid{n-1,1}, 'reduce');
    pyramid{n,2} = imresize(pyramid{n-1,2},0.5,'nearest');
    if max(size(pyramid{n,2}))<100
        break;
    end
end
disp(size(pyramid{n,2}));
%% 逐一优化金字塔中的图像
%第一层
img = pyramid{n,1};
trimap = pyramid{n,2};
[F_ind,B_ind,U_ind,F_rgb,B_rgb,U_rgb,F_s,B_s,U_s,F_mindist,B_mindist] = ...
    GetMattingInfo(img,trimap);
map_g2l_ind = zeros(size(trimap));
map_g2l_ind(F_ind) = 1:length(F_ind);
map_g2l_ind(B_ind) = 1:length(B_ind);
bw = trimap == 128;
bw = imdilate(bw,strel('disk',1));
F_samples = map_g2l_ind(bw&(trimap==255));
B_samples = map_g2l_ind(bw&(trimap==0));
[X,Y] = meshgrid(F_samples,B_samples);
FB_pairs = zeros(length(U_ind),2);
alpha = zeros(length(U_ind),1);
for k = 1:length(U_ind)
    U_ind_k = k;
    U_s_k = U_s(U_ind_k,:);
    U_rgb_k = U_rgb(U_ind_k,:);
    F_mindist_k = F_mindist(U_ind_k,:);
    B_mindist_k = B_mindist(U_ind_k,:);
    
    %% options for ga
    FitnessFcn = @(x) FuzzyCostFunc(x,F_rgb,B_rgb,U_rgb_k,F_s,B_s,U_s_k,F_mindist_k,B_mindist_k);
    numberOfVariables = 2; % Number of decision variables
    lb = [ones(1,numberOfVariables)]; % Lower bound
    ub = [repmat(size(F_rgb,1),1,numberOfVariables/2) repmat(size(B_rgb,1),1,numberOfVariables/2)]; % Upper bound
    
    %% run GA
    [ x,fval] = MyCSO( FitnessFcn,numberOfVariables,lb,ub,FE_per_pixel);
    
    FB_pairs(k,:) = round(x);
    [~,est_alpha_k] = FuzzyCostFunc(x,F_rgb,B_rgb,U_rgb_k,F_s,B_s,U_s_k,F_mindist_k,B_mindist_k);
    alpha(k) = est_alpha_k;
    fitness(k) = fval;
    
end
FE_counter = FE_counter + FE_per_pixel*length(U_ind);
tmp = im2double(trimap);    tmp(U_ind) = alpha;    alpha= tmp;
FB_pairs_info = [F_rgb(FB_pairs(:,1),:),F_s(FB_pairs(:,1),:),...
    B_rgb(FB_pairs(:,2),:),B_s(FB_pairs(:,2),:)];
FB_pairs(:,1) = F_ind(FB_pairs(:,1));       %将FB_pairs坐标转为全局坐标
FB_pairs(:,2) = B_ind(FB_pairs(:,2));
pyramid{n,3} = FB_pairs_info;
pyramid{n,4} = FB_pairs;
pyramid{n,5} = alpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%其余层使用启发式优化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = length(pyramid)-1:-1:1
    img = pyramid{n,1};
    trimap = pyramid{n,2};
    [F_ind,B_ind,U_ind,F_rgb,B_rgb,U_rgb,F_s,B_s,U_s,F_mindist,B_mindist] = ...
        GetMattingInfo(img,trimap);
    FB_pairs = MatchFBPairs(pyramid,n);
    if FE_counter>MFE_p
        break;
    end
    [ fitness,est_alpha,cost_c] = FuzzyCostFunc(FB_pairs,F_rgb,B_rgb,U_rgb,F_s,B_s,U_s,F_mindist,B_mindist);
    FE_counter = FE_counter+length(U_ind);
    %% 根据适应值针对性优化像素对
    [~,cost_ind] = sort(fitness,'descend');
    
    opt_num = ceil((MFE_p-FE_counter)/FE_per_pixel);
    opt_num = min(ceil(length(cost_ind)/2),opt_num);
    roi_U_ind = cost_ind(1:opt_num);
    FB_pairs_roi = FB_pairs(roi_U_ind,:);
    est_alpha_roi = est_alpha(roi_U_ind,:);
    fitness_roi = fitness(roi_U_ind,:);
    
    for k = 1:length(roi_U_ind)
        
        U_ind_k = roi_U_ind(k);
        U_s_k = U_s(U_ind_k,:);
        U_rgb_k = U_rgb(U_ind_k,:);
        F_mindist_k = F_mindist(U_ind_k,:);
        B_mindist_k = B_mindist(U_ind_k,:);
        
        %% options for EA
        FitnessFcn = @(x) FuzzyCostFunc(x,F_rgb,B_rgb,U_rgb_k,F_s,B_s,U_s_k,F_mindist_k,B_mindist_k);
        numberOfVariables = 2; % Number of decision variables
        lb = [ones(1,numberOfVariables)]; % Lower bound
        ub = [repmat(size(F_rgb,1),1,numberOfVariables/2) repmat(size(B_rgb,1),1,numberOfVariables/2)]; % Upper bound
        
        %% run EA
        [ x,fval] = MyCSO_with_initX( FitnessFcn,numberOfVariables,lb,ub,FE_per_pixel,FB_pairs_roi(k,:));
        if fval<fitness_roi(k)
            FB_pairs_roi(k,:) = round(x);
            [~,est_alpha_k] = FuzzyCostFunc(x,F_rgb,B_rgb,U_rgb_k,F_s,B_s,U_s_k,F_mindist_k,B_mindist_k);
            est_alpha_roi(k) = est_alpha_k;
            fitness_roi(k) = fval;
        end
    end
    FB_pairs(roi_U_ind,:) = FB_pairs_roi;
    est_alpha(roi_U_ind,:) = est_alpha_roi;
    fitness(roi_U_ind,:) = fitness_roi;
    
    
    %% 记录该层的信息
    FB_pairs_info = [F_rgb(FB_pairs(:,1),:),F_s(FB_pairs(:,1),:),...
        B_rgb(FB_pairs(:,2),:),B_s(FB_pairs(:,2),:)];
    FB_pairs(:,1) = F_ind(FB_pairs(:,1)); 
    FB_pairs(:,2) = B_ind(FB_pairs(:,2));
    pyramid{n,3} = FB_pairs_info;
    pyramid{n,4} = FB_pairs;
    tmp = im2double(trimap);    tmp(U_ind) = est_alpha;    alpha= tmp;
    pyramid{n,5} = alpha;
    
    FE_counter = FE_counter+opt_num*FE_per_pixel;
    if FE_counter>MFE_p
        break;
    end
end
%检查是否有未优化的层
for i = size(pyramid,1):-1:1
    if ~isempty(pyramid{i,3})
        continue;
    end
    [F_ind,B_ind] = GetMattingInfo(pyramid{i,1},pyramid{i,2});
    FB_pairs = MatchFBPairs(pyramid,i);
    FB_pairs(:,1) = F_ind(FB_pairs(:,1));       
    FB_pairs(:,2) = B_ind(FB_pairs(:,2));
    pyramid{i,4} = FB_pairs;
end
FB_pairs = pyramid{1,4};
map_g2l_ind = zeros(size(trimap_g));
map_g2l_ind(trimap_g ==255) = 1:nnz(trimap_g==255);
map_g2l_ind(trimap_g ==0) = 1:nnz(trimap_g==0);
FB_pairs(:,1) = map_g2l_ind(FB_pairs(:,1));
FB_pairs(:,2) = map_g2l_ind(FB_pairs(:,2));

alpha = FB2alpha(FB_pairs,img_g,trimap_g,0);
imwrite(alpha,'alpha.png');
imshow(alpha);
