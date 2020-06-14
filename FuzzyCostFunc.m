function [ fitness,est_alpha,cost_c,cost_sF,cost_sB,cost_cd] = ...
FuzzyCostFunc( x,F_rgb,B_rgb,U_rgb,F_s,B_s,U_s,F_mindist,B_mindist)
%COSTFUNC Fitness function of alpha matting
%   Detailed explanation goes here
    if size(x,2)~=2
        error('x must be nx2 matrix');
    end
    if size(U_rgb,1) == 1
        U_rgb = repmat(U_rgb,size(x,1),1);
        U_s = repmat(U_s,size(x,1),1);
    end
    x = round(x);
    x_F = x(:,1); x_B = x(:,2);
    
    Fx_rgb = F_rgb(x_F,:);
    Bx_rgb = B_rgb(x_B,:);
    Fx_s =  F_s(x_F,:);
    Bx_s =  B_s(x_B,:);
    % Alpah
    Fx_Bx_rgb = Fx_rgb - Bx_rgb;
    est_alpha = sum((U_rgb - Bx_rgb).*Fx_Bx_rgb,2)./(sum(Fx_Bx_rgb.*Fx_Bx_rgb,2)+1);
    est_alpha(est_alpha>1) = 1;
    est_alpha(est_alpha<0) = 0;
    % Chromatic distortion
    cost_c  = norm2(U_rgb-(est_alpha.*Fx_rgb+(1-est_alpha).*Bx_rgb));
    % Spatial cost
    cost_sF = norm2(Fx_s-U_s)./F_mindist;
    cost_sB = norm2(Bx_s-U_s)./B_mindist;
    cost_cd = norm2(Fx_Bx_rgb);
    
    mem_sF = exp(-0.17*(cost_sF-1));
    mem_sB = exp(-0.17*(cost_sB-1));
    mem_cd = exp(-0.11*(cost_c));
    fitness = -EinsteinIntersect(mem_cd,0.5*(mem_sF+mem_sB));
end

