function [ output ] = norm2( x )
%NORM2 Summary of this function goes here
%   Detailed explanation goes here
    output = sum(x.^2,2).^(1/2);

end

