function [ d ] = mahalanobisDist( x1,x2,Hinv)
%MAHALANOBISDIST Summary of this function goes here
%   Detailed explanation goes here
    dif = (x1 -x2);
    d =  dif * Hinv * dif';
    
end

