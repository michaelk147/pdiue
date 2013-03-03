clear all;
close all;
dim = 3;

searchspace = [0, 0, 0; 1, 1, 1];



N = size(ps,1);


figure(1);
hold off;
plot3(ps(:,2),ps(:,3),ps(:,1),'xb');
hold on;

