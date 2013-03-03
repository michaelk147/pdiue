clear all;
close all;
dim = 3;


pddebugdata
pddebugdata_ssp


N = size(ps,1);
markerscale = 3.0;
gfxlinewd  = 4;

figure(1);
plot3(ps(:,2),ps(:,3),ps(:,1),'+b','MarkerSize', 10, 'linewidth', gfxlinewd);

hold on;
dts = 2;
for i = 1:size(sigmas,1)
    plot(ps(i,3),ps(i,dts),'ob','MarkerSize',sigmas(i,dts),'MarkerEdgeColor',[0.7 0.7 1], 'linewidth', gfxlinewd)
end



figure(2)
%for i = 1:size(sigmas,1)
%    plot(ps(i,3),ps(i,1),'ob','MarkerSize',sigmas(i,2),'MarkerEdgeColor',[1 0.7 0.7])
%end

for i = 1:size(sigmas,1)
    plot(ps(i,3),ps(i,dts),'ob','MarkerSize',sigmas(i,3) * max(ps(:,dts)) / max(ps(:,3)),'MarkerEdgeColor',[0.7 1 0.7], 'linewidth', gfxlinewd)
end

for i = 1:size(sigmas,1)
    plot(mse(i,3),mse(i,dts),'+b','MarkerSize',weights(1,i)*15,'MarkerEdgeColor',[1 0.2 0.2])
end

for i = 1:size(fusedMSE,1)
    plot(fusedMSE(i,3),fusedMSE(i,dts),'ob','MarkerSize',fusedMSE(i,4)*30,'MarkerEdgeColor',[1 0 0])
end

hold off;


figure(7); scatter(ps(:,3),ps(:,1));

hold on;

dts = 1;
for i = 1:size(sigmas,1)
    plot(ps(i,3),ps(i,dts),'ob','MarkerSize',sigmas(i,dts),'MarkerEdgeColor',[0.7 0.7 1])
end

%for i = 1:size(sigmas,1)
%    plot(ps(i,3),ps(i,1),'ob','MarkerSize',sigmas(i,2),'MarkerEdgeColor',[1 0.7 0.7])
%end

for i = 1:size(sigmas,1)
    plot(ps(i,3),ps(i,1),'ob','MarkerSize',sigmas(i,3) * max(ps(:,1)) / max(ps(:,3)),'MarkerEdgeColor',[0.7 1 0.7])
end

for i = 1:size(sigmas,1)
    plot(mse(i,3),mse(i,1),'xb','MarkerSize',weights(1,i)*10,'MarkerEdgeColor',[1 0.2 0.2])
end

for i = 1:size(fusedMSE,1)
    plot(fusedMSE(i,3),fusedMSE(i,1),'ob','MarkerSize',fusedMSE(i,4)*30,'MarkerEdgeColor',[1 0 0])
end

hold off;
