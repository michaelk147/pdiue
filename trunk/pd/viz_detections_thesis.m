clear all;
close all;
dim = 3;


pddebugdata
pddebugdata_ssp


N = size(ps,1);

markerscale = 1.0;
gfxlinewd  = 3;
ps = [ps(:,1) 480 .- ps(:,2) ps(:,3)];
fusedMSE = [fusedMSE(:,1) 480 .- fusedMSE(:,2) fusedMSE(:,3)];

figure(1);
plot3(ps(:,1),ps(:,3),ps(:,2),'+b','MarkerSize', 10 * markerscale, 'linewidth', gfxlinewd);

hold on;

useweights = true;
if ( useweights )
  for i = 1:size(weights,2)
	% plot3(ps(i,1),ps(i,3),ps(i,2),'ob','MarkerSize',weights(1,i)*30)
  end
end

plot3(fusedMSE(:,1),fusedMSE(:,3),fusedMSE(:,2),'xr','MarkerSize', 35* markerscale, 'linewidth', gfxlinewd+2);


%axis([0 640 0 .8 0 480],"ji")


view(-23,35);

hold off;

print -dtex -color -F:10 'vizdetections.tex'
