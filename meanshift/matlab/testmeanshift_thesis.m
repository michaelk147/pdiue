clear all;

randn('seed',5)
dim = 2;
N = 30;


mu = zeros(1,dim) + [0 0];
mu2 = mu + 2 /1;
sigma = 0.1;


S1 = eye(dim) * sigma; 
ps = mvnrnd(mu,S1,N/2);

S2 = eye(dim) * sigma; 
S2(dim,dim) = sigma/2;
ps = [ps; mvnrnd(mu2,S2,N/2)];
ps = [ps; [0 3] ];

N = size(ps,1);




mus = [mu;mu2];



maxs = [min(ps) max(ps)];

extraspace = [maxs(1,3) - maxs(1,1);maxs(1,4) - maxs(1,2)] * 0.1;
axisspecs = [maxs(1,1)-extraspace(1) maxs(1,3)+extraspace(1) maxs(1,2)-extraspace(2) maxs(1,4)+extraspace(2)];





% test meanShiftEstimate
esigma = 0.31;
esigmas = [ esigma esigma ];
useweights = false;
if (useweights )
weights = normrnd(2,1,1,N);
else
weights = ones(1,N);
end
[em flows] = meanShiftEstimate(ps,weights,esigmas,0.0001,0.001);

disweights = weights*10;
disweights = max(0.01,disweights);



markerscale = 3.0;
gfxlinewd  = 4;
for f = 1:2

	figure(f);
	hold off;
	plot(ps(:,1),ps(:,2),'+b;data;','MarkerSize',8*markerscale, 'linewidth', gfxlinewd);
	axis(axisspecs);
	hold on;

	for i = 1:size(flows,2)
           if (f == 1) break end
	   flow = flows{i};
	   plot(flow(:,1),flow(:,2),'b', 'linewidth', gfxlinewd);
	  % plot(flow([2,size(flow,1)],1),flow([2,size(flow,1)],2),'b<','Markersize',8);
	end

 	if (f > 0 )
	  plot(mus(:,1),mus(:,2),'og;real means;','MarkerSize',15*markerscale, 'linewidth', gfxlinewd);
	  plot(em(:,1),em(:,2),'rx;estimates;','MarkerSize',10*markerscale, 'linewidth', gfxlinewd);
	end

	if ( useweights )
  	  for i = 1:size(weights,2)
		 plot(ps(i,1),ps(i,2),'ob','MarkerSize',disweights(1,i),'MarkerEdgeColor',[0.7 0.7 1])
	  end
	end
	legend("location","southeast");
	legend("hide");
	legoff = 0.1;
	range = [(axisspecs(2)-axisspecs(1)) (axisspecs(3)-axisspecs(4))  ];
	legrb = [axisspecs(2) axisspecs(3)] - range*legoff;
	txtrb = [axisspecs(2) axisspecs(3)] - range*legoff - [range(1)*0.05 0];
	ydisplacement = range(2)*0.08;
	

		plot(legrb(1),legrb(2)		          ,'xr','MarkerSize',10*markerscale, 'linewidth', gfxlinewd);
		text(txtrb(1),txtrb(2)			  ,'Sch{\"{a}}tzung','horizontalalignment', 'right');
		plot(legrb(1),legrb(2) - ydisplacement * 1,'og','MarkerSize',15*markerscale, 'linewidth', gfxlinewd);
		text(txtrb(1),txtrb(2) - ydisplacement * 1,'Echter Mittelwert','horizontalalignment', 'right');
		plot(legrb(1),legrb(2) - ydisplacement * 2,'+b','MarkerSize',8*markerscale, 'linewidth', gfxlinewd);
		text(txtrb(1),txtrb(2) - ydisplacement * 2,'Datenpunkt','horizontalalignment', 'right');
		if ( useweights ) 
			plot(legrb(1),legrb(2) - ydisplacement * 3,'ob','MarkerSize',10*markerscale, 'linewidth', gfxlinewd);
			text(txtrb(1),txtrb(2) - ydisplacement * 3,'Gewicht','horizontalalignment', 'right');
		end
	
	if (f == 2)
		plot([legrb(1) - range(1)*1/80; legrb(1) + range(1)*1/80],[legrb(2) - ydisplacement * 3;legrb(2) - ydisplacement * 3],'b', 'linewidth', gfxlinewd);
		text(txtrb(1),txtrb(2) - ydisplacement * 3,'Iterationsverlauf','horizontalalignment', 'right');
	end
	%axis("equal");
	if ( f == 1 )
	  print -dtex -color -F:10 'meanshiftexaus1.tex'
	else
	  print -dtex -color -F:10 'meanshiftexaus2.tex'
	end
end
