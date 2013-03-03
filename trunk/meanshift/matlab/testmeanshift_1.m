clear all;

randn('seed',5)
dim = 2;
N = 30;

mu = zeros(1,dim) + [0 0];
mu2 = mu +2;
sigma = 0.2;


S1 = eye(dim) * sigma; 
ps = mvnrnd(mu,S1,N/2);

S2 = eye(dim) * sigma; 
S2(dim,dim) = sigma/2;
ps = [ps; mvnrnd(mu2,S2,N/2)];

N = size(ps,1);


figure(1);
hold off;
plot(ps(:,1),ps(:,2),'+b;data;','MarkerSize',8);

hold on;
mus = [mu;mu2];

plot(mus(:,1),mus(:,2),'og;real means;','MarkerSize',15);
%plot(mu2(:,1),mu2(:,2),'og','MarkerSize',15);

maxs = [min(ps) max(ps)];

extraspace = [maxs(1,3) - maxs(1,1);maxs(1,4) - maxs(1,2)] * 0.1;
axis([maxs(1,1)-extraspace(1) maxs(1,3)+extraspace(1) maxs(1,2)-extraspace(2) maxs(1,4)+extraspace(2)]);


% test meanShiftEstimate
esigma = 0.5;
esigmas = [ esigma esigma ];
%weights = normrnd(1,0.5,1,N);
weights = ones(1,N);
em = meanShiftEstimate(ps,weights,esigmas,0.001,0.1);

disweights = weights*10;
disweights = max(0.01,disweights);

for i = 1:size(weights,2)
    % plot(ps(i,1),ps(i,2),'ob','MarkerSize',disweights(1,i),'MarkerEdgeColor',[0.7 0.7 1])
end
plot(em(:,1),em(:,2),'rx;estimates;','MarkerSize',10);
legend("location","southeast");
print -deps -color foo.eps
print -dtex -color foo.tex
