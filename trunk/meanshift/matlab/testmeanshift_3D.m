clear all;

randn('seed',42)
dim = 3;
N = 30;

mu = zeros(1,dim) + [0 0 0];
mu2 = mu +2;
sigma = sqrt(0.5);


S = eye(dim) * sigma;

ps = mvnrnd(mu,S,N/2);

ps = [ps; mvnrnd(mu2,S,N/2)];

ps = [ps; mvnrnd(mu - 2,S,1)];

N = size(ps,1);


figure(1);
hold off;
plot3(ps(:,1),ps(:,2),ps(:,3),'xb');

hold on;
plot3(mu(:,1),mu(:,2),mu(:,3),'og');
plot3(mu2(:,1),mu2(:,2),mu2(:,3),'og');

%maxs = sqrt(sigma)*eye(dim) * 2*sqrt(2*log(2))*4;
maxs = [min(ps) max(ps)];

%extraspace = [maxs(1,3) - maxs(1,1);maxs(1,4) - maxs(1,2)] * 0.1;

%axis([maxs(1,1)-extraspace(1) maxs(1,3)+extraspace(1) maxs(1,2)-extraspace(2) maxs(1,4)+extraspace(2)]);
%axis([mu(1)-maxs(1,1),mu(1) + maxs(1,1),mu(2)-maxs(2,2),mu(2)+maxs(2,2)]);

% test meanShiftEstimate
esigmas = [ sigma sigma sigma ];
weights = normrnd(1,0.5,1,N);
%weights = ones(1,N);
em = meanShiftEstimate(ps,weights,esigmas,0.001,0.01);

disweights = weights*10;
disweights = max(0.01,disweights);

for i = 1:size(weights,2)
    plot3(ps(i,1),ps(i,2),ps(i,3),'ob','MarkerSize',disweights(1,i),'MarkerEdgeColor',[0.7 0.7 1])
end
plot3(em(:,1),em(:,2),em(:,3),'ro');