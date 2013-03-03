function [ est_means flows] = meanShiftEstimate( ys, weights, sigmas, distthresh, collapsepointdist)
%MEANSHIFTESTIMATE Summary of this function goes here
%   Detailed explanation goes here
    viz = true;
     n = size(ys,1);
    if nargin < 4
        distthresh = 0.01;
    end
    if nargin < 5
        collapsepointdist = 0.1;
    end
    sigmas = sigmas .* sigmas;
    % clipping function of weights 
    ws = max(0,weights);
    
    for i = 1:n
        Hs(:,:,i) = diag(sigmas);
        invHs(:,:,i) = diag(1./sigmas);
        sqrtinvHs(:,:,i) = diag(sqrt(1./sigmas));
        sqrtinvHsWs(:,:,i) = diag(sqrt(1./sigmas)*ws(i));
    end;
    
    for j = 1:n
       y = ys(j,:);
       if (viz)
          hold on;
          %plot(y(1,1),y(1,2),'b+');
	  flow = [y(1,:)];
       end
       for t= 1:1000
           %Hh = calcHh(y);
           omega2 = 0;
           for k = 1:n
               omega2 = omega2 + sqrtinvHsWs(:,:,k) * exp(-mahalanobisDist(y,ys(k,:),invHs(:,:,k))/2);
           end
           for k = 1:n
               omega1 = sqrtinvHsWs(:,:,k) * exp(-mahalanobisDist(y,ys(k,:),invHs(:,:,k))/2);
               omiInvHs(:,:,k) = omega1 / omega2 * invHs(:,:,i);
           end
           Hh = zeros(size(Hs(:,:,1)));
           for k = 1:n
               Hh = Hh + omiInvHs(:,:,k);
           end
           Hh = inv(Hh);
           
           sumpart = zeros(size(y))';
           
           for i = 1:n
                sumpart = sumpart + omiInvHs(:,:,i) * ys(i,:)';
           end
           ym = Hh * sumpart;
           
           d = norm(ym'-y);
           y = ym';
           if (viz)
              %plot(y(1,1),y(1,2),'b+');
	      flow = [flow; y(1,:)];
           end
           if ( d < distthresh )
               break;
           end
       end
       flows{j} = flow;
       %plot(flow(:,1),flow(:,2),'b');
       %plot(flow(:,1),flow(:,2),'r+','Markersize',8);
       est_means(j,:) = y;
    end
    
    dthresh = 0.1;
    stripped(1,:) = est_means(1,:);
    for i = 1:size(est_means,1)
       m = est_means(i,:);
       kind = size(stripped,1)+1;
       for k = 1:size(stripped,1)
	  v = stripped(k,:);
          dist = norm(v-m);
          if ( dist < dthresh )
             kind = k;
             break;
          end
       end
       stripped(kind,:) = m;
    end
  est_means = stripped;
end

