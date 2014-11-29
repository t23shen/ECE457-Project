function [gFbest,outCluster,geoCenter,runtimes] = gaclustering(data, k, pop,Iter)
% data: data points need to be clustered
% k: number of clusters
% pop: population size
% maximum: number of iterations
%   Detailed explanation goes here
    [num, d] = size(data);     % num: number of points; d: dimensions
    F = zeros(1,pop);
    Flist = Inf.*ones(1,pop);  
	nsite=3;    % number of mutation sites
	pc=0.95;    % Crossover probability
	pm=0.2;    % Mutation probability
	nsbit=6;   % String length (bits)
    initpop = cell(1,k);
    count = 0;
    maxfit = 0;
    avgfit = 0;
    runtimes = 0;
    gFbest = Inf;
    nonimproviter = 0;
    centers = cell(k,1);
    outCluster = cell(1,k);
    for p = 1:pop
        for a = 1:k
            initpop{p,a}= randi([0,9],[1,nsbit*d]);
        end
    end
    tic
    % Initialize solution <- initial population
    for t = 1:Iter
    %while gFbest >330
        popnew=init_gen(pop,k,initpop);
        clear clusters
        clusters = cell(pop,k);
        clear sumDistinCluster;
        sumDistinCluster =zeros(pop,k);
        clear F
        F = zeros(1,pop);  
        for i = 1:pop                                          % number of population
            for p = 1:num                                       % number of points
                mindSqr = Inf;
                minIndex = 1;
                for j = 1:k                                     % number of clusters
                    dSqr = (norm(popnew{i,j} - data(p,:)))^2;
                    if (dSqr < mindSqr)
                        mindSqr = dSqr;
                        minIndex = j;
                    end
                end
                clusters{i,minIndex} = [clusters{i,minIndex}, p];
            end
            geoCenter = zeros(k, d);
            for j = 1:k                                         % number of clusters
               for q = 1:length(clusters{i,j})                  % number of points in each cluster
                   geoCenter(j,:) = geoCenter(j,:) + data(clusters{i,j}(q),:);
               end
               geoCenter(j,:) = geoCenter(j,:) ./ length(clusters{i,j});
               
               for q = 1:length(clusters{i,j})                  % number of points in each cluster
                   sumDistinCluster(i,j) = sumDistinCluster(i,j) + (norm(data(clusters{i,j}(q),:) - geoCenter(j,:)))^2;
               end
               sumDistinCluster(i,j) = sumDistinCluster(i,j);
            end

            % calculate the fitness function
            for j = 1:k
                F(i) = F(i) + sumDistinCluster(i,j);
            end
            F(i) = (F(i)/num)/d; 
            
            % update Fbest and pbest
            if (F(i) < Flist(i))
                Flist(i) = F(i);
            end
        end
       % Record the fittest
       [bstFthisIter, bstIndex] = min(Flist);
       if (bstFthisIter < gFbest)
            gFbest = bstFthisIter;
            outCluster = clusters(bstIndex,:);
       end
       avgfit = mean(F);
      
       % Crossover pair
       ii=floor(pop*rand)+1; jj=floor(pop*rand)+1;
       if(F(ii)> gFbest && F(jj)> gFbest)
           count=count+2;
           % Cross over
           for c = 1:k
                [initpop{ii,c},initpop{jj,c}]= crossover(initpop{ii,c},initpop{jj,c});
           end
       end
       % Mutation at n sites
       kk=floor(pop*rand)+1;
       if(F(kk) > avgfit)
            count=count+1;
            for c = 1:k
                initpop{kk,c}= mutate(initpop{kk,c},nsite);
            end
       end
       runtimes = runtimes+1;
%       if(nonimproviter >= 50)
 %          break;
 %      end
    end
cputime = toc

    geoCenter = zeros(k, d);
    for i = 1:length(outCluster)
       for j = 1:length(outCluster{i})
           geoCenter(i,:) = geoCenter(i,:) + data(outCluster{i}(j),:);
       end
       geoCenter(i,:) = geoCenter(i,:) ./ length(outCluster{i});
    end
    
    
    for i = 1:length(outCluster{1})
        plot(data(outCluster{1}(i),1), data(outCluster{1}(i),2), 'rx')
        hold on
    end
    
    for i = 1:length(outCluster{2})
        plot(data(outCluster{2}(i),1), data(outCluster{2}(i),2), 'bx')
        hold on
    end
    
    for i = 1:length(outCluster{3})
        plot(data(outCluster{3}(i),1), data(outCluster{3}(i),2), 'gx')
        hold on
    end
    
%     for i = 1:length(outCluster{4})
 %        plot(data(outCluster{4}(i),1), data(outCluster{4}(i),2), 'kx')
%         hold on
%     end
%     for i = 1:length(outCluster{5})
%         plot(data(outCluster{5}(i),1), data(outCluster{5}(i),2), 'yx')
%         hold on
%     end
%     for i = 1:length(outCluster{6})
%         plot(data(outCluster{6}(i),1), data(outCluster{6}(i),2), 'cx')
%         hold on
%     end

    plot(geoCenter(1,1), geoCenter(1,2), 'ro');
    hold on
 
    plot(geoCenter(2,1), geoCenter(2,2), 'bo');
    hold on
    
    plot(geoCenter(3,1), geoCenter(3,2), 'go');
    hold on
	% All the sub functions-------------------------------------
	% generation of the initial population
	function [pop]=init_gen(np,k,initpop)
        pop = cell(np,k);
        for b = 1:np
            for c=1:k
                pop{b,c} = dectodec(initpop{b,c});
            end
        end
    end
    
	function [coord] = dectodec(dec)  
        coord = [];
        psize = length(dec);
        for i= 1:(psize/6)
            x = 0;
            for j = 1:6
                x = x + dec(j+(i-1)*6)*10^(6-j);
            end
            coord = [coord x/(10^4)];
        end
    end
    % Crossover operator
    function [c,d]=crossover(a,b)
        nn=length(a)-1;
        % generating a random crossover point
        cpoint=floor(nn*rand)+1;
        c=[a(1:cpoint) b(cpoint+1:end)];
        d=[b(1:cpoint) a(cpoint+1:end)];
    end

    % Mutatation operator
    function anew=mutate(a,nsite)
        nn=length(a); anew=a;
        for i=1:nsite,
            j=floor(rand*nn)+1;
            anew(j)=randi([0,9],1);
        end
    end
end

