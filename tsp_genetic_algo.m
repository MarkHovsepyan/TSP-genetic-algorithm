%*************************************************************************%
% Project:  Travelling Salesman Problem Using Genetic Algorithm
% Course: Optimization
% By Mark Hovsepyan
%*************************************************************************%

% function TSP: main function that performs genetic algorithm to solve travelling salesman problem.

function TSP = tsp_genetic_algo()
    
%% Initialization
    
    xy_data            = csvread('armenia_xy.csv'); % for data descrition look for "armenia_xy info"
    xy_coord           = xy_data;
    
    %[8 27; 13 20; 14 29; 35 31; 17 10; 1 42; 21 18; 28 2; 29 21; 33 11; 36 3; 41 26; 43 9; 12 20; 15 59; 40 7]; % City Coordinates
    
    pop_size           = 60; % population size
    gen_num = 100 ;          % number of generations
    [num_of_cities, ~] = size(xy_coord); % get the number of cities
    global_min         = Inf; % positive integers
    dist_matrix        = [];  % point a to b distances/cost
    pop                = zeros(pop_size, num_of_cities); % initialize pop matrix of zeroes
    total_dist         = zeros(1,pop_size); % total distance of the best path
    distHistory = zeros(1,gen_num);   % initialize a matrix of zeroes for ploting
    tmp_pop = zeros(4,num_of_cities); % initialize a matrix of zeroes
    new_pop = zeros(pop_size,num_of_cities); % initialize a matrix of zeroes
    current_gen = 0;
    gen = 0; % iterator for generations
   
    
%main body%    
    pop_size = 4*ceil(pop_size/4);
    nPoints = size(xy_coord,1);
    display(nPoints);
    a = meshgrid(1:nPoints);
    dist_matrix = reshape(sqrt(sum((xy_coord(a,:)-xy_coord(a',:)).^2,2))...
    ,nPoints,nPoints);   % compute distances from each points
    display(dist_matrix);
  
    
    % Generate a population using permutation encoding
    
    pop(1,:) = (1:num_of_cities);
    for k = 2:pop_size
        pop(k,:) = randperm(num_of_cities);
    end
    
    % FITNESS EVALUATION: After initializing the population, fitness value for 
    % each individual is calculated. The fitness values are generated using
    % a fitness function. The function provides largest and smallest values
    % for each of the individuals. If the individual has a larger fitness
    % value then result will be a better solution but if a smaller value is
    % obtained then the obtained solution isn't better.

    figure('Name','Best Solution');
    for gen = 1:gen_num
        % Evaluate each population member (Calculate Total Distance)
        for p = 1:pop_size
            tmp_tot_dist = dist_matrix(pop(p,num_of_cities), pop(p, 1)); % temp variable for total dstance distance
            for k = 2:num_of_cities
                tmp_tot_dist = tmp_tot_dist + dist_matrix(pop(p,k-1),pop(p,k)); % summation of distances of pop(k)
            end
            total_dist(p) = tmp_tot_dist; %select as elite
        end
        
        % Find the least cost if the new computed total distance is less than 
        % the current global_min => the new global minimum distance is the min_dist
        
        [min_dist,index] = min(total_dist);
        distHistory(gen) = min_dist;
        if min_dist < global_min
            global_min = min_dist;
            optimal_route = pop(index, :);
            path = optimal_route([1 : num_of_cities 1]);
            % Plot the Best Route
            %pause(1);
            plot(xy_coord(path, 1), xy_coord(path, 2), 'r.-'); 
            title(sprintf('Distance = %1.4f',min_dist));     
            drawnow;
            current_gen = gen;
        end
        
        % CROSSOVER AND MUTATION: 
        randomOrder = randperm(pop_size);
        for p = 4:4:pop_size
            paths = pop(randomOrder(p-3:p),:); % select parent 1 from the random population
            dists = total_dist(randomOrder(p-3:p)); % select parent 2 from the matrix of elite  population
            [~,idx] = min(dists);       % get the index of the path with the min distance from the elite population
            possible_route = paths(idx,:);  %get the values/permutation encoding
            routeInsertionPoints = sort(ceil(num_of_cities*rand(1,2)));
            I = routeInsertionPoints(1); 
            J = routeInsertionPoints(2);
       
            for k = 1:4                         % Mutate and create a 4 new set of population
                tmp_pop(k,:) = possible_route;  
                display(tmp_pop);
                display(size(tmp_pop));
                
                switch k
                    case 2                   
                        tmp_pop(k,I:J) = tmp_pop(k,J:-1:I); % mutate by flipping
                    case 3 
                        tmp_pop(k,[I J]) = tmp_pop(k,[J I]); % mutate by swapping
                    case 4 
                        tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]); % mutate by sliding 
                    otherwise % default
                end
            end
            new_pop(p-3:p,:) = tmp_pop;
    
        end
        pop = new_pop; % new generation  
        
    end  
    
    %*********************************************************************%
    %% Plotting
    figure('Name','RESULTS','Numbertitle','off');
    title('Best Solution History');
    plot(distHistory,'b','LineWidth',3);
    set(gca,'XLim',[0 gen_num+1],'YLim',[0 1.1*max([1 distHistory])]);
    xlabel('Generation');
    ylabel('Distance');
    
    %% Displaying
    display(optimal_route);
    display(min_dist);
    display(current_gen);
    %*********************************************************************%
end    


