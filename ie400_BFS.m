%Read data
blocks = readmatrix("blocks.txt");
points = readmatrix("points.txt");

idxs = nchoosek(1:size(points,1),2);
paths = cell(size(points,1),1);
pathLengths = zeros(length(idxs),1);
counter = 1;
for i = 1:size(points,1)
    for j = i+1:size(points,1)
        [~,paths{counter},pathLengths(counter)] = findPath(points(i,:),points(j,:),blocks);
        counter = counter + 1;
    end
end

idxs(~pathLengths,:)= [];
paths(~pathLengths) = [];
pathLengths(~pathLengths) = [];
dist = pathLengths;
lendist = length(dist);

Aeq = spones(1:length(idxs)); % Adds up the number of trips
beq = size(points,1);


Aeq = [Aeq;spalloc(size(points,1),length(idxs),size(points,1)*(size(points,1)-1))]; % allocate a sparse matrix
for ii = 1:size(points,1)
    whichIdxs = (idxs == ii); % find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
    Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
end
beq = [beq; 2*ones(size(points,1),1)];

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

tours = detectSubtours(x_tsp,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,size(points,1))]; % a guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1)+1; % Counter for indexing
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
    
    
    % How many subtours this time?
    tours = detectSubtours(x_tsp,idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
end

optimizedPaths = paths(find(x_tsp));
%disp(size(optimizedPaths));
fprintf('The absolute gap of the output is:');
disp(output.absolutegap);
fprintf('The optimal cost is: ');
format shortG
disp(costopt);
figure;
hold on

for i = 1:size(points,1)
    xs = [];
    ys = [];
    currpath = optimizedPaths{i};
    
    for j = 1:size(currpath,1)
        xs = [xs currpath(j,1)];
        ys = [ys currpath(j,2)];
        legend
        plot(xs,ys);
    end
end
plot(points(:,1), points(:,2), 'o');
for i = 1:size(blocks,1)
    rectangle('Position', blocks(i,:), 'EdgeColor', 'b', 'FaceColor', 'r', 'LineWidth', 2);
end


function [hasPath,path,pathLength] = findPath(S,F,blocks)
    import java.util.LinkedList
    q = LinkedList();
    visitedMat = zeros(abs(S-F)+1);
    previousMat = cell(abs(S-F)+1);
    q.add(S);
    index = calculateOriginIndex(S,F,S);
    visitedMat(index(1),index(2)) = 1;
    while(~q.isEmpty())
        N = q.remove()';
        if(N == F)
            break;
        end
        possMoves = [N + [0 1] ; N + [1 0] ; N + [0 -1] ; N + [-1 0]];
        for i = 1:size(possMoves,1)
            index = calculateOriginIndex(S,F,possMoves(i,:));
            
            if(~isBlocked(blocks,possMoves(i,:)) && isInBoundary(S,F,possMoves(i,:)) && ~visitedMat(index(1),index(2)))
                
                q.add(possMoves(i,:));
                visitedMat(index(1),index(2)) = 1;
                previousMat{index(1),index(2)} = N;
            end
        end
    end
   [hasPath,pathLength,path] = findPathHelper(S,F,previousMat);
end

function x = isBlocked(blocks,N)
    for i = 1:size(blocks,1)
        xmin = blocks(i,1);
        xmax = blocks(i,1) + blocks(i,3);
        ymin = blocks(i,2);
        ymax = blocks(i,2) + blocks(i,4);

        if(xmin <= N(1) && N(1) <= xmax && ymin <= N(2) && N(2) <= ymax)
            x = 1;
            return
        end
    end
    x = 0;
end
function x = isInBoundary(S,F,N)
    xmin = min(S(1),F(1));
    xmax = max(S(1),F(1));

    ymin = min(S(2),F(2));
    ymax = max(S(2),F(2));
    if(xmin <= N(1) && N(1) <= xmax && ymin <= N(2) && N(2) <= ymax)
        x = 1;
    else
        x = 0;
    end
end

function index = calculateOriginIndex(S,F,N)
    xmin = min(S(1),F(1));
    ymin = min(S(2),F(2));
    index = N -[xmin ymin]+[1 1];
end

function [hasPath,pathLength,path] = findPathHelper(S,F,previousMat)
    path = [];
    pathLength = 0;
    hasPath = 1;
    if (S == F)
        path = [path ; S];
        return
    end
    index = calculateOriginIndex(S,F,F);
    if(isempty(previousMat{index(1),index(2)}))
        hasPath = 0;
        return
    end
    path = [path ; F];
    desiredDistance = sum(abs(S-F));
    while(~isempty(previousMat{index(1),index(2)}))
        N = previousMat{index(1),index(2)};
        path = [path ; N];
        pathLength = pathLength + 1;
        index = calculateOriginIndex(S,F,N);
        if(pathLength > desiredDistance)
            path = [];
            pathLength = 0;
            hasPath = 0;
            break;
        end
    end
end

function subTours = detectSubtours(x,idxs)
% Returns a cell array of subtours. The first subtour is the first row of x, etc.

%   Copyright 2014 The MathWorks, Inc. 

x = round(x); % correct for not-exactly integers
r = find(x); % indices of the trips that exist in the solution
substuff = idxs(r,:); % the collection of node pairs in the solution
unvisited = ones(length(r),1); % keep track of places not yet visited
curr = 1; % subtour we are evaluating
startour = find(unvisited,1); % first unvisited trip
    while ~isempty(startour)
        home = substuff(startour,1); % starting point of subtour
        nextpt = substuff(startour,2); % next point of tour
        visited = nextpt; unvisited(startour) = 0; % update unvisited points
        while nextpt ~= home
            % Find the other trips that starts at nextpt
            [srow,scol] = find(substuff == nextpt);
            % Find just the new trip
            trow = srow(srow ~= startour);
            scol = 3-scol(trow == srow); % turn 1 into 2 and 2 into 1
            startour = trow; % the new place on the subtour
            nextpt = substuff(startour,scol); % the point not where we came from
            visited = [visited,nextpt]; % update nodes on the subtour
            unvisited(startour) = 0; % update unvisited
        end
        subTours{curr} = visited; % store in cell array
        curr = curr + 1; % next subtour
        startour = find(unvisited,1); % first unvisited trip
    end
end
