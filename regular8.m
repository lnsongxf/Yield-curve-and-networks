clear; clc

%define parameters
n = 8; % # nodes
T = 50; % # periods forward running the risk free rate
risk_ratio = 4;
m = .05;
lambda= m; %internsity of distress
eta=m*risk_ratio; %probability of recovery
delta = .05; %time prefrence
gamma = 1; %RRA parameter of CRRA
mu = 0.015; %growth rate in the economy
sigma = 0.1; %volatility of the dividend

%Dividend at each state
epsilon = .2;
x(2) = 1; %x_bar, value of idiosyncratic component at boom state
x(1) = x(2)-epsilon; %x_lowerbar, value of idiosyncratic component at busat state

reg = 4;
Delta = regular_net2powern(n,reg); %reg-regular network with n nodes

H = NaN(n,1+1+3+5+7+5+3+1+1); %total of 27 states
H(:,1) = zeros(n,1);
H(:,2) = [1; zeros(n-1,1)];
H(:,3:6) = [[1;1;zeros(n-2,1)],[1;0;0;1;0;0;0;0],[1;0;1;0;0;0;0;0],[1;0;0;0;1;0;0;0]];
H(:,7:11) = [[1;1;1;0;0;0;0;0],[1;0;1;1;0;0;0;0],[1;0;1;0;0;1;0;0],[1;0;0;0;1;1;0;0],[1;0;1;0;1;0;0;0]];
H(:,12:19) = [[1;1;1;1;0;0;0;0],[1;0;1;0;0;1;0;1],[1;1;1;0;0;1;0;0],[1;0;1;1;0;0;0;1],[1;1;1;0;1;0;0;0],[1;0;1;0;1;0;0;1],[1;0;1;0;1;0;1;0],[1;0;0;1;1;0;0;1]];
H(:,20:24) = [fliplr(1-H(:,6:10))];
H(:,25:27) = [fliplr(1-H(:,3:5))];
H(:,28) = [fliplr(1-H(:,2))];
H(:,29) = [fliplr(1-H(:,1))];
%{
H = [zeros(n,1),[1,zeros(1,n-1)]',...
    [1,1,zeros(1,n-2)]',[1,0,0,1,0,0,0,0]',[1,0,1,0,0,0,0,0]',...
    [1,1,1,0,0,0,0,0]',[1,0,1,1,0,0,0,0]',[1,0,1,0,0,1,0,0]',[1,0,0,0,1,1,0,0]',[1,0,1,0,1,0,0,0]',...
    [1,1,1,1,0,0,0,0]',[1,0,1,0,0,1,0,1]',[1,1,1,0,0,1,0,0]',[1,0,1,1,0,0,0,1]',[1,1,1,0,1,0,0,0]',[1,0,1,0,1,0,0,1]',[1,0,1,0,1,0,1,0]'];
H = [H, 1-H];
%}
H = fliplr(H')';
num_states_per_unique = [1,8,8,8,8,4,8,16,8,16,8,8,8,8,8,16,16,2,4];
dist_num = [0,1,2*ones(1,4),3*ones(1,5),4*ones(1,8)];

num_states_per_unique = [num_states_per_unique,fliplr(num_states_per_unique(1:11))];
dist_num = [dist_num, n-fliplr(dist_num(1:11))];
dim = size(H,2);

counter = 1;
% 0 dist
pairs = [1,NaN(1),dist_num(counter)];
counter = counter+1;
% 1 dist
aux = NaN(num_states_per_unique(counter),1);
for i = 1:num_states_per_unique(2)
    aux(i) = 2^(i-1)+1;
end

pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 2 dist 1 in each, and connected
%run on odd nodes for i, and for each take the indice of even that are
%connected to it for j, then use formulae i j 2^(i-1)+1+2^(j-1)+1-1
%1:2 1+3-1=5 1:8 2+129-1=130

aux = NaN(num_states_per_unique(counter)*2,1);
%eventually what this does is 18 12 23 34 45 ... i i+1
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1);
    j = [i*ones(2,1), search_on_ind(ind)'];
    aux(2*i-1:2*i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 2 dist 1 in each, and not connected
aux = NaN(num_states_per_unique(counter)*2,1);
%eventually what this does is 14 16 38 36 52 58 etc
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = find(Delta(i,search_on_ind) == 0);
    j = [i*ones(2,1), search_on_ind(ind)'];
    aux(2*i-1:2*i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;


% 2 dist in one subgroup (ring n/2), neighbours
aux = NaN(num_states_per_unique(counter),1);
%eventually what this does is 14 16 38 36 52 58 etc
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group
    ind = find(Delta(i,search_on_ind) == 1); %neigbours
    j = [i*ones(2,1), search_on_ind(ind)'];
    aux(2*i-1:2*i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;


% 2 dist in one subgroup (ring n/2), not neigbours
aux = NaN(num_states_per_unique(counter),1);
%eventually what this does is 14 16 38 36 52 58 etc
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group
    search_on_ind(search_on_ind == i) = []; %don't search on self    
    
    ind = find(Delta(i,search_on_ind) == 0); %no-neigbours
    j = [i, search_on_ind(ind)];
    aux(2*i-1:2*i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 3 dist: 2 next to each other and one in the middle
aux = NaN(num_states_per_unique(counter),1);

for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search other sub-group
    ind = find(Delta(i,search_on_ind) == 1); %neigbours
    j = [i, search_on_ind(ind)];
    aux(i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 3 dist: 2 next to each other and one in the side (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1); %search for neighbours
    j = search_on_ind(ind);
    
    for jj = 1:2
        search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
        ind = find((Delta(i,search_on_ind) == 0)&... search for no neighbour to i
            (Delta(j(jj),search_on_ind) == 1)); % but is a neighbour to j(jj)
        jjj = [i, j(jj), search_on_ind(ind)];
        aux(2*i-jj+1) = sum(2.^(jjj-1),2)+1;
    end
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 3 dist: 2 next to each other and one opposite
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 0);
    j = [i, search_on_ind(ind)];
    aux(i) = sum(2.^(j-1),2)+1;

end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 3 dist: 2 not next to each other and one in the other subgroup (4 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
	search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1); %search for neighbours
    ii = search_on_ind(ind);
    
    for iii = 1:2
        search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
        search_on_ind(search_on_ind == ii(iii)) = []; %don't search on self
        ind = find(Delta(ii(iii),search_on_ind) == 0); %find no neighbours
        j = [i, ii(iii), search_on_ind(ind)];
        aux(2*i-iii+1) = sum(2.^(j-1),2)+1;
    end
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 3 dist: 3 in one subgroup (4 options) and 0 at other (1 option)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = find(~isnan(Delta(i,search_on_ind))); %search all other subgroup
    j = search_on_ind(ind);
	aux(i) = sum(2.^(j-1),2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 2 neighbours in one subgroup (4 options) and 2 at neigbours other with one at middle (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = Delta(i,search_on_ind) == 1; %find neighbours
    j = [i, search_on_ind(ind)];
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = find(Delta(i,search_on_ind) == 1); %find neighbours    
    ii = search_on_ind(ind);

	aux(2*i-1:2*i) = [2.^(ii-1)+sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 2 neighbours in one subgroup (4 options) and 2 at neigbours other with one opposite (2 options)
aux = NaN(num_states_per_unique(counter),1);

for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = Delta(i,search_on_ind) == 0; %find non-neighbours
    j = [i, search_on_ind(ind)];
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = find(Delta(i,search_on_ind) == 1); %find neighbours    
    ii = search_on_ind(ind);

	aux(2*i-1:2*i) = [2.^(ii-1)+sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 2 non-neighbours in one subgroup (2 options) and 2 neigbours at other with one at middle (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = Delta(i,search_on_ind) == 1; %find neighbours
    j = [i, search_on_ind(ind)];
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = find(Delta(i,search_on_ind) == 0); %find non-neighbours
    j = [j, search_on_ind(ind)];

	aux(i) = [sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 2 non-neighbours in one subgroup (4 options) and 2 at neigbours other with one at middle (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = Delta(i,search_on_ind) == 0; %find non-neighbours
    j = [i, search_on_ind(ind)];

    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1); %find neighbours
    jjjj = search_on_ind(ind);
    
    for jj=1:2
        jjj = jjjj(jj);
        search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
        %find jjj neighbours that is not a neighbour of i
        ind = find((Delta(jjj,search_on_ind) == 1)&...
                   (Delta(i,search_on_ind) == 0));
        j = [j, jjjj(jj), search_on_ind(ind)];
        aux(2*i+jj-2) = [sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
        j(3:4) = [];
    end
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 3 neighbours in one subgroup (4 options) and 1 at neigbours sub-g not at opposite (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = ~isnan(Delta(i,search_on_ind)); %search all other subgroup
    j = search_on_ind(ind);
    
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 0); %find neighbours
    ii = search_on_ind(ind);

	aux(2*i-1:2*i) = [2.^(ii-1)+sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 3 neighbours in one subgroup (4 options) and 1 at neigbours sub-g at opposite (2 options)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = ~isnan(Delta(i,search_on_ind)); %search all other subgroup
    j = search_on_ind(ind);
    
    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1); %find neighbours
    ii = search_on_ind(ind);

	aux(2*i-1:2*i) = [2.^(ii-1)+sum(2.^(j-1))+size(j,2)-size(j,2)+1]';
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 4 in one subgroup (1 option) and 0 at other (1 option)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group    
    j = search_on_ind;

	aux(i) = sum(2.^(j-1))+size(j,2)-size(j,2)+1;
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];
counter = counter+1;

% 4 dist: 2 non-neighbours in each subgroup (2 options each)
aux = NaN(num_states_per_unique(counter),1);
for i = 1:n
    search_on_ind = (2-mod(i,2)):2:n; %search on own sub-group
    search_on_ind(search_on_ind == i) = []; %don't search on self
    ind = Delta(i,search_on_ind) == 0; %find non-neighbours
    j = [i, search_on_ind(ind)];

    search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
    ind = find(Delta(i,search_on_ind) == 1); %find neighbours
    ii = search_on_ind(ind);
    
    for iii = 1:2
        iiii = ii(iii);
        search_on_ind = (1+mod(i,2)):2:n; %search on neigbour sub-group
        search_on_ind(search_on_ind == iiii) = []; %don't search on iiii
        ind = find(Delta(iiii, search_on_ind) == 0); %find non-neighbours
        j = [j, iiii, search_on_ind(ind)];
        aux(i) = sum(2.^(j-1))+size(j,2)-size(j,2)+1;        
        j(3:4) = [];
    end
end
aux = unique(aux);
pairs = [pairs; [aux(1)*ones(size(aux,1)-1,1),aux(2:end),ones(size(aux,1)-1,1)*dist_num(counter)]];

pairs = [pairs; fliplr([2^n+1-pairs(pairs(:,3) ~= 4,1:2),n-pairs(pairs(:,3) ~= 4,3)]')'];
state_pairs = unique(pairs,'rows');
state_num_dist = fliplr(unique(pairs(:,[1,3]),'rows'));
state_num_dist = sortrows(state_num_dist);

num_state = size(state_num_dist,1);

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==k-1);
end

dividend = x(2)*(n-state_num_dist(:,1))+x(1)*state_num_dist(:,1);

Pi = NaN(size(state_num_dist,1));
%{
Constructing the transition matrix
The way to do so is to start at a unique state -->
then compute probabilities of transition from that unique state to all
other states - but! do it in groups of equiv states to keep track of 
the probability from one unique state to the other -->
at the end sum up the appropriate probabilities and we have an entry in Pi
%}
for k_origin = 1:size(state_num_dist,1)
    clc
    sprintf('%d/%d',k_origin,size(state_num_dist,1))

str_ = dec2bin(state_num_dist(k_origin,2)-1,n); %generate binary number k-1
digit = str_-'0'; %arrange each digit in a row vector
H_origin = logical(digit'); %define this vector as state of the world

Lambda = lambda*Delta*H_origin;
Eta = eta*ones(n,1);
%Eta = eta*Delta*(1-H_origin); %network effects for recovery

%Natural rate of distress/recovery
%Lambda = lambda*Delta*H_origin+lambda;
%Eta = eta*Delta*(1-H_origin)+eta;

for k_dest = 1:size(state_num_dist,1)

% indices find the numbers of states that are equiv to state_num_dist(k_dest,2)
indices = [state_pairs(state_pairs(:,1) == state_num_dist(k_dest,2),2); state_num_dist(k_dest,2)];
%remove NaN: states that are unique have only 1 destination to count
indices = indices(isnan(indices)==0);  
aux_ = NaN(size(indices,1),1);

% then run over each of those states and find the prob to get to that state
% from k_origin, and then sum these probs up.
for ind = 1:size(indices,1)

str_ = dec2bin(indices(ind)-1,n);
digit = str_-'0';
H_dest = logical(digit');
aux = NaN(n,1);

for i = 1:n
aux(i) = (1-Eta(i))*H_origin(i)*H_dest(i)+Eta(i)*H_origin(i)*(1-H_dest(i))+...
    Lambda(i)*(1-H_origin(i))*H_dest(i)+(1-Lambda(i))*(1-H_origin(i))*(1-H_dest(i));
%{
aux(i) = (1-Eta(i,k_origin))*H(i,k_origin)*H(i,k_dest)+Eta(i,k_origin)*H(i,k_origin)*(1-H(i,k_dest))+...
         Lambda(i,k_origin)*(1-H(i,k_origin))*H(i,k_dest)+(1-Lambda(i,k_origin))*(1-H(i,k_origin))*(1-H(i,k_dest));
%}
end
         aux_(ind) = prod(aux);
end
         
Pi(k_origin,k_dest) = sum(aux_);
        
end
end


%{
*************************************************************************
ARROW-DEBREU MATRIX
*************************************************************************
%}


MU = @(u) u.^(-gamma); %marginal utility function
%alpha_RA = .5;
%MU = @(u) exp(-alpha_RA*u); %marginal utility function
U = MU(dividend); %define marginal utility vector curly U
D = diag(U);
delta_star = delta - gamma*mu+.5*gamma*(gamma+1)*sigma^2;
A = exp(-delta_star)*((D\eye(num_state))*Pi*D);

y_t = NaN(T,num_state);
B_t = NaN(T,num_state);
e = ones(num_state,1);
R = NaN(T,T,num_state);
for t = 1:T
    B_t(t,:) = mpower(A,t)*e;
    y_t(t,:) = -log(B_t(t,:))/t;
    for tau = 1:t-1
        R(t,tau,:) = mpower(Pi,tau)*...
        log(mpower(A,t-tau)*e)-...
        log(mpower(A,t)*e)+log(mpower(A,tau)*e);
        R(t,tau,:) = R(t,tau,:)/tau;
    end
end

%Cov
Return = NaN(num_state,num_state,T,T);
E_Return_cond = NaN(num_state,T,T);
Y = NaN(num_state);
E_Y_cond = NaN(n,T);
Cov = NaN(num_state,T,T);
for i = 1:num_state
    for j = 1:num_state
        Y(i,j) = dividend(i)/dividend(j);
    end
end

for t = 2:T
    clc
    sprintf('%d/%d',t,T)
    for tau = 1:t-1
        aux = mpower(Pi,tau);
        for i = 1:num_state
            for j = 1:num_state
                Return(i,j,tau,t) = B_t(t,i)/B_t(t-tau,j);
            end
            E_Return_cond(i,tau,t) = Pi(i,:)*Return(i,:,tau,t)';
            E_Y_cond(i,tau) = aux(i,:)*Y(i,:)';
            
            Cov(i,tau,t) = ...
            sum(Pi(i,:).*(Return(i,:,tau,t)-E_Return_cond(i,tau,t)).*(Y(i,:)-E_Y_cond(i,tau)));
            
        end
    end
end


%%
%{
**************************************************************************
                                MISC
**************************************************************************
%}

line_style = {':','--','-','-.'};
%line_color = {'k','b','g','m'}; %,'g','y','m','c'};
line_color = {'k','b','g','m','g','y','m','c'};
Alphabet=char('a'+(1:26)-1)';

border = 2e-2;

fs1 = 12;
fs2 = 13;


%%
%{
**************************************************************************
                                FIGURES
                     BOND PREMIA AND YIELD CURVES
**************************************************************************
%}

tau = 1;
t = 50;
highlight = 19;
dist_num = 4;
%choose_states = sort(unique([3, 5, 7, highlight]));
[choose_states, aux] = find(state_num_dist(:,1)==dist_num);
choose_states = choose_states';


graph_title = strcat(sprintf(...
    '                                       4-regular network, n=%d',n),...
    sprintf('\nParametrization:'),...
    ' \lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon)...
    );

figure('Position', [10 50 1350 600],'Visible','on');

h1 = subplot(1,3,1);
shade1 = [min(y_t,[],2), max(y_t,[],2)-min(y_t,[],2)];
a1 = area(shade1,'FaceColor','k','edgealpha',0.1);
alpha .2
a1(1).FaceColor = [ 1 1 1];
hold on

leg = cell(size(choose_states,2)+2,1);
leg{1} = ''; leg{2} = 'Spectrum';

p = matlab.graphics.axis.Axes(num_state); % handles for subplots

for k = 1:size(choose_states,2)
    leg{k+2} = dec2bin(state_num_dist(choose_states(k),2)-1,n);
    if choose_states(k) == highlight
        p(k) = plot(y_t(:,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1},'k*'));
    else
        p(k) = plot(y_t(:,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1}));
        
    end
    %p(k).LineWidth = 1;

    hold on
end

t1 = title('a) Yield curve');
xl = xlabel('T (bond maturity)','FontSize',fs1);
legend(leg,...
    'Location','northeast','Orientation','vertical');
legend('boxoff');

YLim(1) = min(min(y_t(:,choose_states)));
YLim(2) =  max(max(y_t(:,choose_states)));
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])


h2 = subplot(1,3,2);

shade2 = [min(R(:,tau,:),[],3), max(R(:,tau,:),[],3)];
a2 = area(shade2,'FaceColor','k','edgealpha',0.1);
hold on
alpha .2

for k = 1:size(choose_states,2)
    if choose_states(k) == highlight
        p(k) = plot(R(:,tau,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1},'k*'));
    else
        p(k) = plot(R(:,tau,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1}));
    end
    hold on
end

t2 = title(strcat(sprintf('b) Bond premium,\nkeeping holding period'), ' \tau=',...
    sprintf('%d',tau)));
xlabel('T (underlying bond maturity)','FontSize',fs1)
YLim(1) = min(min(R(:,tau,:)));
YLim(2) = max(max(R(:,tau,:)));
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])


h3 = subplot(1,3,3);

shade = [min(R(t,:,:),[],3)', max(R(t,:,:),[],3)'];
a3 = area(shade,'FaceColor','k','edgealpha',0.1);
hold on
alpha .2


for k = 1:size(choose_states,2)
    if choose_states(k) == highlight
        p(k) = plot(R(t,:,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1},'k*'));
    else
        p(k) = plot(R(t,:,choose_states(k)),strcat(line_style{mod(k-1,size(line_style,2))+1}));
    end
    hold on
end

t3 = title(strcat(sprintf('c) Bond premium,\nkeeping bond maturity T='),...
    sprintf('%d',T)));

xlabel('\tau (holding period)','FontSize',fs1)
YLim(1) = min(min(R(t,:,:)));
    YLim(2) = max(max(R(t,:,:)));
%get(h2,'YLim')
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])

fig_height = .55; push_up = .05;
set(h1,'Position',get(h1,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h2,'Position',get(h2,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h3,'Position',get(h3,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])



set(t1,'Position',get(t1,'Position').*[1 1 1]+[0 -.0359 0])
set(t2,'Position',get(t2,'Position').*[1 1 1]+[0 0 0])
set(t3,'Position',get(t3,'Position').*[1 1 1]+[0 0 0])

set(t1,'FontSize',fs2)
set(t2,'FontSize',fs2)
set(t3,'FontSize',fs2)



text_pos = get(xl,'Position').*[1 1 0]+[0 -.001 0];
text_pos = text_pos(1:2);
text(-97, -0.00018, graph_title,'FontSize',16);


return





%%
%{
**************************************************************************
                                FIGURES
                            Single picture
**************************************************************************
%}

newhandle = subplot(1,3,3,'Visible', 'off');
pic_pos = get(newhandle,'Position');
r_pic = [4,4]';

%clc
h = figure('Position', [50 100 1024 280],'Visible', 'off');
h1 = subplot(1,3,1);
tau = 2;
highlight = 4; %num_state
%set(gcf, 'Visible', 'on')
for k = 1:num_state
plot(R(:,tau,k))
hold on
end
plot(R(:,tau,highlight),'r','LineWidth',3)
title(sprintf('tau = %d',tau))
%axis([0 50 -1e-4 5e-4])


subplot(1,3,2)
for k = 1:num_state
plot(y_t(:,k))
hold on
end
plot(y_t(:,highlight),'r','LineWidth',3)

reg = 2;
s_number = state_num_dist(highlight,2);
pic = k_reg_pic(n,reg,s_number,r_pic);
c=copyobj(pic,get(h1,'Parent'));
set(c,'Position',pic_pos);
axis off

%axis([0 50 -5 .5])

%%
%{
**************************************************************************
                                FIGURES
                       Conditional covariance...
**************************************************************************
%}

highlight = 19; %show figure for number of distressed s_number
s_number = state_num_dist(highlight,2);


graph_title1 = strcat(sprintf(...
    'Ring network, n=%d, %d distressed',n,state_num_dist(highlight,1)),...
    sprintf('\nParametrization:'),...
    ' \lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon)...
    );


newhandle = subplot(1,2,2,'Visible', 'off');
pic_pos = get(newhandle,'Position');


h = figure('Position', [50 100 1500 600],'Visible', 'on');
h1 = subplot(1,3,1);

aux = NaN(T);
aux(:,:) = Cov(highlight,:,:);
contourf(aux)

xlabel('Bond maturity T')
ylabel('Holding period \tau')
title(graph_title1)
colorbar


reg = 2;
r_pic = [4,4]';
pic = k_reg_pic(n,reg,s_number,r_pic);
c=copyobj(pic,get(h1,'Parent'));
set(c,'Position',pic_pos);
axis off

set(h1,'Position',get(h1,'Position').*[1 1 1.5 .85]+[0 0 0 0])
set(c,'Position',get(c,'Position').*[1 1 1 .9]+[0 0 0 0])


return


%%

%{
**************************************************************************
                                FIGURES
                       Spectral gap with a picture
**************************************************************************
%}

s_number = 4; %show figure for number of distressed s_number
figure('Position', [50 100 1300 400],'Visible','off');
subplot(1,2,1)
contourf(R(:,:,s_number))
colorbar
xlabel('Holding period, \tau')
ylabel('t-period bond')

h1 = subplot(1,2,2);
axis off
reg = 2;
s_number = state_num_dist(highlight,2);
pic = k_reg_pic(n,reg,s_number,r_pic);
axis off
c=copyobj(pic,get(h1,'Parent'));
set(c,'Position',pic_pos);

return
%%
%{
**************************************************************************
                                FIGURES
                               With pics
**************************************************************************
%}
newhandle = subplot(1,3,3,'Visible', 'off');
pic_pos = get(newhandle,'Position');

%clc
h = figure('Position', [50 100 1024 280]);
h1 = subplot(1,3,1);
t = 50;
highlight = 2;
%set(gcf, 'Visible', 'on')
for k = 1:8
plot(R(t,:,k))
hold on
end
plot(R(t,:,highlight),'r','LineWidth',3)
title(graph_title)



subplot(1,3,2)
for k = 1:8
plot(y_t(:,k))
hold on
end
plot(y_t(:,highlight),'r','LineWidth',3)
%title(direct)

s_number = state_num_dist(highlight,2);
r_pic = [4,3]';
pic = k_reg_pic(n,reg,s_number,r_pic);
c=copyobj(pic,get(h1,'Parent'));
set(c,'Position',pic_pos);





%%
%{
**************************************************************************
                                FIGURES
    All yield curves in one figure, grouped by number of distressed
**************************************************************************
%}

%legend record
leg_record = cell(n-1,max(unique_state_x_dist));


%there are n-1 figures, each depicts 1,...,n-1 dist nodes curves for unique
%states, with 0,n attached to each for comparison

%need to arrange n-1 subfigures in a rectangle
k = 1;
fig_dim = [1,1];

while k < n-1

    
    if fig_dim(1) == fig_dim(2)
        fig_dim = [fig_dim(1)+1,fig_dim(2)];
    else
        fig_dim = [fig_dim(1),fig_dim(2)+1];
    end
    
    k = prod(fig_dim);
    
end

cyc_counter = 0;
ring_fig = figure('Visible', 'off');%('Position', [0, 0, 1365, 690]);
set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]);

h = matlab.graphics.axis.Axes(n-1); % handles for subplots
p = matlab.graphics.axis.Axes(size(state_num_dist,1)); % handles for subplots

for k=2:size(state_num_dist,1)-1
       
    %h is subfigure handle
    h(state_num_dist(k,1)) =...
        subplot(fig_dim(2),fig_dim(1),state_num_dist(k,1));
    %p is the curve handle
    p(k) = plot(y_t(:,k));
    hold on

    cyc_counter = cyc_counter+1;
    % Restart legend for every subplot
    if cyc_counter == 1
        leg = cell(1,unique_state_x_dist(1+state_num_dist(k,1))+2);
    end
    
    cyc_1_size = size(line_color,2);
    cyc_2_size = size(line_style,2);
   
    choose_color = cyc_counter-cyc_1_size*...
        floor((cyc_counter-1)/cyc_1_size);
    choose_line = choose_color-cyc_2_size*...
        floor((choose_color-1)/cyc_2_size);
    %choose_line = 1+floor((cyc_counter-1)/cyc_1_size);
    %[state_num_dist(k,1),cyc_counter,choose_color,choose_line]
    
    p(k).Color = line_color{choose_color};
    p(k).LineStyle = line_style{choose_line};
    
    %title(sprintf('%s) %d distressed',...
    %    Alphabet(state_num_dist(k,1)),state_num_dist(k,1)))
    xlabel(sprintf('Maturity (t)'))
    ylabel(sprintf('Yield (R_f)'))

%Legend construction
%number of distressed firms
     leg{cyc_counter} = dec2bin(state_num_dist(k,2)-1,n);
     

    %When reaching the end of cycle 2 restart the counter
    if cyc_counter == unique_state_x_dist(1+state_num_dist(k,1))
        hold on
        p(1) = plot(y_t(:,1),'r');
        leg{cyc_counter+1} = dec2bin(0,n);
        hold on
        p(size(state_num_dist,1)) = plot(y_t(:,size(state_num_dist,1)),'r--');
        leg{cyc_counter+2} = dec2bin(2^n-1,n);
        cyc_counter = 0;
%Legend is differemt for each n. Create legend for n=4
        legend(leg,...
            'Location','southeast','Orientation','vertical');
        legend('boxoff');
        leg_record(state_num_dist(k),1:size(leg,2)) = leg;
        %axis([0,T,min(y_t(1,:)-1),max(y_t(1,:)-1)])
        axis([0,T,min(y_t(1,1:end-1)),max(y_t(1,1:end-1))])        
        
    end
    
end


%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(ring_fig,'-dpng','-r80',...
    sprintf('ring_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}


%%

%{
**************************************************************************
                                FIGURES
    All yield curves in one figure, grouped by number of distressed:
        choose one subplot to display as a figure on its own
                              with pictures
**************************************************************************
%}

clc

subfig = 4; % Choose figure to dispaly

tau = 1;
t = T;
highlight = 19;
dist_num = subfig;
%choose_states = sort(unique([3, 5, 7, highlight]));
[choose_states, aux] = find(state_num_dist(:,1)==dist_num);
choose_states = choose_states';

if unique_state_x_dist(subfig+1) <= 9
    aux = 5;
elseif unique_state_x_dist(subfig+1) > 9 && unique_state_x_dist(subfig+1) <= 11
    aux = 6;
elseif unique_state_x_dist(subfig+1) > 11 && unique_state_x_dist(subfig+1) <= 13
    aux = 7;
else
    disp('Error: too many unique states to represent in the figure')
    return
end
newhandle = get(subplot(aux,aux,(aux-2)*aux+1),'Position');

YLim(1) = min(min(y_t(:,choose_states)));
YLim(2) =  max(max(y_t(:,choose_states)));

height_fig = 512+75;
width_fig = 1.12*height_fig;
x_pos = 100;
y_pos = 75;


reg_subfig = figure;
set(reg_subfig, 'Position', [x_pos y_pos width_fig height_fig]);

c=copyobj(h(subfig),get(h(subfig),'Parent'));
%newhandle = get(subplot(5,5,(5-2)*5+1),'Position');
set(c,'Parent',reg_subfig);
set(c,'Position',newhandle.*[1,1,aux,aux]);
axis(c,[0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])
%axis(c,[0 50 0.04 0.055])



ha=gca;
hay=get(ha,'Ylabel');
hax=get(ha,'Xlabel');
NewFontSize=10;
set(hay,'Fontsize',NewFontSize);
set(hax,'Fontsize',NewFontSize);

%subplot(1,3,1)

graph_title = strcat(sprintf('%d-regular network, n=%d (%d distressed);',reg,n,subfig),...
    sprintf(' Parametrization:\n --'),...
    '\lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon),' --');
title(c,graph_title);%,'Position',[1.2 2.35]);

leg = leg_record(subfig,1:unique_state_x_dist(subfig+1));
leg{end+1} = dec2bin(0,n);
leg{end+1} = dec2bin(2^n-1,n);
legend(c,leg,...
            'Location','northeast','Orientation','vertical');
legend(c,'boxoff');

s_number = state_num_dist(state_num_dist(:,1)==subfig,2);


if unique_state_x_dist(subfig+1) == 1
    pic_position = 10;
    addition = [.015,-.1,0,0; .015,-.1,0,0; .015,-.1,0,0;...
                .015,-.1,0,0; .015,-.1,0,0 ];
    x_mult = 1.1; y_mult = 1.2;
    multiplication = [ 1,1,x_mult,y_mult ];    
elseif unique_state_x_dist(subfig+1) == 4
    pic_position = [5,15,21,23];
    addition = [.015,-.1,0,0; .015,-.1,0,0; .015,-.1,0,0;...
                .015,-.1,0,0; .015,-.1,0,0 ];
    x_mult = 1.1; y_mult = 1.2;
    multiplication = [ 1,1,x_mult,y_mult; 1,1,x_mult,y_mult;...
                       1,1,x_mult,y_mult; 1,1,x_mult,y_mult ];
elseif unique_state_x_dist(subfig+1) == 5
    pic_position = [5,15,21,23,25];
    addition = [.015,-.1,0,0; .015,-.1,0,0; .015,-.1,0,0;...
                .015,-.1,0,0; .015,-.1,0,0 ];
    x_mult = 1.1; y_mult = 1.2;            
    multiplication = [ 1,1,x_mult,y_mult; 1,1,x_mult,y_mult; 1,1,x_mult,y_mult;...
                       1,1,x_mult,y_mult; 1,1,x_mult,y_mult; ];
elseif unique_state_x_dist(subfig+1) == 8
    pic_position = [5,10,15,20,21,22,23,24];
    addition = [.025,-.05,0,0; .025,-.075,0,0; .025,-.1,0,0; .025,-.125,0,0;...
                -.025,-.1,0,0; 0,-.1,0,0; .025,-.1,0,0; .05,-.1,0,0 ];
            scale = 1;
    x_mult = 1.1*scale; y_mult = 1.2*scale;
    multiplication = [ 1,1,x_mult,y_mult; 1,1,x_mult,y_mult; 1,1,x_mult,y_mult; 1,1,x_mult,y_mult;...
                       1,1,x_mult,y_mult; 1,1,x_mult,y_mult; 1,1,x_mult,y_mult; 1,1,x_mult,y_mult; ];
end

r_pic = [4,2.85]';
for i = 1:unique_state_x_dist(subfig+1)

    newhandle = subplot(5,5,pic_position(i));
    set(newhandle,'Position',get(newhandle,'Position').*multiplication(i,:)+addition(i,:),...
        'Visible', 'off');
    %set(newhandle)
    pic = k_reg_pic(n,reg,s_number(i),r_pic);
    c=copyobj(pic,get(h(subfig),'Parent'));
    set(c,'Parent',reg_subfig,'Position',get(newhandle,'Position'));
    title(c,leg{i});
    
end


%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(reg_subfig,'-dpng','-r80',...
sprintf('ring_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
%}




