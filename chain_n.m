clear; clc

%define parameters
n = 6; % # nodes
T = 50; % # periods forward running the risk free rate
risk_ratio = 2.5;
m = .1;
lambda= m; %internsity of distress
eta=m*risk_ratio; %probability of recovery
delta = .95; %time prefrence
gamma = 1; %RRA parameter of CRRA
mu = 0.015; %growth rate in the economy
sigma = 0.1; %volatility of the dividend

%Dividend at each state
epsilon = .2;
x(2) = 1; %x_bar, value of idiosyncratic component at boom state
x(1) = x(2)-epsilon; %x_lowerbar, value of idiosyncratic component at busat state

Delta = chain_net(n); %
%%
H = zeros(n,2^n); %States of the world (H)

for k=1:2^n %state of the world
    str = dec2bin(k-1,n); %generate binary number k-1
    digit = str-'0'; %arrange each digit in a row vector
    H(:,k) = digit'; %define this vector as state of the world k
end



state_pairs01 = []; % 2^n-1 by 2 matrix of state-pairs
counter = 0;
state_pairs01 = NaN(2^(n-1),4);
for dn = 1:floor(n/2)+1
    % dn is a dummy var
    dist_num = dn-1; % dist_num is # distressed
    aux0_H_index = find(sum(H) == dist_num); % find indices of dist_num
    % Find for each state it's equivalent
    if dist_num == n/2
        for i=1:size(aux0_H_index,2)/2
            state_pairs01(i+counter,:) = [aux0_H_index(i),2^n+1-aux0_H_index(i),dist_num,n-dist_num];
        end
    else
        for i=1:size(aux0_H_index,2)
            state_pairs01(i+counter,:) = [aux0_H_index(i),2^n+1-aux0_H_index(i),dist_num,n-dist_num];
        end        
    end
    
    counter = counter+size(aux0_H_index,2);
end

state_pairs01 = sortrows(state_pairs01,1);
%{
*************************************************************************
            Second symmetry: mirror image
*************************************************************************
%}


state_pairs_ = []; %Matrix of state-pairs
for dist_num = 1:floor(n/2) % can use only until n/2, the rest is by 0-1

    state_pairs = []; %Matrix of state-pairs
    
% if n is odd or dist_num not n/2 can use
    aux_H_index = state_pairs01(find(state_pairs01(:,3) == dist_num),1);
% if n is even and dist_num=n/2 we need to use also
if dist_num == n/2
    aux_H_index = [aux_H_index,...
        state_pairs01(find(state_pairs01(:,3) == dist_num),2)];
end
    
    aux = H(:,aux_H_index);
    
%Find pairs of states which are identical up to rotation and mirror
    for i=1:size(aux,2)
        %clc
        %sprintf('%d/%d, dist_num=%d',i,size(aux,2),dist_num)
        flag = 0;
        for j=1:size(aux,2)
            if i < j
                if isequal(aux(:,i)', fliplr(aux(:,j)'))
                    state_pairs = [state_pairs; ...
                        aux_H_index(i),aux_H_index(j)];
                    flag = 1;
                end
            end
        end

        if flag == 0 % the state is unique
            state_pairs = [state_pairs; aux_H_index(i), NaN(1)];
        end
    end
        
        

    
    state_pairs_ = [state_pairs_; unique(state_pairs,'rows')]; %Matrix of state-pairs


end

% remove duplicates
row2remove = [];
for i=1:size(state_pairs_,1)
    [row,col] = find(state_pairs_(:,1) == state_pairs_(i,2));
    row2remove = [row2remove; row];
end

row2remove = unique(row2remove)';
state_pairs_(row2remove,:)=[];
state_pairs_ = sortrows(state_pairs_);


%{
So far we have eliminated rows that are duplicates. But some states on the
left column of state_pairs_ may have an NaN on the right -- the following
part deals with those.
%}

%find indices of unique states that have no match
ind = find(isnan(state_pairs_(:,2))); %ind in state_pairs_
ind_ = state_pairs_(find(isnan(state_pairs_(:,2))),1); %num of state

row2remove = [];
for i = 1:size(ind,1)
%if there is more than one entry of that unique state then remove the row
%with the NaN on the right
    if sum(state_pairs_(:,1) == ind_(i)) > 1 
        row2remove = [row2remove; ind(i)];
    end    
end

state_pairs_(row2remove,:)=[];

%{
we are now left with a list of state_pairs where on the left col we
have unique set of states and on the right col states that correspond to
them in H. But this is done only for half of the matrix -- now we find the
list of complements. Need to be careful of states with dist=n/2
%}

% find the translation of the latter list into their complements
state_pairs__ = [];
for i = 1:size(state_pairs_,1)
	if sum(H(:,state_pairs_(i,1))) ~= n/2
        state_pairs__(i,1) = ...
        state_pairs01(find(state_pairs01(:,1) == state_pairs_(i,1)),2);

    
%some states are unique (have NaN on the right), so then when trying to
%find a complement there would not be one. could just add NaN(1,2) to the
%list of complements - both would work
        if isnan(state_pairs_(i,2)) == 0
            state_pairs__(i,2) = ...
        state_pairs01(find(state_pairs01(:,1) == state_pairs_(i,2)),2);
        else
            state_pairs__(i,2) = NaN(1);
        end
    end
end

state_pairs__(find(state_pairs__(:,1) == 0),:) = [];

state_pairs = [state_pairs_; state_pairs__]; % full list of state pairs
unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n]; % full list of unique states

%arrange unique states by number of distressed
aux = NaN(size(unique_states,1),2);
for k = 1:size(unique_states,1)
    aux(k,:) = [sum(H(:,unique_states(k))),unique_states(k)];
end


[aux_,index]=sort(aux(:,1)); %something is wrong here sortrows perhaps

%use the column indices from sort() to sort all columns of aux.
state_num_dist = aux(index,:)
%unique_states = state_num_dist(:,2)

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    kk=k-1;
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==kk);
end


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
indices = [state_pairs(find(state_pairs(:,1) == state_num_dist(k_dest,2)),2); state_num_dist(k_dest,2)];
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

%%
%{
a = mpower(Pi,1e10);
plot(a(1,:),'-ro')
hold on
plot(sum(a,2),'-r')

axis([1,size(state_num_dist,1),0,1.1])

return
%}


%%

%{
*************************************************************************
ARROW-DEBREU MATRIX
*************************************************************************
%}

dividend = x(2)*(n-state_num_dist(:,1))+x(1)*state_num_dist(:,1);

MU = @(u) u.^(-gamma); %marginal utility function
U = MU(dividend); %define marginal utility vector curly U
D = diag(U);
delta_star = delta - gamma*mu+.5*gamma*(gamma+1)*sigma^2;
%A = delta*Pi.*(MU(dividend)*(1./MU(dividend))')'; %A-D matrix
A = exp(-1+delta_star)*(inv(D)*Pi*D);

y_t = NaN(T,size(state_num_dist,1));
B_t = NaN(T,size(state_num_dist,1));
e = ones(size(state_num_dist,1),1);
R = NaN(T);
for t = 1:T
    %B_t(t,:) = exp(-t*(1-delta_star))*mpower(inv(D)*Pi*D,t)*e;
    B_t(t,:) = mpower(A,t)*e;
    y_t(t,:) = -log(B_t(t,:))/t+1;
    for tau = 1:T
        %R(tau,t) = log((mpower(Pi,tau).*mpower(A,T-tau)*e)/(mpower(A,T)*e))...
        %    +log(mpower(A,tau)*e);
    end
    %y_t(t,:) = (1./(sum(mpower(A,t),2))).^(1/t);
    %y_t(t,:) = (1./(sum(mpower(A,t),2))-1)/t+1;
    %y_t(t,:) = (-log(sum(mpower(A,t),2)))/t+1;    
end


line_style = {'-','--',':','-.'};
%line_color = {'k','b','g','m'}; %,'g','y','m','c'};
line_color = {'k','b','g','m','g','y','m','c'};
Alphabet=char('a'+(1:26)-1)';

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
ring_fig = figure;%('Position', [0, 0, 1365, 690]);
set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]);

set(gcf, 'Visible', 'off')

for k=2:size(state_num_dist,1)-1
       
    %h is subfigure handle
    h(state_num_dist(k,1)) =...
        subplot(fig_dim(2),fig_dim(1),state_num_dist(k,1));
    %p is the curve handle
    p(k) = plot(y_t(:,k)-1);
    hold on
    
%{

    
    
%}
    cyc_counter = cyc_counter+1;
    % Restard legend for every subplot
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
    
    title(sprintf('%s) %d distressed',...
        Alphabet(state_num_dist(k,1)),state_num_dist(k,1)))
    xlabel(sprintf('Maturity (t)'))
    ylabel(sprintf('Yield (R_f)'))

%Legend construction
%number of distressed firms
     leg{cyc_counter} = dec2bin(state_num_dist(k,2)-1,n);
     

    %When reaching the end of cycle 2 restart the counter
    if cyc_counter == unique_state_x_dist(1+state_num_dist(k,1))
        %hold on
        %plot(y_t(:,1)-1,'r');
        leg{cyc_counter+1} = dec2bin(0,n);
        %hold on
        %plot(y_t(:,size(state_num_dist,1))-1,'r--');
        leg{cyc_counter+2} = dec2bin(2^n-1,n);
        cyc_counter = 0;
%Legend is differemt for each n. Create legend for n=4
        legend(leg{1:end},...
            'Location','southeast','Orientation','vertical');
        legend('boxoff');
        leg_record(state_num_dist(k),1:size(leg,2)) = leg;
        %axis([0,T,min(y_t(1,:)-1),max(y_t(1,:)-1)])
        %axis([0,T,min(y_t(1,1:end-1)-1),max(y_t(1,1:end-1)-1)])        
        
    end
    
    
end


%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(ring_fig,'-dpng','-r80',...
    sprintf('ring_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}


%%

width_fig = 700;
height_fig = 600;
x_pos = 100;
y_pos = 50;

% Choose figure to dispaly
subfig = 2;
c=copyobj(h(subfig),get(h(subfig),'Parent'));
chain_subfig = figure;
set(c,'Parent',chain_subfig);
set(chain_subfig, 'Position', [x_pos y_pos width_fig height_fig])
ha=gca;
set(c,'Position','default');
hay=get(ha,'Ylabel');
hax=get(ha,'Xlabel');
NewFontSize=10;
set(hay,'Fontsize',NewFontSize);
set(hax,'Fontsize',NewFontSize);

graph_title = strcat(sprintf('Chain network, n=%d (%d distressed)',n,subfig),...
    sprintf('\nParametrization:'),...
    ' \lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon)...
    );
title(graph_title)


leg = leg_record(subfig,1:unique_state_x_dist(subfig+1));
leg{end+1} = dec2bin(0,n);
leg{end+1} = dec2bin(2^n-1,n);
legend(leg,...
            'Location','southeast','Orientation','vertical');
legend('boxoff');


%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(chain_subfig,'-dpng','-r80',...
sprintf('ring_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
%}






















%{

line_shape = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
chain_fig = figure;
unique_states = [1,2,3,4,6,7,10,8,12,16];
for k=1:10%size(unique_states,1)
    plot(r(:,unique_states(k))-1,line_shape{k}) 
    hold on

end

title(sprintf('Chain network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))

%Legend is differemt for each n. Create legend for n=4
legend('No distress','1 outer distressed','1 inner distressed',...
    '2 distressed: in-out neighb','2 distressed: in-out non-neighb',...
    '2 distressed: in-in','2 distressed: out-out',...
    '3 distressed: out','3 distressed: in','4 distressed',...
    'Location','northeast','Orientation','vertical')
axis([0,50,0,0.13])


%%

cd C:\Users\Oren\Documents\MATLAB\Network
print(chain_fig,'-djpeg',...
    sprintf('chain-eta%.2f-lambda%.2f',eta,lambda))




%}