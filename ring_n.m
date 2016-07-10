clear; clc

%define parameters
n = 8; % # nodes
T = 50; % # periods forward running the risk free rate
risk_ratio = 4;
m = .2;
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


Delta = ring_net(n); %Ring network

H = zeros(n,2^n); %States of the world (H)
Lambda = zeros(n,2^n); %collect all lambda_i for each node at each state
Eta = zeros(n,2^n); %collect all lambda_i for each node at each state
%Pi = zeros(2^n,2^n); %markov matrix
dividend = NaN(2^n,1); %total dividend at each state

str = cell(2^n,1);

for k=1:2^n %state of the world
    if mod(k,100) == 0
    clc
    sprintf('%d/%d',k,2^n)    
    end
    
    str{k} = dec2bin(k-1,n); %generate binary number k-1
    digit = str{k}-'0'; %arrange each digit in a row vector
    H(:,k) = logical(digit'); %define this vector as state of the world k
    dividend(k) = x(1)*sum(H(:,k))+x(2)*(n-sum(H(:,k))); %total dividend at k

    Lambda(:,k) = lambda*Delta*H(:,k);
    Eta(:,k) = eta;
    %Eta(:,k) = eta*Delta*(1-H(:,k)); %network effects for recovery
    
    %Natural rate of distress/recovery
    %Lambda(:,k) = lambda*Delta*H(:,k)+lambda;
    %Eta(:,k) = eta*Delta*(1-H(:,k))+eta;
end

H



return
%Find pairs of states which are identical up to rotation and mirror
state_pairs = []; %Matrix of state-pairs
count=0;
for i=1:2^n
    clc
    count = count+1;
    sprintf('%d/%d',count,2^n)
    for j=1:2^n
        if i < j
            for k=1:n-1
                if H(:,j) == circshift(H(:,i),k)
                    state_pairs = [state_pairs; i,j];
                elseif H(:,j) == circshift(fliplr(H(:,i)')',k) %mirror image
                    state_pairs = [state_pairs; i,j];
                end
            end
        end
    end
end

row2remove = [];
for i=1:size(state_pairs,1)
    [row,col] = find(state_pairs(i,2) == state_pairs(:,1));
    row2remove = [row2remove; row];
end
row2remove = row2remove';
state_pairs(row2remove,:)=[];
unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n]; %Unique states of the world

%arrange unique states by number of distressed
aux = NaN(size(unique_states,1),2);
for k = 1:size(unique_states,1)
    aux(k,:) = [sum(H(:,unique_states(k))),unique_states(k)];
end

[aux_,index]=sort(aux(:,1));

%use the column indices from sort() to sort all columns of aux.
state_num_dist = aux(index,:);
unique_states = state_num_dist(:,2);

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    kk=k-1;
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==kk);    
end

unique_states1 = unique_states;

%Constructing the transition matrix

for k_origin = 1:2^n
    for k_dest = 1:2^n

        aux = NaN(n,1);
        for i = 1:n
aux(i) = (1-Eta(i,k_origin))*H(i,k_origin)*H(i,k_dest)+Eta(i,k_origin)*H(i,k_origin)*(1-H(i,k_dest))+...
         Lambda(i,k_origin)*(1-H(i,k_origin))*H(i,k_dest)+(1-Lambda(i,k_origin))*(1-H(i,k_origin))*(1-H(i,k_dest));
        end
        Pi(k_origin,k_dest) = prod(aux);
        
    end
end



%%
%{
*************************************************************************
CENTRALITY
*************************************************************************
%}
clc
tau_horizon = 2; %measure centrality up to tau_horizon periods ahead

[eigvec,eigval] = eig(Pi);
Prob = eigvec(:,2^(n-1));
Joint_Prob = NaN(n);
DC = NaN(n,n,tau_horizon);

%{
Pi_tau = NaN(2^n,2^n,tau_horizon);

for tau = 1:tau_horizon
    Pi_tau(:,:,tau) = mpower(A,tau)
    
end
%}
%different properties, can't be applied...
indices = NaN(size(find(H(1,:)==1),2),n);
for i = 1:n
    indices(:,i) = find(H(i,:)==1);
end


for i = 1:n
    for j = 1:n
        
        if i ~= j
        
        %Pi(indices(:,i),indices(:,j))'
        %Pi(indices(:,i),indices(:,j))'*Prob(indices(:,i))
        
        %assuming tau=1
        aux = 0;
        for ii = 1:size(indices(:,i),1);
            aux = aux + sum(Pi(indices(ii,i),indices(:,j))*Prob(indices(ii,i)));
        end
        Joint_Prob(i,j) = aux;
        
        end
        
    end
end
%Joint_Prob





%%

MU = @(u) u.^(-gamma); %marginal utility function
U = MU(dividend); %define marginal utility vector curly U
D = diag(U);
delta_star = delta - gamma*mu+.5*gamma*(gamma+1)*sigma^2;
%A = delta*Pi.*(MU(dividend)*(1./MU(dividend))')'; %A-D matrix

y_t = NaN(T,2^n);
B_t = NaN(T,2^n);
e = ones(2^n,1);
for t = 1:T
    B_t(t,:) = exp(-t*(1-delta_star))*mpower(inv(D)*Pi*D,t)*e;
    y_t(t,:) = -log(B_t(t,:))/t+1;
    %y_t(t,:) = (1./(sum(mpower(A,t),2))).^(1/t);
    %y_t(t,:) = (1./(sum(mpower(A,t),2))-1)/t+1;
    %y_t(t,:) = (-log(sum(mpower(A,t),2)))/t+1;    
end




%%

%line_style = {'-','--',':','-.','.-','.-.','.'};
line_style = {'-','--',':','-.'};
line_color = {'k','b','g','m'}; %,'g','y','m','c'};
Alphabet=char('a'+(1:26)-1)';

%legend record
leg_record = cell(n-1,max(unique_state_x_dist));


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
for k=2:size(unique_states,1)-1
       
    h(state_num_dist(k)) =...
        subplot(fig_dim(2),fig_dim(1),state_num_dist(k));
    p(k) = plot(y_t(:,unique_states(k))-1);
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
    choose_line = 1+floor((cyc_counter-1)/cyc_1_size);
    p(k).Color = line_color{choose_color};
    p(k).LineStyle = line_style{choose_line};

    
    title(sprintf('%s) %d distressed',...
        Alphabet(state_num_dist(k,1)),state_num_dist(k,1)))
    xlabel(sprintf('Maturity (t)'))
    ylabel(sprintf('Yield (R_f)'))

%Legend construction
%number of distressed firms
     leg{cyc_counter} = str{state_num_dist(k,2)};

    %When reaching the end of cycle 2 restart the counter
    if cyc_counter == unique_state_x_dist(1+state_num_dist(k,1))
        hold on
        plot(y_t(:,1)-1,'r');
        leg{cyc_counter+1} = str{1};
        hold on
        plot(y_t(:,2^n)-1,'r--');
        leg{cyc_counter+2} = str{2^n};
        cyc_counter = 0;
%Legend is differemt for each n. Create legend for n=4
        legend(leg,...
            'Location','southeast','Orientation','vertical');
        legend('boxoff');
        leg_record(state_num_dist(k),1:size(leg,2)) = leg;
        axis([0,T,min(y_t(1,1:end-1)-1),max(y_t(1,1:end-1)-1)])
        
    end
    
    
end

%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(ring_fig,'-dpng','-r80',...
    sprintf('ring_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}

%%

width_fig = 1024;
height_fig = 1024;
x_pos = 100;
y_pos = 50;

% Choose figure to dispaly
subfig = 4;
c=copyobj(h(subfig),get(h(subfig),'Parent'));
ring_subfig = figure;
set(c,'Parent',ring_subfig);
set(ring_subfig, 'Position', [x_pos y_pos width_fig height_fig])
ha=gca;
set(c,'Position','default');
hay=get(ha,'Ylabel');
hax=get(ha,'Xlabel');
NewFontSize=10;
set(hay,'Fontsize',NewFontSize);
set(hax,'Fontsize',NewFontSize);

graph_title = strcat(sprintf('Ring network, n=%d (%d distressed)',n,subfig),...
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
leg{end+1} = str{1};
leg{end+1} = str{2^n};
legend(leg,...
            'Location','southeast','Orientation','vertical');
legend('boxoff');


%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(ring_subfig,'-dpng','-r80',...
sprintf('ring_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
%}

%%

%{
Create a figure of concentrated shocks and disperssed shocks
For 7!
%}

%find concentrated and unique states
states_concentrated = NaN(n,2);
states_dispersed = NaN(floor(n/2)+1,2);
states_dispersed(1,:) = [1, 0];
kk = 0;
%set the second column as the number of distressed
for k = 1:size(states_dispersed,1)-1    
    kk = kk+2^(2*(k-1));
    states_dispersed(k+1,:) = [kk+1, k-1];
end
for k = 1:n
    states_concentrated(k,:) = [2^(k-1), k-1];
end



leg = cell(size(states_concentrated,1),1);
%leg{1} = dec2bin(0,n);

fig = figure;
set(gcf, 'Visible', 'off')

cyc_counter = 0;
for k=1:size(states_concentrated,1)
    cyc_counter = cyc_counter + 1;
    p(k) = plot(y_t(:,states_concentrated(k,1))-1);
    hold on
    
    %legend
    leg{k} = dec2bin(states_concentrated(k,1)-1,n);

    
    cyc_1_size = size(line_color,2);
    cyc_2_size = size(line_style,2);
   
    choose_color = cyc_counter-cyc_1_size*...
        floor((cyc_counter-1)/cyc_1_size);
    choose_line = 1+floor((cyc_counter-1)/cyc_1_size);
    p(k).Color = line_color{choose_color};
    p(k).LineStyle = line_style{choose_line};
    
end


graph_title = strcat('Concentrated shocks ( ',...
    '\lambda=',sprintf('%.2f; ',lambda),...
    '\eta=',sprintf('%.2f; ',eta),...
    'n=',sprintf('%d)',n));
title(graph_title)
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))

legend(leg,...
    'Location','southeast','Orientation','vertical');
legend('boxoff');







%{
cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(ring_fig,'-dpng','-r80',...
    sprintf('ring_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}






