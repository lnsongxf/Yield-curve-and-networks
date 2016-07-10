clear; clc

%define parameters
n = 8; % # nodes
load(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\chain_%d',n))



T = 50; % # periods forward running the risk free rate
delta = 0.05; %time prefrence
gamma = 1; %RRA parameter of CRRA
mu = 0.015; %growth rate in the economy
sigma = 0.1; %volatility of the dividend

%Dividend at each state
epsilon = .2;
x(2) = 1; %x_bar, value of idiosyncratic component at boom state
x(1) = x(2)-epsilon; %x_lowerbar, value of idiosyncratic component at busat state

Delta = chain_net(n); %Ring network

dividend = x(2)*(n-state_num_dist(:,1))+x(1)*state_num_dist(:,1);
Pi = NaN(size(state_num_dist,1));

%spectral gap record parameters
lambda_start_pt = .01;
lambda_end_pt = .49;
eta_start_pt = .01;
eta_end_pt = .49;
steps = 9;
LAMBDA = linspace(lambda_start_pt,lambda_end_pt,steps);
ETA = linspace(eta_start_pt,eta_end_pt,steps);

spectral_gap_record = NaN(steps);
for lambda_ind = 1:steps
    for eta_ind = 1:steps
        
    clc
    sprintf('lambda: %d/%d, eta: %d/%d',lambda_ind,steps,eta_ind,steps)
        
        lambda= (LAMBDA(lambda_ind)); %internsity of distress
        eta= (ETA(eta_ind)); %internsity of distress
        
        risk_ratio = eta/lambda;

%{
Constructing the transition matrix
The way to do so is to start at a unique state -->
then compute probabilities of transition from that unique state to all
other states - but! do it in groups of equiv states to keep track of 
the probability from one unique state to the other -->
at the end sum up the appropriate probabilities and we have an entry in Pi
%}
for k_origin = 1:size(state_num_dist,1)


str_ = dec2bin(state_num_dist(k_origin,2)-1,n); %generate binary number k-1
digit = str_-'0'; %arrange each digit in a row vector
H_origin = logical(digit'); %define this vector as state of the world

Lambda = lambda*Delta*H_origin;
Eta = eta*ones(n,1);
%Eta = eta*Delta*(1-H_origin); %network effects for recovery

%Natural rate of distress/recovery
%Lambda = lambda*Delta*H(:,k)+lambda;
%Eta = eta*Delta*(1-H(:,k))+eta;

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

        eigval = sort(real(eig(Pi)),'descend');
        spectral_gap_record(lambda_ind,eta_ind) = eigval(1)-eigval(2);

    end
end
%%


width_fig = 1024+325;
height_fig = 512+70;
x_pos = 10;
y_pos = 100;

%set(gcf, 'Visible', 'off')
fig_spectral_gap = figure;
subfig(1) = subplot(1,2,1);
contourf(ETA,LAMBDA,spectral_gap_record)
colorbar

graph_title = strcat(sprintf('Spectral gap, ring with %d nodes',n));
title(graph_title)
xlabel('\eta')
ylabel('\lambda','rot',0)

xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 -.015 0])
set(xlabh,'FontSize', 14)

ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(xlabh,'Position') - [.289 -.243 0])
set(ylabh,'FontSize', 14)

subfig(2) = subplot(1,2,2);
line_style = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
begin_pt = steps-5;
leg = cell(steps-begin_pt+1,1);

for i = begin_pt:steps
    choose_line_style = i-begin_pt+1;
    plot(LAMBDA/ETA(i),spectral_gap_record(:,i),line_style{choose_line_style})
    hold on
    
    leg{i-begin_pt+1} = strcat('\eta=',sprintf('%.2f',ETA(i)));
end

graph_title = strcat(sprintf('Spectral gap, ring with %d nodes, varrying',n),' \eta/\lambda');

title(graph_title)
%title(strcat(sprintf('Spectral gap, ring with %d nodes, varrying',n),' \eta/\lambda'))



%title('ddd')
xlabel('\lambda/\eta')
ylabel(sprintf('Spectral gap'))

xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [.2 -.012 0])
set(xlabh,'FontSize', 14)

ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [-.2 .02 0])
set(ylabh,'FontSize', 14)

legend(leg,...
            'Location','northeast','Orientation','vertical');
legend('boxoff');

set(fig_spectral_gap, 'Position', [x_pos y_pos width_fig height_fig])
set(subfig(1),'position',[.08 .1 .3375 .73]) % sets figure size
set(subfig(2),'position',[.58 .1 .35 .73]) % sets figure size

%{
%small eta
upto = 1;
plot(LAMBDA/ETA(upto),spectral_gap_record(:,upto))
hold on
%medium eta
upto = floor(steps/2);
plot(LAMBDA/ETA(upto),spectral_gap_record(:,upto))
hold on
%large eta
upto = steps;
plot(LAMBDA/ETA(upto),spectral_gap_record(:,upto))
%}

%{
print(fig_spectral_gap,'-dpng','-r100',...
    sprintf('spectral_gap_n=%d_rr=%.2f.jpg',...
    n,risk_ratio))
%}






