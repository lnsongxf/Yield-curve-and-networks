clear; clc

%define parameters

n = 8; % # nodes
load(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\star_%d',n))
num_state = size(state_num_dist,1);

T = 50; % # periods forward running the risk free rate
risk_ratio = 1;
m = .1;
lambda = m; %internsity of distress
eta = m*risk_ratio; %probability of recovery
delta = .05; %time prefrence
gamma = 4; %RRA parameter of CRRA
mu = 0.015; %growth rate in the economy
sigma = 0.1; %volatility of the dividend

%Dividend at each state
epsilon = .5;
x(2) = 1; %x_bar, value of idiosyncratic component at boom state
x(1) = x(2)-epsilon; %x_lowerbar, value of idiosyncratic component at busat state


%Directionality of the network
%0 non directional star network, two way eff
%1     directional star network, center contaminates but cannot be contaminated
directional = 0;
Delta = star_net(n,directional); %firm number 1 is the center

if directional == 0 && lambda > 1/n
    disp('Error: In a non-directional star network largest lambda is 1/n')
    return
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
H_origin = logical(digit)'; %define this vector as state of the world

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
indices = indices(~isnan(indices));
aux_ = NaN(size(indices,1),1);

% then run over each of those states and find the prob to get to that state
% from k_origin, and then sum these probs up.
for ind = 1:size(indices,1)
  
str_ = dec2bin(indices(ind)-1,n);
digit = str_-'0';
H_dest = logical(digit)';

aux = NaN(n,1);
for i = 1:n
aux(i) = (1-Eta(i))*H_origin(i)*H_dest(i)+Eta(i)*H_origin(i)*(1-H_dest(i))+...
    Lambda(i)*(1-H_origin(i))*H_dest(i)+(1-Lambda(i))*(1-H_origin(i))*(1-H_dest(i));
end
         aux_(ind) = prod(aux);
end
         
Pi(k_origin,k_dest) = sum(aux_);
        
end
end


%{
a = mpower(Pi,1e10);
plot(a(1,:),'-bo')
hold on
plot(sum(a,2),'-r')

axis([1,size(state_num_dist,1),0,1.1])
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
A = exp(-delta_star)*(D\eye(size(state_num_dist,1))*Pi*D);

y_t = NaN(T,size(state_num_dist,1));
B_t = NaN(T,size(state_num_dist,1));
e = ones(size(state_num_dist,1),1);
R = NaN(T,T,size(state_num_dist,1));
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

%{
**************************************************************************
                                FIGURES
                     BOND PREMIA AND YIELD CURVES
**************************************************************************
%}

tau = 10;
t = T;
highlight = 2*n-1;
%dist_num = 1;
choose_states = sort(unique([2, 3, 2*n-1, 2*n-2, highlight]));
%[choose_states, aux] = find(state_num_dist(:,1)==dist_num);
%choose_states = choose_states';

if directional == 1
    graph_title = strcat(sprintf('                                        Star network, n=%d',n),...
        sprintf('\n            directional: center contaminates but is not contaminated'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.2f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
elseif directional == 0
    graph_title = strcat(sprintf('                                        Star network, n=%d',n),...
        sprintf('\n             non-directional: contagion works on center as well'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.4f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
    
end






figure('Position', [10 50 1350 600]);

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

%shade2 = [min(R(:,tau,:),[],3), max(R(:,tau,:),[],3)];
shade2 = [min(R(:,tau,:),[],3), max(R(:,tau,:),[],3)-min(R(:,tau,:),[],3)];
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

t2 = title(strcat('b) Bond premium, keeping holding period \tau=',...
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

t3 = title(strcat('c) Bond premium, keeping holding period T=',...
    sprintf('%d',t)));

xlabel('\tau (holding period)','FontSize',fs1)
YLim(1) = min(min(R(t,:,:)));
YLim(2) = max(max(R(t,:,:)));
%get(h2,'YLim')
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])

fig_height = .55; push_up = .05;
set(h1,'Position',get(h1,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h2,'Position',get(h2,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h3,'Position',get(h3,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])



set(t1,'Position',get(t1,'Position').*[1 1 1]+[0 1.6e-3 0])
set(t2,'Position',get(t2,'Position').*[1 1 1]+[-3 1e-7 0])
set(t3,'Position',get(t3,'Position').*[1 1 1]+[0 1e-7 0])

set(t1,'FontSize',fs2)
set(t2,'FontSize',fs2)
set(t3,'FontSize',fs2)

text_pos = get(xl,'Position').*[1 1 0]+[0 -.001 0];
text_pos = text_pos(1:2);
text(-95,-2.9e-4, graph_title,'FontSize',16);


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

if directional == 1
    graph_title = strcat(sprintf('Star network, n=%d',n),...
        sprintf('\ndirectional: center contaminates but is not contaminated'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.2f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
elseif directional == 0
    graph_title = strcat(sprintf('Star network, n=%d',n),...
        sprintf('\nnon-directional: contagion works on center as well'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.4f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
    
end


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
pic = star_pic(n,s_number);
c=copyobj(pic,get(h1,'Parent'));
set(c,'Position',pic_pos);




%%

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
star_fig = figure;%('Position', [0, 0, 1365, 690]);
set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]);

set(gcf, 'Visible', 'off')

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
print(star_fig,'-dpng','-r80',...
    sprintf('star_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}


%%

width_fig = 512+100;
height_fig = 512;
x_pos = 100;
y_pos = 50;

% Choose figure to dispaly
subfig = 1;
c=copyobj(h(subfig),get(h(subfig),'Parent'));
star_subfig = figure('Visible','on');
set(c,'Parent',star_subfig);
set(star_subfig, 'Position', [x_pos y_pos width_fig height_fig])
ha=gca;
set(c,'Position','default');
hay=get(ha,'Ylabel');
hax=get(ha,'Xlabel');
NewFontSize=10;
set(hay,'Fontsize',NewFontSize);
set(hax,'Fontsize',NewFontSize);

graph_title = strcat(sprintf('Star network, n=%d (%d distressed)',n,subfig),...
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
print(star_subfig,'-dpng','-r80',...
sprintf('star_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
%}

%%

%center distressed in full line -, and if not --
%number if distressed in different color. 0 distressed in green, all
%distressed in red, the rest in various colors.

%Colors = {'g','b','c','k','r'};
star_fig = figure('Visible','on');
leg = cell(1,6);
%Center not distressed, i-1 periphery distressed

% there are four colors, so grouping n states to 4 groups we have x in each
% group: x = ceil((n-2)/4)
state_num_in_group = ceil((n-2)/4);
counter = 1;
for i = [1,2,n]
    if i == 1
        color_used = 'g-';
    elseif i > 1 && i <= state_num_in_group+1
        color_used = 'c--';
    elseif i > 1*state_num_in_group+1 && i <= 2*state_num_in_group+1
        color_used = 'b--';
    elseif i > 2*state_num_in_group+1 && i <= 3*state_num_in_group+1
        color_used = 'k--';
    elseif i > 3*state_num_in_group+1 && i <= 4*state_num_in_group+1
        color_used = 'm--';
    end
    
    leg{counter} = strcat('0',num2str(i-1));
    plot(y_t(:,state_num_dist(:,2) == 2^(i-1)),color_used) % i-1 periphery firm is distressed; center not
    hold on
    counter = counter + 1;
end


%Center distressed, i-1 periphery distressed
for i = [1,n-1,n] %[1,2,n,n+1]
    if i <= 1*state_num_in_group
        color_used = 'c-';    
    elseif i > 1*state_num_in_group && i <= 2*state_num_in_group
        color_used = 'b-';
    elseif i > 2*state_num_in_group && i <= 3*state_num_in_group
        color_used = 'k-';
    elseif i > 3*state_num_in_group && i <= 4*state_num_in_group && i ~= n
        color_used = 'm-';
    elseif i == n
        color_used = 'r-';    
    end
    
    leg{counter} = strcat('1',num2str(i-1));    
    plot(y_t(:,state_num_dist(:,2) == 2^(n-1)+2^(i-1)),color_used) % 0 periphery firm is distressed; center yes
    hold on
    counter = counter + 1;    
end


if directional == 1
    graph_title = strcat(sprintf('Star network, n=%d',n),...
        sprintf('\ndirectional: center contaminates but is not contaminated'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.2f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
elseif directional == 0
    graph_title = strcat(sprintf('Star network, n=%d',n),...
        sprintf('\nnon-directional: contagion works on center as well'),...
        sprintf('\nParametrization:'),...
        ' \lambda=',sprintf('%.2f ',lambda),...
        ' \eta=',sprintf('%.4f',eta),...
        ' \gamma=',sprintf('%.2f',gamma),...
        ' \mu=',sprintf('%.3f',mu),...
        ' \sigma=',sprintf('%.2f',sigma),...
        ' \delta^*=',sprintf('%.2f',delta_star),...
        ' \epsilon=',sprintf('%.2f',epsilon)...
        );
    
end


title(graph_title)
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))
%axis([0,50,0,0.13]) %for comparison!



legend(leg,'Location','northeast','Orientation','vertical')

cd C:\Users\Oren\Documents\MATLAB\Network\figures
if directional == 1
    print(star_fig,'-dpng','-r80',...
    sprintf('star_directional_n%d_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
    n,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
elseif directional == 0
    print(star_fig,'-dpng','-r80',...
    sprintf('star_non-directional_n%d_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
    n,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
end




