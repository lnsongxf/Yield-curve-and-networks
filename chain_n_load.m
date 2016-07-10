clear; clc

clear; clc

%define parameters

n = 8; % # nodes
load(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\chain_%d',n))

num_state = size(state_num_dist,1);

T = 50; % # periods forward running the risk free rate
risk_ratio = 2;
m = .2;
lambda = m; %internsity of distress
eta = m*risk_ratio; %probability of recovery
delta = .05; %time prefrence
gamma = 1; %RRA parameter of CRRA
mu = 0.015; %growth rate in the economy
sigma = 0.1; %volatility of the dividend

%Dividend at each state
epsilon = .2;
x(2) = 1; %x_bar, value of idiosyncratic component at boom state
x(1) = x(2)-epsilon; %x_lowerbar, value of idiosyncratic component at busat state


Delta = chain_net(n);

%%


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

%%
%{
*************************************************************************
                        State probability distribution
*************************************************************************
%}

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
A = exp(-delta_star)*(D\eye(num_state)*Pi*D);

y_t = NaN(T,num_state);
B_t = NaN(T,num_state);
e = ones(num_state,1);
R = NaN(T,T,num_state);
for t = 1:T
    B_t(t,:) = mpower(A,t)*e;
    y_t(t,:) = -log(B_t(t,:))/t;
    for tau = 1:t-1
        R(t,tau,:) = log(((mpower(Pi,tau)*(mpower(A,t-tau))*e))./(mpower(A,t)*e))+log(mpower(A,tau)*e);%tau*y_t(tau,:)';
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

tau = 6;
t = T;
highlight = 19;
dist_num = 1;
%choose_states = sort(unique([3, 5, 7, highlight]));
[choose_states, aux] = find(state_num_dist(:,1)==dist_num);
choose_states = choose_states';



graph_title = strcat(sprintf(...
    '                                         Chain network, n=%d',n),...
    sprintf('\n Parametrization:'),...
    ' \lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon)...
    );

figure('Position', [10 50 1350 600]);

h1 = subplot(1,3,1);
shade1 = [min(y_t,[],2), max(y_t,[],2)-min(y_t,[],2)];
a1 = area(shade1,'FaceColor','k','edgealpha',0.1);
alpha .2
a1(1).FaceColor = [ 1 1 1];
hold on
YLim(1) = min(min(y_t(:,choose_states)));
YLim(2) =  max(max(y_t(:,choose_states)));
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])

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

t2 = title(strcat('b) Bond premium,',sprintf('\nkeeping holding period at'),' \tau=',...
    sprintf('%d',tau)));
xlabel('T (underlying bond maturity)','FontSize',fs1)
YLim(1) = min(min(R(:,tau,choose_states)));
YLim(2) = max(max(R(:,tau,choose_states)));
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

t3 = title(strcat('c) Bond premium,',sprintf('\nkeeping underlying bond maturity period at T='),...
    sprintf('%d',t)));

xlabel('\tau (holding period)','FontSize',fs1)
YLim(1) = min(min(R(t,:,choose_states)));
YLim(2) = max(max(R(t,:,choose_states)));
%get(h2,'YLim')
axis([0, 50, YLim(1)-abs(YLim(1))*border, YLim(2)+abs(YLim(2))*border])

fig_height = .55; push_up = .05;
set(h1,'Position',get(h1,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h2,'Position',get(h2,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])
set(h3,'Position',get(h3,'Position').*[1 2 1 fig_height]+[0 push_up 0 0])



set(t1,'Position',get(t1,'Position').*[1 1 1]+[0 1e-3 0])
set(t2,'Position',get(t2,'Position').*[1 0 1]+[1 4e-4 0])
get(t2,'Position')
set(t3,'Position',get(t3,'Position').*[1 0 1]+[0 4.05e-4 0])
%get(t3,'Position')


set(t1,'FontSize',fs2)
set(t2,'FontSize',fs2)
set(t3,'FontSize',fs2)

text_pos = get(xl,'Position').*[1 1 0]+[0 -.001 0];
text_pos = text_pos(1:2);
text(-97,-0.00012, graph_title,'FontSize',16);

return

%%
%{
**************************************************************************
                                FIGURES
YIELD CURVES: 1 FIGURE WITH SUBFIGURES ACCORDING TO NUMBER OF DISTRESSED
**************************************************************************
%}

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
%%
cyc_counter = 0;
chain_fig = figure('Visible', 'off');
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
print(chain_fig,'-dpng','-r80',...
    sprintf('chain_n%d_eta%.2f_lambda%.2f_gamma%.1f.jpg',...
    n,eta,lambda,gamma))
%}

%%

%{
**************************************************************************
                                FIGURES
                             With pictures
**************************************************************************
%}

clc

height_fig = 512+75;
width_fig = 1.12*height_fig;
x_pos = 100;
y_pos = 75;



% Choose figure to dispaly
subfig = 1;

graph_title = strcat(sprintf('Chain network, n=%d (%d distressed);',n,subfig),...
    sprintf(' Parametrization:\n --'),...
    '\lambda=',sprintf('%.2f ',lambda),...
    ' \eta=',sprintf('%.2f',eta),...
    ' \gamma=',sprintf('%.2f',gamma),...
    ' \mu=',sprintf('%.3f',mu),...
    ' \sigma=',sprintf('%.2f',sigma),...
    ' \delta^*=',sprintf('%.2f',delta_star),...
    ' \epsilon=',sprintf('%.2f',epsilon),' --');

if unique_state_x_dist(subfig+1) < 14


    c=copyobj(h(subfig),get(h(subfig),'Parent'));

    chain_subfig = figure('Visible', 'on');
    set(chain_subfig, 'Position', [x_pos y_pos width_fig height_fig]);


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

    set(c,'Parent',chain_subfig);
    set(c,'Position',newhandle.*[1,1,aux,aux]);
    axis off

    ha=gca;
    hay=get(ha,'Ylabel');
    hax=get(ha,'Xlabel');
    NewFontSize=10;
    set(hay,'Fontsize',NewFontSize);
    set(hax,'Fontsize',NewFontSize);
    axis(c,[0 50 0.044 0.05])

    %subplot(1,3,1)

    title(c,graph_title);

    leg = leg_record(subfig,1:unique_state_x_dist(subfig+1));
    leg{end+1} = dec2bin(0,n);
    leg{end+1} = dec2bin(2^n-1,n);
    legend(c,leg,...
                'Location','northeast','Orientation','vertical');
    legend(c,'boxoff');

    s_number = state_num_dist(state_num_dist(:,1)==subfig,2);
    addition = repmat([-.37,-.05,0,0], [unique_state_x_dist(subfig+1),1]);
    x_mult = 2; y_mult = 1;
    multiplication = repmat([1.45,1,x_mult,y_mult], [unique_state_x_dist(subfig+1),1]);

    pic_position = [];
    aux_ = min(unique_state_x_dist(subfig+1),aux);
    for i = 1:aux_
        pic_position = [pic_position, aux*i];
    end

    if unique_state_x_dist(subfig+1) > aux
        for i = 1:unique_state_x_dist(subfig+1)-aux
            pic_position = [pic_position, aux*(aux-1)+i];
            
        end
    end
if n == 7   
    for i = aux+1:size(addition,1)
        addition(i,:) = addition(i,:)+[.3 0 0 0];
    end

    for i = 1:aux
        addition(i,:) = addition(i,:)+[-.03 .125 0 0];
    end
end    
    
    for i = 1:unique_state_x_dist(subfig+1)
        newhandle = subplot(aux,aux,pic_position(i));
        set(newhandle,'Position',get(newhandle,'Position').*multiplication(i,:)+addition(i,:),...
            'Visible', 'off');
        pic = chain_pic(n,s_number(i));
        c=copyobj(pic,get(h(subfig),'Parent'));
        set(c,'Parent',chain_subfig,'Position',get(newhandle,'Position'));
        title(c,leg{i});
    end


    
    cd C:\Users\Oren\Documents\MATLAB\Network\figures
    print(chain_subfig,'-dpng','-r80',...
    sprintf('chain_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
    n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
    



    else

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
    sprintf('chain_n%d_%ddist_lambda0%.0f_eta0%.0f_gamma%.0f_mu00%.0f_sigma0%.0f_delta0%.0f_epsilon0%.0f.png',...
    n,subfig,lambda*100,eta*100,gamma,mu*1000,sigma*10,delta*100,epsilon*10))
    %}


    %}

end




