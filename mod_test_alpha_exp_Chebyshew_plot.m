clc, clf

t_min=0;t_max=20; m=500; x_initial = [5*10^6; 10^3; 10^3];%För att ha startvärden nära ett steady state sätt in [3.1298490038*10^9; 10^-5; 402531911.894]; 
time_mesh =[];
for k = 1:m
    time_mesh(end+1) = -(t_max-t_min)/2*cos((2*k-1)*pi/(2*m)) + (t_max+t_min)/2;
end


alpha =  [10^-11 10^-12 10^-10 10^-12 10^-12]*100;  %[1,1,1,1,1];%
alpha_1=alpha_vec(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),time_mesh);
x_23s = ForwardODE23s(alpha_1,time_mesh,x_initial); 
x_Newton = ForwardNewton(alpha_1,time_mesh,x_initial);  
x_45 = ForwardODE45(alpha_1,time_mesh,x_initial);
t_plot = [1:500];


gron = [102,194,165]/255;
orange = [252,141,98]/255;
lila = [141,160,203]/255;

alpha_unknown=4;
alpha_exp_23s=mod_calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton=mod_calculate_alpha_exp(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45=mod_calculate_alpha_exp(alpha,alpha_unknown,x_45,t_min,t_max);

%orginal
alpha_exp_23s_org=calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton_org=calculate_alpha_exp(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45_org=calculate_alpha_exp(alpha,alpha_unknown,x_45,t_min,t_max);

max_45_org = max(alpha_exp_45_org);
max_Newton_org = max(alpha_exp_Newton_org);
max_23s_org = max(alpha_exp_23s_org);

max_45_mod = max(alpha_exp_45);
max_Newton_mod = max(alpha_exp_Newton);
max_23s_mod = max(alpha_exp_23s);


max_org = max([max_45_org max_Newton_org max_23s_org]);
max_mod = max([max_45_mod max_Newton_mod max_23s_mod]);

min_45_org = min(alpha_exp_45_org);
min_Newton_org = min(alpha_exp_Newton_org);
min_23s_org = min(alpha_exp_23s_org);

min_45_mod= min(alpha_exp_45);
min_Newton_mod = min(alpha_exp_Newton);
min_23s_mod = min(alpha_exp_23s);

min_org = min([min_45_org min_Newton_org min_23s_org]);
min_mod = min([min_45_mod min_Newton_mod min_23s_mod]);

figure('name','Chebyshew')
subplot(2,1,1)
plot(time_mesh(2:end-1),alpha_exp_23s,'color', gron ,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','True value of parameter')
title('Modifierad')
ylim([min_mod-max_mod/3,max_mod*1.3])


subplot(2,1,2)
plot(time_mesh(2:end-1),alpha_exp_23s_org,'color',gron,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','True value of parameter')
title('Orginal')
ylim([min_org-max_org/3,max_org*1.3])
hold on
%%
figure



gron = [102,194,165]/255;
orange = [252,141,98]/255;
lila = [141,160,203]/255;

pT = plot(time_mesh,x_45(1,:), 'Color',gron);
hold on
pM1 = plot( time_mesh,x_45(2,:),'color', orange);
hold on
pM2 = plot(time_mesh,x_45(3,:),  'color',lila );
hold on

pT.LineWidth = 1;
pM1.LineWidth = 1;
pM2.LineWidth = 1;

legend('x_T', 'x_{M1}', 'x_{M2}');
xlabel('Mätpunkt');
ylabel('Densitet')



%% 
function alpha = alpha_vec(dm1,dm2,at1,at2,k12,time_mesh)
scaling_factor_dm1 = dm1;
scaling_factor_dm2 = dm2;
scaling_factor_at1 = at1;
scaling_factor_at2 = at2;
scaling_factor_k12 = k12;

function_flag = 0; % constant

exact_dm1 = ExactParameter(scaling_factor_dm1,function_flag,time_mesh); %Exact profile for dm1 to produce data.
exact_dm2 = ExactParameter(scaling_factor_dm2,function_flag,time_mesh); %Exact profile for dm2 to produce data.
exact_at1 = ExactParameter(scaling_factor_at1,function_flag,time_mesh); %Exact profile for at1 to produce data.
exact_at2 = ExactParameter(scaling_factor_at2,function_flag,time_mesh); %Exact profile for at2 to produce data.
exact_k12 = ExactParameter(scaling_factor_k12,function_flag,time_mesh); %Exact profile for k12 to produce data.

alpha = [exact_dm1; exact_dm2; exact_at1; exact_at2; exact_k12];

end