clear; close all; clc;

%% Initial Conditions vector

V_E = 0; 
V_0 = 0; % No virus input
V_I = 0;
R_cyt = 0;
R_CM = 0;
P_S = 0;
P_NS = 0;
RC_CM = 0;
R_ds = 0;
RIGI = 5.3422;
aRIGI = 0;
MAVS = 277.78;
aMAVS = 0;
IKKe = 3.0832;
pIKKe = 0;
TBK1 = 97.169;
pTBK1 = 0;
IRF3 = 37.862;
pIRF3 = 0;
IKK = 37.969;
aIKK = 0;
NFkB_IkBac = 11.3587;
pNFkBn = 0;
NFkBn = 0;
NFkBc = 101.735;
IkBac = 0;
IRF7 = 24.92;
pIRF7 = 0;
IFNb_mRNA = 0;
IFNa_mRNA = 0;
IFNl_mRNA = 0;
IFN_c = 0;
IFNl_c = 0;
JAK = 151.88;
RJC = 0;
STAT1c = 1114.8;
CP = 20;
ISGn = 0;
IFNex = 0; % No IFN input
STAT2c = 6.5019;
TYK = 20.701;
RTC = 0;
ARC = 0;
IFNAR1 = 1000;
IFNAR2 = 1000;
IFNAR_d = 0;
IRF9c = 45;
ARC_STAT2c = 0;
ARC_STAT12c = 0;
STAT2c_IRF9 = 0;
ISGF3c = 0;
PSC_c = 0;
ISGF3_CP = 0;
PSC_CP = 0;
NP = 40;
STAT1n = 0;
STAT2n = 0;
PIAS = 41.96;
PSC_n = 0;
IRF9n = 0;
ISGF3n = 0;
PSC_NP = 0;
B_u = 500;
B_o_NP = 0;
B_o = 0;
ISGF3_PIAS = 0;
STAT2n_IRF9 = 0;
ISGF3_NP = 0;
ISGav_mRNA = 0;
ISGav = 0;
ISGn_mRNA_n = 0;
IRF9_mRNA_n = 0;
IRF7_mRNA = 0;
ISGn_mRNA_c = 0;
IRF9_mRNA_c = 0;

Init_Cond = [V_E V_0 V_I R_cyt R_CM P_S P_NS RC_CM R_ds RIGI aRIGI MAVS aMAVS IKKe pIKKe TBK1 ...
    pTBK1 IRF3 pIRF3 IKK aIKK NFkB_IkBac pNFkBn NFkBn NFkBc IkBac IRF7 pIRF7 IFNb_mRNA IFNa_mRNA IFNl_mRNA IFN_c IFNl_c JAK RJC STAT1c CP ISGn IFNex STAT2c TYK RTC ARC IFNAR1 IFNAR2 IFNAR_d IRF9c ARC_STAT2c ARC_STAT12c STAT2c_IRF9 ISGF3c PSC_c ISGF3_CP PSC_CP NP STAT1n STAT2n PIAS PSC_n IRF9n ISGF3n PSC_NP B_u B_o_NP B_o ISGF3_PIAS STAT2n_IRF9 ISGF3_NP ISGav_mRNA ISGav ISGn_mRNA_n IRF9_mRNA_n IRF7_mRNA ISGn_mRNA_c IRF9_mRNA_c];
%%
load('param.mat')

I_n = 1; I_a = 0; VC = 0;

% Simulate for 120 hours
t_start_ss = 0;
t_end_ss = 120 * 60; % 120 hours × 60 min/hour
tspan_ss = linspace(t_start_ss, t_end_ss);

% Solve ODE
[Tss, Yss] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan_ss, Init_Cond);

% Save results
save('SteadyState_120h.mat', 'Tss', 'Yss');

%%
clear; close all; clc;

load('SteadyState_120h.mat', 'Tss', 'Yss');
load('param.mat')

y0 = Yss(end,:);
y0(2) = 100;           % MOI =10;  

tend = 48*60;
tspan = linspace(0,tend);

c = [0, 1];

var_names = {'ExtVirus', 'VirusInit', 'IntVirus', 'R_{cyt}', '(+)RNA_{CM}', 'SP', 'NSP', 'RC_CM', 'dsRNA', 'RIGI','aRIGI','MAVS','aMAVS',...
    'IKKe','pIKKe','TBK1','pTBK1', 'IRF3','pIRF3','IKK','aIKK','NFkBIkBac','pNFkBn','NFkBn','NFkBc', 'IkBac', 'IRF7', 'pIRF7', 'IFNbmRNA',...
    'IFNamRNA','IFNlmRNA', 'IFN_c', 'IFNl_c', 'JAK','RJC', 'STAT1c','CP', 'ISGn','IFNex','STAT2c','TYK','RTC','ARC', 'Rec1','Rec2',...
    'IFNARd','IRF9_c','ARC-STAT2_c', 'ARC-STAT12_c','STAT2-IRF9_c','ISGF3_c', 'PSC_c','ISGF3-CP','PSC-CP','NP','STAT1_n','STAT2_n','PIAS','PSC_n',...
    'IRF9_n','ISGF3_n','PSC-NP','B_u','B_o-NP','B_o','ISGF3-PIAS','STAT2-IRF9_n','ISGF3n-NP', 'ISGavmRNA','ISGav', 'ISGnmRNA_n', 'IRF9mRNA_n',...
    'IRF7mRNA', 'ISGnmRNA_c', 'IRF9mRNA_c'};

tic
for i = 1:size(c,2)

   I_n = c(i);

    for j = 1: size(c,2)

        VC = c(j); 
        
        for k = 1:size(c,2)

            I_a = c(k);

            [T,Y] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan,y0);
            disp(1)
            save(['HCV_48h_VC', num2str(VC), '_Ia', num2str(I_a),'_In',num2str(I_n), '.mat'], 'T', 'Y');

            T = T/60;

            for m = [1, 70] 

                f = figure(m);
                % f.WindowState = 'maximized';
                set(f,'units','points','position',[0,0,600,400])
                plot(T, (Y(:,m)),'LineWidth',1.5,'DisplayName',  ['isgt(' num2str(I_a) ')VC(' num2str(VC) ')isgn(' num2str(I_n) ')']);
                hold on
                ylabel(var_names{1, m},'FontSize',15, 'Interpreter', 'tex')
                xlabel('time [h]','FontSize',15,'Interpreter','tex')
                set(gca,'FontSize',15, 'Yscale', 'log')
                set(gca,'YMinorTick','off')
                set(gca,'LineWidth',1.5)
                set(gca, 'Color', 'none')
                set(gca,'TickLabelInterpreter','Latex')
                xlim([4 48])
                xticks([4 8 16 24 36 48])
                legend('Location','bestoutside')
                % saveas(gcf, var_names{1,m}, 'fig')
            end
        end 
    end 
end 
toc