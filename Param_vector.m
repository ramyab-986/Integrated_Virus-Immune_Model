%% Import data from spreadsheet


%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 129);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:DY2";

% Specify column names and types
opts.VariableNames = ["Vc2n", "Vn2c", "a_RC", "b_IRF3", "b_MAVS", "b_RIGI", "b_prot", "deg", "degARCBySOCS", "degRecBySOCS", "k1", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k2", "k20", "k21", "k22", "k23", "k27", "k29", "k3", "k31", "k32", "k34", "k35", "k36", "k37", "k38", "k39", "k4", "k40", "k41", "k42", "k43", "k44", "k45", "k46", "k47", "k48", "k49", "k5", "k50", "k51", "k52", "k53", "k54", "k56", "k57", "k58", "k59", "k6", "k60", "k61", "k62", "k63", "k64", "k65", "k66", "k67", "k69", "k7", "k70", "k71", "k72", "k73", "k74", "k75", "k76", "k77", "k78", "k79", "k8", "k9", "k_IFN", "k_IKK", "k_IKKe_TBK1", "k_IRF3_IKKe_TBK1", "k_MAVS", "k_RIGI", "k_TFBS_IFNa", "k_TFBS_IFNb", "k_TFBS_IFNl", "k_a", "k_act", "k_c", "k_deph", "k_e", "k_en", "k_expr_IkBa", "k_f", "k_inh_p65", "k_mRNA_IFNb", "k_mRNA_IFNl", "k_mRNA_MX1", "k_mx1_mRNA", "k_r", "k_rigi_synt", "k_t", "k_transISG", "k_transISGn", "k_trans_IFNl", "k_transp_NFkB", "kinhBySOCS", "ks", "mu_IFN", "mu_IFNl", "mu_IkBa", "mu_V_E", "mu_V_I", "mu_mRNA_IFNa", "mu_mRNA_IFNb", "mu_mRNA_IFNl", "mu_mRNA_mx1", "mu_mx1", "mu_p", "mu_r", "mu_rigi", "nSP", "rcsat", "tau", "tau_5"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
tbl = readtable("Param.xlsx", opts, "UseExcel", false);

%% Convert to output type
Vc2n = tbl.Vc2n;
Vn2c = tbl.Vn2c;
a_RC = tbl.a_RC;
b_IRF3 = tbl.b_IRF3;
b_MAVS = tbl.b_MAVS;
b_RIGI = tbl.b_RIGI;
b_prot = tbl.b_prot;
deg = tbl.deg;
degARCBySOCS = tbl.degARCBySOCS;
degRecBySOCS = tbl.degRecBySOCS;
k1 = tbl.k1;
k10 = tbl.k10;
k11 = tbl.k11;
k12 = tbl.k12;
k13 = tbl.k13;
k14 = tbl.k14;
k15 = tbl.k15;
k16 = tbl.k16;
k17 = tbl.k17;
k18 = tbl.k18;
k19 = tbl.k19;
k2 = tbl.k2;
k20 = tbl.k20;
k21 = tbl.k21;
k22 = tbl.k22;
k23 = tbl.k23;
k27 = tbl.k27;
k29 = tbl.k29;
k3 = tbl.k3;
k31 = tbl.k31;
k32 = tbl.k32;
k34 = tbl.k34;
k35 = tbl.k35;
k36 = tbl.k36;
k37 = tbl.k37;
k38 = tbl.k38;
k39 = tbl.k39;
k4 = tbl.k4;
k40 = tbl.k40;
k41 = tbl.k41;
k42 = tbl.k42;
k43 = tbl.k43;
k44 = tbl.k44;
k45 = tbl.k45;
k46 = tbl.k46;
k47 = tbl.k47;
k48 = tbl.k48;
k49 = tbl.k49;
k5 = tbl.k5;
k50 = tbl.k50;
k51 = tbl.k51;
k52 = tbl.k52;
k53 = tbl.k53;
k54 = tbl.k54;
k56 = tbl.k56;
k57 = tbl.k57;
k58 = tbl.k58;
k59 = tbl.k59;
k6 = tbl.k6;
k60 = tbl.k60;
k61 = tbl.k61;
k62 = tbl.k62;
k63 = tbl.k63;
k64 = tbl.k64;
k65 = tbl.k65;
k66 = tbl.k66;
k67 = tbl.k67;
k69 = tbl.k69;
k7 = tbl.k7;
k70 = tbl.k70;
k71 = tbl.k71;
k72 = tbl.k72;
k73 = tbl.k73;
k74 = tbl.k74;
k75 = tbl.k75;
k76 = tbl.k76;
k77 = tbl.k77;
k78 = tbl.k78;
k79 = tbl.k79;
k8 = tbl.k8;
k9 = tbl.k9;
k_IFN = tbl.k_IFN;
k_IKK = tbl.k_IKK;
k_IKKe_TBK1 = tbl.k_IKKe_TBK1;
k_IRF3_IKKe_TBK1 = tbl.k_IRF3_IKKe_TBK1;
k_MAVS = tbl.k_MAVS;
k_RIGI = tbl.k_RIGI;
k_TFBS_IFNa = tbl.k_TFBS_IFNa;
k_TFBS_IFNb = tbl.k_TFBS_IFNb;
k_TFBS_IFNl = tbl.k_TFBS_IFNl;
k_a = tbl.k_a;
k_act = tbl.k_act;
k_c = tbl.k_c;
k_deph = tbl.k_deph;
k_e = tbl.k_e;
k_en = tbl.k_en;
k_expr_IkBa = tbl.k_expr_IkBa;
k_f = tbl.k_f;
k_inh_p65 = tbl.k_inh_p65;
k_mRNA_IFNb = tbl.k_mRNA_IFNb;
k_mRNA_IFNl = tbl.k_mRNA_IFNl;
k_mRNA_MX1 = tbl.k_mRNA_MX1;
k_mx1_mRNA = tbl.k_mx1_mRNA;
k_r = tbl.k_r;
k_rigi_synt = tbl.k_rigi_synt;
k_t = tbl.k_t;
k_transISG = tbl.k_transISG;
k_transISGn = tbl.k_transISGn;
k_trans_IFNl = tbl.k_trans_IFNl;
k_transp_NFkB = tbl.k_transp_NFkB;
kinhBySOCS = tbl.kinhBySOCS;
ks = tbl.ks;
mu_IFN = tbl.mu_IFN;
mu_IFNl = tbl.mu_IFNl;
mu_IkBa = tbl.mu_IkBa;
mu_V_E = tbl.mu_V_E;
mu_V_I = tbl.mu_V_I;
mu_mRNA_IFNa = tbl.mu_mRNA_IFNa;
mu_mRNA_IFNb = tbl.mu_mRNA_IFNb;
mu_mRNA_IFNl = tbl.mu_mRNA_IFNl;
mu_mRNA_mx1 = tbl.mu_mRNA_mx1;
mu_mx1 = tbl.mu_mx1;
mu_p = tbl.mu_p;
mu_r = tbl.mu_r;
mu_rigi = tbl.mu_rigi;
nSP = tbl.nSP;
rcsat = tbl.rcsat;
tau = tbl.tau;
tau_5 = tbl.tau_5;

%% Clear temporary variables
clear opts tbl

%%
% Define parameter order
    param = [
        Vc2n Vn2c b_IRF3 b_MAVS b_RIGI...
        b_prot k1 k10 k11 k12 k13 k14 k15 k16 k17 k18...
        k19 k2 k20 k21 k22 k23 k27 k29 k3 k31 k32 k34...
        k35 k36 k37 k38 k39 k4 k40 k41 k42 k43 k44 k45 k46 k47 k48...
        k49 k5 k50 k51 k52 k53 k54 k56 k57 k58...
        k59 k6 k60 k61 k62 k63 k64 k65 k66 k67...
        k7 k8 k9 k_IFN k_IKK k_IKKe_TBK1 k_IRF3_IKKe_TBK1 ...
        k_MAVS k_RIGI k_TFBS_IFNa k_TFBS_IFNb...
        k_TFBS_IFNl k_act k_deph k_expr_IkBa k_inh_p65... 
        k_mRNA_IFNb k_mRNA_IFNl ...
        k_mRNA_MX1 k_mx1_mRNA k_rigi_synt k_transISG k_trans_IFNl ...
        k_transp_NFkB mu_IFN mu_IFNl mu_IkBa mu_mRNA_IFNa ...
        mu_mRNA_IFNb mu_mRNA_IFNl mu_mRNA_mx1 mu_mx1 mu_rigi ...
        k_en k_f mu_V_I k_a k_e k_r k_c k_t tau a_RC ...
        rcsat mu_r mu_p mu_V_E nSP ks k71 k79 k76 ...
        deg k72 k74 k73 k75 tau_5 k78 k77 k70... 
        k_transISGn k69 degRecBySOCS degARCBySOCS kinhBySOCS
    ];

save('param.mat', 'param')