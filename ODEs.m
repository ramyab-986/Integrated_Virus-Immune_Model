function deriv = ODEs(t,y,Param, I_n, I_a, VC)

% Define Parameters
Vc2n             = Param(1);
Vn2c             = Param(2);
b_IRF3           = Param(3);
b_MAVS           = Param(4);
b_RIGI           = Param(5);
b_prot           = Param(6);
k1               = Param(7);
k10              = Param(8);
k11              = Param(9);
k12              = Param(10);
k13              = Param(11);
k14              = Param(12);
k15              = Param(13);
k16              = Param(14);
k17              = Param(15);
k18              = Param(16);
k19              = Param(17);
k2               = Param(18);
k20              = Param(19);
k21              = Param(20);
k22              = Param(21);
k23              = Param(22);
k27              = Param(23);
k29              = Param(24);
k3               = Param(25);
k31              = Param(26);
k32              = Param(27);
k34              = Param(28);
k35              = Param(29);
k36              = Param(30);
k37              = Param(31);
k38              = Param(32);
k39              = Param(33);
k4               = Param(34);
k40              = Param(35);
k41              = Param(36);
k42              = Param(37);
k43              = Param(38);
k44              = Param(39);
k45              = Param(40);
k46              = Param(41);
k47              = Param(42);
k48              = Param(43);
k49              = Param(44);
k5               = Param(45);
k50              = Param(46);
k51              = Param(47);
k52              = Param(48);
k53              = Param(49);
k54              = Param(50);
k56              = Param(51);
k57              = Param(52);
k58              = Param(53);
k59              = Param(54);
k6               = Param(55);
k60              = Param(56);
k61              = Param(57);
k62              = Param(58);
k63              = Param(59);
k64              = Param(60);
k65              = Param(61);
k66              = Param(62);
k67              = Param(63);
k7               = Param(64);
k8               = Param(65);
k9               = Param(66);
k_IFN            = Param(67);
k_IKK            = Param(68);
k_IKKe_TBK1      = Param(69);
k_IRF3_IKKe_TBK1 = Param(70);
k_MAVS           = Param(71);
k_RIGI           = Param(72);
k_TFBS_IFNa      = Param(73);
k_TFBS_IFNb      = Param(74);
k_TFBS_IFNl      = Param(75);
k_act            = Param(76);
k_deph           = Param(77);
k_expr_IkBa      = Param(78);
k_inh_p65        = Param(79);
k_mRNA_IFNb      = Param(80);
k_mRNA_IFNl      = Param(81);
k_mRNA_MX1       = Param(82);
k_mx1_mRNA       = Param(83);
k_rigi_synt      = Param(84);
k_transISG       = Param(85);
k_trans_IFNl     = Param(86);
k_transp_NFkB    = Param(87);
mu_IFN           = Param(88);
mu_IFNl          = Param(89);
mu_IkBa          = Param(90);
mu_mRNA_IFNa     = Param(91);
mu_mRNA_IFNb     = Param(92);
mu_mRNA_IFNl     = Param(93);
mu_ISG_RNA      = Param(94);
mu_mx1           = Param(95);
mu_rigi          = Param(96);
k_en             = Param(97);
k_f              = Param(98);
mu_V_I           = Param(99);
k_a              = Param(100);
k_e              = Param(101);
k_r              = Param(102);
k_c              = Param(103);
k_t              = Param(104);
tau              = Param(105);
a_RC             = Param(106);
rcsat            = Param(107);
mu_r             = Param(108);
mu_p             = Param(109);
mu_V_E           = Param(110);
nSP              = Param(111);
ks               = Param(112);
k71             = Param(113);
k79             = Param(114);
k76              = Param(115);
deg              = Param(116);
k72             = Param(117);
k74             = Param(118);
k73             = Param(119);
k75             = Param(120);
tau_5            = Param(121);
k78             = Param(122);
k77             = Param(123);
k70            = Param(124);
k_transISGn             = Param(125);
k69            = Param(126); % 5.05e-04 units = 1/(molecules.min)
degRecBySOCS     = Param(127);
degARCBySOCS     = Param(128);
kinhBySOCS       = Param(129);


% Define species 

V_E = y(1);
V_0 = y(2);
V_I = y(3);
R_cyt = y(4);
R_CM = y(5);
P_S = y(6);
P_NS = y(7);
RC_CM = y(8);
R_ds = y(9);
RIGI = y(10);
aRIGI = y(11);
MAVS = y(12);
aMAVS = y(13);
IKKe = y(14);
pIKKe = y(15);
TBK1 = y(16);
pTBK1 = y(17);
IRF3 = y(18);
pIRF3 = y(19);
IKK = y(20);
aIKK = y(21);
NFkB_IkBac = y(22);
pNFkBn = y(23);
NFkBn = y(24);
NFkBc = y(25);
IkBac = y(26);
IRF7 = y(27);
pIRF7 = y(28);
IFNb_mRNA = y(29);
IFNa_mRNA = y(30);
IFNl_mRNA = y(31);
IFN_c = y(32);
IFNl_c = y(33);
JAK = y(34);
RJC = y(35);
STAT1c = y(36);
CP = y(37);
ISGn = y(38);
IFNex = y(39);
STAT2c = y(40);
TYK = y(41);
RTC = y(42);
ARC = y(43);
IFNAR1 = y(44);
IFNAR2 = y(45);
IFNAR_d = y(46);
IRF9c = y(47);
ARC_STAT2c = y(48);
ARC_STAT12c = y(49);
STAT2c_IRF9 = y(50);
ISGF3c = y(51);
PSC_c = y(52);
ISGF3_CP = y(53);
PSC_CP = y(54);
NP = y(55);
STAT1n = y(56);
STAT2n = y(57);
PIAS = y(58);
PSC_n = y(59);
IRF9n = y(60);
ISGF3n = y(61);
PSC_NP = y(62);
B_u = y(63);
B_o_NP = y(64);
B_o = y(65);
ISGF3_PIAS = y(66);
STAT2n_IRF9 = y(67);
ISGF3_NP = y(68);
ISGav_mRNA = y(69);
ISGav = y(70);
ISGn_mRNA_n = y(71);
IRF9_mRNA_n = y(72);
IRF7_mRNA = y(73);
ISGn_mRNA_c = y(74);
IRF9_mRNA_c = y(75);

% System of ODEs

% Virallifecycle                                                        
   
f_CM   = 1 - exp(-(t/tau)^4);
 
RC_form_rate = (k_c) * R_cyt * P_NS *( (f_CM) - ( RC_CM/rcsat ));
    
dV_E    =  (k_a/(1 + I_a*ISGav))*P_S*R_cyt - mu_V_E*V_E;

dV_0    =  (-k_en/(1 + I_a*ISGav))*V_0;

dV_I    =  (k_en/(1 + I_a*ISGav))*V_0 - (k_f/(1 + I_a*ISGav))*V_I - mu_V_I*V_I;

dR_cyt    =  (k_e)*R_CM - (k_a/(1 + I_a*ISGav))*P_S*R_cyt - mu_r*(1 + (I_a*ISGav))*R_cyt - RC_form_rate + (k_f/(1 + I_a*ISGav))*V_I;

dR_CM    =  (k_r/(1 + I_a*ISGav))*RC_CM  - (k_e)*R_CM;

dP_S    =  (k_t/(1 + I_a*ISGav))*R_cyt - (k_a/(1 + I_a*ISGav))*nSP*P_S*R_cyt - mu_p*(1 + (I_a*ISGav))*P_S;

dP_NS    =  (k_t/(1 + I_a*ISGav))*R_cyt - RC_form_rate - mu_p*(1 + (I_a*ISGav))*P_NS;

dRC_CM    =  RC_form_rate - a_RC*RC_CM;

dR_ds    =  a_RC*RC_CM + b_RIGI*aRIGI - mu_r*(1 + (I_a*ISGav))*R_ds - k_RIGI*RIGI*R_ds;


% NF-kB pathway & IRFs Pathway 
  
dRIGI 	 = k_rigi_synt - RIGI*mu_rigi + aRIGI*b_RIGI - RIGI*R_ds*k_RIGI + ISGav_mRNA*k_transISG;  %keep the basal expression at 5% for KO simulations of RIGI

daRIGI 	 = RIGI*R_ds*k_RIGI - aRIGI*mu_rigi - aRIGI*b_RIGI ;

dMAVS 	 = aMAVS*b_MAVS - MAVS*aRIGI*(k_MAVS/(1 + (VC*P_NS)));
 
daMAVS 	 = MAVS*aRIGI*(k_MAVS/(1 + (VC*P_NS))) - aMAVS*b_MAVS ;

dIKKe 	 = b_prot*pIKKe - IKKe*aMAVS*k_IKKe_TBK1 ;

dpIKKe 	 = IKKe*aMAVS*k_IKKe_TBK1 - b_prot*pIKKe ;

dTBK1 	 = b_prot*pTBK1 - TBK1*aMAVS*k_IKKe_TBK1 ;

dpTBK1 	 = TBK1*aMAVS*k_IKKe_TBK1 - b_prot*pTBK1 ;

dIRF3 	 = Vn2c*b_IRF3*pIRF3 - (IRF3*k_IRF3_IKKe_TBK1*(pIKKe + pTBK1))*(1/(1 + VC*P_NS)) ;

dpIRF3 	 = (IRF3*Vc2n*k_IRF3_IKKe_TBK1*(pIKKe + pTBK1))*(1/(1 + VC*P_NS)) - b_IRF3*pIRF3 ;

dIKK 	 = aIKK*b_prot - IKK*aMAVS*k_IKK;

daIKK 	 = IKK*aMAVS*k_IKK - aIKK*b_prot ;

dNFkB_IkBac 	 = IkBac*NFkBc*k_inh_p65 - NFkB_IkBac*aIKK*k_act ;

dpNFkBn 	 = Vc2n*aIKK*k_act*(NFkBc + NFkB_IkBac) - k_deph*pNFkBn ;

dNFkBn 	 = k_deph*pNFkBn - NFkBn*k_transp_NFkB ;

dNFkBc 	 = NFkBn*Vn2c*k_transp_NFkB - IkBac*NFkBc*k_inh_p65 - NFkBc*aIKK*k_act ;

dIkBac 	 = k_expr_IkBa*(NFkBn + pNFkBn) - IkBac*mu_IkBa - IkBac*NFkBc*k_inh_p65 ;

dIRF7   = k71*Vn2c*pIRF7 - IRF7*k_IRF3_IKKe_TBK1*(pIKKe + pTBK1)*(1/(1 + VC*P_NS)) + k79*IRF7_mRNA - deg*IRF7;

dpIRF7   = IRF7*Vc2n*k_IRF3_IKKe_TBK1*(pIKKe + pTBK1)*(1/(1 + VC*P_NS)) - k71*pIRF7;

dIFNb_mRNA 	 = B_o*Vn2c*k_TFBS_IFNb - IFNb_mRNA*mu_mRNA_IFNb + (Vn2c*k_mRNA_IFNb*(pIRF3 + pIRF7)*(NFkBn + pNFkBn)) ;

dIFNa_mRNA 	 = B_o*Vn2c*k_TFBS_IFNa - IFNa_mRNA*mu_mRNA_IFNa + (Vn2c*k_mRNA_IFNb*(pIRF3 + pIRF7)) ;

dIFNl_mRNA 	 = B_o*Vn2c*k_TFBS_IFNl - IFNl_mRNA*mu_mRNA_IFNl + Vn2c*k_mRNA_IFNl*pIRF3 ;

dIFN_c = (IFNa_mRNA + IFNb_mRNA)*k_IFN - ks*IFN_c - mu_IFN*IFN_c;

dIFNl_c 	 = IFNl_mRNA*k_trans_IFNl - IFNl_c*mu_IFNl ;


% JAK|STAT Pathway

dJAK 	 = ARC*k32 + RJC*k4 + ARC*I_n*ISGn*k69 - JAK*IFNAR2*k3 ;

dRJC 	 = IFNAR_d*k6 - RJC*k4 + JAK*IFNAR2*k3 - IFNex*RJC*RTC*k5 ;

dSTAT1c 	 = ARC_STAT12c*k12 + ISGF3_CP*k40 - STAT1c*k56 + PSC_CP*k43 - ARC_STAT2c*STAT1c*k11 + STAT1n*Vn2c*k57 ;

dCP 	 = ISGF3_CP*k39 + ISGF3_CP*k40 + PSC_CP*k42 + PSC_CP*k43 - CP*ISGF3c*k38 - CP*PSC_c*k41 ;

dISGn 	 = k_transISGn*ISGn_mRNA_c - ISGn*k31 ;

dIFNex 	 = IFNAR_d*k6 + ks*IFN_c - IFNex*RJC*RTC*k5 ;

dSTAT2c 	 = ARC_STAT2c*k10 + ISGF3_CP*k40 - STAT2c*k58 + STAT2c_IRF9*k36 + STAT2c_IRF9*k61 + PSC_CP*k43 - ARC*STAT2c*k9 - IRF9c*STAT2c*k60 + STAT2n*Vn2c*k59 ;

dTYK 	 = ARC*k32 + RTC*k2 + ARC*I_n*ISGn*k69 - IFNAR1*TYK*k1 ;

dRTC 	 = IFNAR_d*k6 - RTC*k2 + IFNAR1*TYK*k1 - IFNex*RJC*RTC*k5 ;

dARC 	 = ARC_STAT2c*k10 - ARC*k34 - ARC*k32 + ARC_STAT12c*(k13)*(1/(1 + VC*P_NS)) + IFNAR_d*(k7/(1+I_n*kinhBySOCS*ISGn)) - ARC*I_n*ISGn*k69 - ARC*STAT2c*k9 - ARC*STAT2c_IRF9*k8 - I_n*degARCBySOCS*ISGn*ARC ;

dIFNAR1 	 = ARC*k32 + RTC*k2 + ARC*I_n*ISGn*k69 - IFNAR1*TYK*k1 - I_n*degRecBySOCS*IFNAR1*ISGn ;

dIFNAR2 	 = ARC*k32 + RJC*k4 + ARC*I_n*ISGn*k69 - JAK*IFNAR2*k3 - I_n*degRecBySOCS*IFNAR2*ISGn ;

dIFNAR_d 	 = ARC*k34 - IFNAR_d*k6 - IFNAR_d*(k7/(1+I_n*kinhBySOCS*ISGn)) + IFNex*RJC*RTC*k5 ;

dIRF9c 	 = k27 - IRF9c*k29 - IRF9c*k66 + ISGF3c*k15 + ISGF3_CP*k40 + STAT2c_IRF9*k61 + k70*IRF9_mRNA_c + ARC*STAT2c_IRF9*k8 - IRF9c*STAT2c*k60 - IRF9c*PSC_c*k14 + IRF9n*Vn2c*k67 ;

dARC_STAT2c 	 = ARC_STAT12c*k12 - ARC_STAT2c*k10 + ARC*STAT2c*k9 + ARC*STAT2c_IRF9*k8 - ARC_STAT2c*STAT1c*k11 ;

dARC_STAT12c 	 = ARC_STAT2c*STAT1c*k11 - ARC_STAT12c*k13*(1/(1 + VC*P_NS)) - ARC_STAT12c*k12 ;

dSTAT2c_IRF9 	 = IRF9c*STAT2c*k60 - STAT2c_IRF9*k61 - STAT2c_IRF9*k64 - ARC*STAT2c_IRF9*k8 - STAT2c_IRF9*k36 + STAT2n_IRF9*Vn2c*k65 ;

dISGF3c 	 = ISGF3_CP*k39 - ISGF3c*k16 - ISGF3c*k15 - CP*ISGF3c*k38 + IRF9c*PSC_c*k14 + ISGF3n*Vn2c*k17 ;

dPSC_c 	 = ARC_STAT12c*k13*(1/(1 + VC*P_NS)) + ISGF3c*k15 - PSC_c*k18 + PSC_CP*k42 - CP*PSC_c*k41 - IRF9c*PSC_c*k14 + PSC_n*Vn2c*k19 ;

dISGF3_CP 	 = CP*ISGF3c*k38 - ISGF3_CP*k40 - ISGF3_CP*k39 ;

dPSC_CP 	 = CP*PSC_c*k41 - PSC_CP*k43 - PSC_CP*k42 ;

dNP 	 = ISGF3_NP*k48 + ISGF3_NP*k49 + B_o_NP*k51 + B_o_NP*k52 + PSC_NP*k45 + PSC_NP*k46 - ISGF3n*NP*k47 - NP*B_o*k50 - NP*PSC_n*k44 ;

dSTAT1n 	 = ISGF3_NP*k49 + B_o_NP*k52 - STAT1n*k57 + PSC_NP*k46 + STAT1c*Vc2n*k56 ;

dSTAT2n 	 = ISGF3_NP*k49 + B_o_NP*k52 - STAT2n*k59 + STAT2n_IRF9*k37 + STAT2n_IRF9*k63 + PSC_NP*k46 - IRF9n*STAT2n*k62 + STAT2c*Vc2n*k58 ;

dPIAS 	 = ISGF3_PIAS*k54 - ISGF3n*PIAS*k53 ;

dPSC_n 	 = ISGF3n*k21 - PSC_n*k19 + PSC_NP*k45 - IRF9n*PSC_n*k20 - NP*PSC_n*k44 + PSC_c*Vc2n*k18 ;

dIRF9n 	 = ISGF3n*k21 - IRF9n*k67 - IRF9n*k35 + ISGF3_NP*k49 + B_o_NP*k52 + STAT2n_IRF9*k63 - IRF9n*STAT2n*k62 - IRF9n*PSC_n*k20 + IRF9c*Vc2n*k66 ;

dISGF3n 	 = ISGF3_NP*k48 - ISGF3n*k21 - ISGF3n*k17 + B_o*k23 + ISGF3_PIAS*k54 - ISGF3n*NP*k47 - ISGF3n*B_u*k22 - ISGF3n*PIAS*k53 + IRF9n*PSC_n*k20 + ISGF3c*Vc2n*k16 ;

dPSC_NP	 = NP*PSC_n*k44 - PSC_NP*k46 - PSC_NP*k45 ;

dB_u	 = B_o*k23 + B_o_NP*k52 - ISGF3n*B_u*k22 ;

dB_o_NP  = NP*B_o*k50 - B_o_NP*k52 - B_o_NP*k51 ;

dB_o 	 = B_o_NP*k51 - B_o*k23 + ISGF3n*B_u*k22 - NP*B_o*k50 ;

dISGF3_PIAS 	 = ISGF3n*PIAS*k53 - ISGF3_PIAS*k54 ;

dSTAT2n_IRF9 	 = IRF9n*STAT2n*k62 - STAT2n_IRF9*k63 - STAT2n_IRF9*k65 - STAT2n_IRF9*k37 + STAT2c_IRF9*Vc2n*k64 ;

dISGF3_NP 	 = ISGF3n*NP*k47 - ISGF3_NP*k49 - ISGF3_NP*k48 ;

dISGav_mRNA 	 = B_o*Vn2c*k_mx1_mRNA - ISGav_mRNA*mu_ISG_RNA + Vn2c*k_mRNA_MX1*pIRF3 ;

dISGav 	 = ISGav_mRNA*k_transISG - ISGav*mu_mx1 ;

dISGn_mRNA_n   = k72*B_o - k73*ISGn_mRNA_n;

dIRF9_mRNA_n   = k74*B_o - k75*IRF9_mRNA_n;

dIRF7_mRNA   = k76*B_o - (log(2)/tau_5)*IRF7_mRNA;

dISGn_mRNA_c  = Vn2c*ISGn_mRNA_n*k73 - k77*ISGn_mRNA_c;

dIRF9_mRNA_c  = Vn2c*IRF9_mRNA_n*k75 - k78*IRF9_mRNA_c;

% Derivatives 

deriv = [dV_E; dV_0; dV_I; dR_cyt; dR_CM; dP_S; dP_NS; dRC_CM; dR_ds; dRIGI; daRIGI; dMAVS; daMAVS; dIKKe; dpIKKe; dTBK1; dpTBK1; dIRF3; dpIRF3; dIKK; daIKK; dNFkB_IkBac; dpNFkBn; dNFkBn; dNFkBc; dIkBac; dIRF7; dpIRF7; dIFNb_mRNA; dIFNa_mRNA; dIFNl_mRNA; dIFN_c; dIFNl_c; dJAK; dRJC; dSTAT1c; dCP; dISGn; ...
         dIFNex; dSTAT2c; dTYK; dRTC; dARC; dIFNAR1; dIFNAR2; dIFNAR_d; dIRF9c; dARC_STAT2c; dARC_STAT12c; dSTAT2c_IRF9; dISGF3c; dPSC_c; dISGF3_CP; dPSC_CP; dNP; dSTAT1n; dSTAT2n; dPIAS; dPSC_n; dIRF9n; dISGF3n; dPSC_NP; dB_u; dB_o_NP; dB_o; dISGF3_PIAS; dSTAT2n_IRF9; dISGF3_NP; dISGav_mRNA; dISGav; dISGn_mRNA_n; dIRF9_mRNA_n; dIRF7_mRNA; dISGn_mRNA_c; dIRF9_mRNA_c];

end