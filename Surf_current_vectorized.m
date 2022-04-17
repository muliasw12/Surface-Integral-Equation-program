clearvars -except tStart...
    am_1_1...
    bm__1...
    outer_pp_1...
    outer_pm_1...
    outer_mp_1...
    outer_mm_1...
    tStart

load parameter_segitiga
load singular_arcioni_tp
load singular_integral

x1 = eV.';
y1_re = n.';
y1_im = k.';
eV_interp = (0.1:0.0001:29.92).';
% eV_interp = (1.5:0.001:6).';
n_re = interp1(x1,y1_re,eV_interp,'pchip');
n_im = interp1(x1,y1_im,eV_interp,'pchip');

%% Material and background parameters (clear)

row         = 40334;
h_bar_eV    = 6.582119e-16;
omega       = eV_interp(row)./h_bar_eV;
lamda       = 1240*10^-9/eV_interp(row);
silver      = n_re(row).^2-n_im(row).^2 + 1i*2.*n_re(row).*n_im(row);
vacuum      = 1;
water       = 1.77;
eps_r       = 11.68;
eps_bg      = 1;

eps_0       = 8.854e-12;
mu_0        = 1.257e-6;
c           = 1/sqrt(eps_0*mu_0);
k_vec       = [0 0 1];
E_0         = [1 0 0];

k_bg        = 2*pi.*sqrt(eps_bg)./lamda;
kv_bg       = k_bg*k_vec;
k           = 2*pi.*sqrt(eps_r)./lamda;
kv          = k*k_vec;

imp_bg      = sqrt(mu_0/(eps_0*eps_bg));
imp         = sqrt(mu_0/(eps_0*eps_r));

%% Green scalar & smooth (clear)

obs1 = titik_tengah(:,1);
obs2 = titik_tengah(:,2);
obs3 = titik_tengah(:,3);
source1 = titik_tengah(:,1).';
source2 = titik_tengah(:,2).';
source3 = titik_tengah(:,3).';

c1 = bsxfun(@minus, obs1, source1);
c2 = bsxfun(@minus, obs2, source2);
c3 = bsxfun(@minus, obs3, source3);
R = sqrt(c1.^2 + c2.^2 + c3.^2);
logicR = R==0;
logicR_ = reshape(repmat(R==0,1,3),jumSegitiga,jumSegitiga,3);

greenskalar = exp(1i*k*R)./(4*pi*R);
greenskalar_bg = exp(1i*k_bg*R)./(4*pi*R);
greensmooth = ((exp(1i*k.*R)-1)./R + k.^2.*R./2)./(4*pi);
greensmooth_bg = ((exp(1i*k_bg.*R)-1)./R + k_bg.^2.*R./2)./(4*pi);

grad_greenskalar_x = (1i).*k.*c1.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c1.*exp((1i).*k.*R)./(4.*pi.*R.^3);
grad_greenskalar_y = (1i).*k.*c2.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c2.*exp((1i).*k.*R)./(4.*pi.*R.^3);
grad_greenskalar_z = (1i).*k.*c3.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c3.*exp((1i).*k.*R)./(4.*pi.*R.^3);
gradientgreen = reshape([grad_greenskalar_x grad_greenskalar_y grad_greenskalar_z],jumSegitiga,jumSegitiga,3);

grad_greenskalar_x = (1i).*k_bg.*c1.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c1.*exp((1i).*k_bg.*R)./(4.*pi.*R.^3);
grad_greenskalar_y = (1i).*k_bg.*c2.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c2.*exp((1i).*k_bg.*R)./(4.*pi.*R.^3);
grad_greenskalar_z = (1i).*k_bg.*c3.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c3.*exp((1i).*k_bg.*R)./(4.*pi.*R.^3);
gradientgreen_bg = reshape([grad_greenskalar_x grad_greenskalar_y grad_greenskalar_z],jumSegitiga,jumSegitiga,3);

grad_greensmooth_x = k.^2.*c1./(8.*pi.*R)+(1i).*k.*c1.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c1.*(exp((1i).*k.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth_y = k.^2.*c2./(8.*pi.*R)+(1i).*k.*c2.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c2.*(exp((1i).*k.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth_z = k.^2.*c3./(8.*pi.*R)+(1i).*k.*c3.*exp((1i).*k.*R)./(4.*pi.*R.^2)...
    -c3.*(exp((1i).*k.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth = reshape([grad_greensmooth_x grad_greensmooth_y grad_greensmooth_z],jumSegitiga,jumSegitiga,3);

grad_greensmooth_x = k_bg.^2.*c1./(8.*pi.*R)+(1i).*k_bg.*c1.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c1.*(exp((1i).*k_bg.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth_y = k_bg.^2.*c2./(8.*pi.*R)+(1i).*k_bg.*c2.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c2.*(exp((1i).*k_bg.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth_z = k_bg.^2.*c3./(8.*pi.*R)+(1i).*k_bg.*c3.*exp((1i).*k_bg.*R)./(4.*pi.*R.^2)...
    -c3.*(exp((1i).*k_bg.*R)-1)./(4.*pi.*R.^3);
grad_greensmooth_bg = reshape([grad_greensmooth_x grad_greensmooth_y grad_greensmooth_z],jumSegitiga,jumSegitiga,3);

temp = ones(jumSegitiga,jumSegitiga).*(1i)*k/(4*pi);
greensmooth(logicR) = temp(logicR);
temp = ones(jumSegitiga,jumSegitiga).*(1i)*k_bg/(4*pi);
greensmooth_bg(logicR) = temp(logicR);
temp = zeros(jumSegitiga,jumSegitiga,3);
grad_greensmooth(logicR_) = temp(logicR_);
temp = zeros(jumSegitiga,jumSegitiga,3);
grad_greensmooth_bg(logicR_) = temp(logicR_);

gradientgreen_source = -gradientgreen;
gradientgreen_source_bg = -gradientgreen_bg;
gradientgreensmooth_source = -grad_greensmooth;
gradientgreensmooth_source_bg = -grad_greensmooth_bg;

%% Incoming wave matrix (clear)

f_plus = repmat(PanjangSisi,1,3).*rho_plus./(2*repmat(luas(SegitigaPlus),1,3));
f_minus = -repmat(PanjangSisi,1,3).*rho_minus./(2*repmat(luas(SegitigaMinus),1,3));
div_f_plus = PanjangSisi./luas(SegitigaPlus);
div_f_minus = -PanjangSisi./luas(SegitigaMinus);

EmPlus = repmat(E_0,jumSisi,1).*repmat(exp((1i)*dot(titik_tengah(SegitigaPlus,:),repmat(kv_bg,jumSisi,1),2)),1,3);
EmMinus = repmat(E_0,jumSisi,1).*repmat(exp((1i)*dot(titik_tengah(SegitigaMinus,:),repmat(kv_bg,jumSisi,1),2)),1,3);
q1 = dot(f_plus,EmPlus,2).*luas(SegitigaPlus)+dot(f_minus,EmMinus,2).*luas(SegitigaMinus);
HmPlus = repmat(cross(k_vec,E_0),jumSisi,1).*repmat(exp((1i)*dot(titik_tengah(SegitigaPlus,:),repmat(kv_bg,jumSisi,1),2)),1,3)./(mu_0*c);
HmMinus = repmat(cross(k_vec,E_0),jumSisi,1).*repmat(exp((1i)*dot(titik_tengah(SegitigaMinus,:),repmat(kv_bg,jumSisi,1),2)),1,3)./(mu_0*c);
q2 = dot(f_plus,HmPlus,2).*luas(SegitigaPlus)+dot(f_minus,HmMinus,2).*luas(SegitigaMinus);

q = [q1 ; q2].';

%% Logical edge matrix (clear)

logic_identical_pp = logical(logic_identical(SegitigaPlus,SegitigaPlus));
logic_identical_pm = logical(logic_identical(SegitigaPlus,SegitigaMinus)); %AWAS KEBALIK PM MP
logic_identical_mp = logical(logic_identical(SegitigaMinus,SegitigaPlus));
logic_identical_mm = logical(logic_identical(SegitigaMinus,SegitigaMinus));

logic_adjacent_pp = logical(logic_adjacent(SegitigaPlus,SegitigaPlus));
logic_adjacent_pm = logical(logic_adjacent(SegitigaPlus,SegitigaMinus)); %AWAS KEBALIK PM MP
logic_adjacent_mp = logical(logic_adjacent(SegitigaMinus,SegitigaPlus));
logic_adjacent_mm = logical(logic_adjacent(SegitigaMinus,SegitigaMinus));

logic_touch_pp = logical(logic_touch(SegitigaPlus,SegitigaPlus));
logic_touch_pm = logical(logic_touch(SegitigaPlus,SegitigaMinus)); %AWAS KEBALIK PM MP
logic_touch_mp = logical(logic_touch(SegitigaMinus,SegitigaPlus));
logic_touch_mm = logical(logic_touch(SegitigaMinus,SegitigaMinus));

logic_far_pp = logical(logic_far(SegitigaPlus,SegitigaPlus));
logic_far_pm = logical(logic_far(SegitigaPlus,SegitigaMinus)); %AWAS KEBALIK PM MP
logic_far_mp = logical(logic_far(SegitigaMinus,SegitigaPlus));
logic_far_mm = logical(logic_far(SegitigaMinus,SegitigaMinus));

%% + + D_matrix region 1

sum_1 = K_2_1_plus(SegitigaPlus,:,:)./(4*pi)-k_bg.^2.*K_2_2_plus(SegitigaPlus,:,:)./(8*pi)+...
    repmat(greensmooth_bg(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_plus(SegitigaPlus,:)./(4*pi)-k_bg.^2.*K_1_2_plus(SegitigaPlus,:)./(8*pi)+...
    greensmooth_bg(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
outer_pp = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;

sum_2_sing = K_1_1_plus(SegitigaPlus,:)./(4*pi)-k_bg.^2.*K_1_2_plus(SegitigaPlus,:)./(8*pi);
sum_2_smooth = greensmooth_bg(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
outer_pp_sing = sum_2_sing.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;
outer_pp_smooth = sum_2_smooth.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;

int_smooth_div = greensmooth_bg(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
int_smooth = repmat(greensmooth_bg(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;
Cmn = bsxfun(@times,div_f_plus,div_f_plus.')./4;
int_arcioni_div = bsxfun(@times,div_f_plus,div_f_plus.').*I_1_pp./(4*pi);
int_arcioni = Cmn.*I_2_pp./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k_bg.^2;
jum = outer_smooth + outer_arcioni;
outer_pp(logic_identical_pp) = jum(logic_identical_pp);

sum_1 = repmat(greenskalar_bg(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar_bg(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
jum = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;
outer_pp(logic_far_pp) = jum(logic_far_pp);

%% + - D_matrix region 1

sum_1 = K_2_1_minus(SegitigaPlus,:,:)./(4*pi)-k_bg.^2.*K_2_2_minus(SegitigaPlus,:,:)./(8*pi)+...
    repmat(greensmooth_bg(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_minus(SegitigaPlus,:)./(4*pi)-k_bg.^2.*K_1_2_minus(SegitigaPlus,:)./(8*pi)+...
    greensmooth_bg(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
outer_pm = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;

int_smooth_div = greensmooth_bg(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
int_smooth = repmat(greensmooth_bg(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;
Cmn = bsxfun(@times,div_f_plus,div_f_minus.')./4;
int_arcioni_div = bsxfun(@times,div_f_plus,div_f_minus.').*I_1_pm./(4*pi);
int_arcioni = Cmn.*I_2_pm./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k_bg.^2;
jum = outer_smooth + outer_arcioni;
outer_pm(logic_identical_pm) = jum(logic_identical_pm);

sum_1 = repmat(greenskalar_bg(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar_bg(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
jum = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k_bg.^2;
outer_pm(logic_far_pm) = jum(logic_far_pm);

%% - + D_matrix region 1

sum_1 = K_2_1_plus(SegitigaMinus,:,:)./(4*pi)-k_bg.^2.*K_2_2_plus(SegitigaMinus,:,:)./(8*pi)+...
    repmat(greensmooth_bg(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_plus(SegitigaMinus,:)./(4*pi)-k_bg.^2.*K_1_2_plus(SegitigaMinus,:)./(8*pi)+...
    greensmooth_bg(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
outer_mp = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;

int_smooth_div = greensmooth_bg(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
int_smooth = repmat(greensmooth_bg(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;
Cmn = bsxfun(@times,div_f_minus,div_f_plus.')./4;
int_arcioni_div = bsxfun(@times,div_f_minus,div_f_plus.').*I_1_mp./(4*pi);
int_arcioni = Cmn.*I_2_mp./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k_bg.^2;
jum = outer_smooth + outer_arcioni;
outer_mp(logic_identical_mp) = jum(logic_identical_mp);

sum_1 = repmat(greenskalar_bg(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar_bg(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
jum = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;
outer_mp(logic_far_mp) = jum(logic_far_mp);

%% - - D_matrix region 1

sum_1 = K_2_1_minus(SegitigaMinus,:,:)./(4*pi)-k_bg.^2.*K_2_2_minus(SegitigaMinus,:,:)./(8*pi)+...
    repmat(greensmooth_bg(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_minus(SegitigaMinus,:)./(4*pi)-k_bg.^2.*K_1_2_minus(SegitigaMinus,:)./(8*pi)+...
    greensmooth_bg(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
outer_mm = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;

int_smooth_div = greensmooth_bg(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
int_smooth = repmat(greensmooth_bg(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;
Cmn = bsxfun(@times,div_f_minus,div_f_minus.')./4;
int_arcioni_div = bsxfun(@times,div_f_minus,div_f_minus.').*I_1_mm./(4*pi);
int_arcioni = Cmn.*I_2_mm./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k_bg.^2;
jum = outer_smooth + outer_arcioni;
outer_mm(logic_identical_mm) = jum(logic_identical_mm);

sum_1 = repmat(greenskalar_bg(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar_bg(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
jum = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k_bg.^2;
outer_mm(logic_far_mm) = jum(logic_far_mm);

am_1 = omega.*mu_0.*(outer_pp+outer_pm+outer_mp+outer_mm)./(1i);

%% + + D_matrix region 2

sum_1 = K_2_1_plus(SegitigaPlus,:,:)./(4*pi)-k.^2.*K_2_2_plus(SegitigaPlus,:,:)./(8*pi)+...
    repmat(greensmooth(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_plus(SegitigaPlus,:)./(4*pi)-k.^2.*K_1_2_plus(SegitigaPlus,:)./(8*pi)+...
    greensmooth(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
outer_pp = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;

int_smooth_div = greensmooth(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
int_smooth = repmat(greensmooth(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;
Cmn = bsxfun(@times,div_f_plus,div_f_plus.')./4;
int_arcioni_div = bsxfun(@times,div_f_plus,div_f_plus.').*I_1_pp./(4*pi);
int_arcioni = Cmn.*I_2_pp./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k.^2;
jum = outer_smooth + outer_arcioni;
outer_pp(logic_identical_pp) = jum(logic_identical_pp);

sum_1 = repmat(greenskalar(SegitigaPlus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar(SegitigaPlus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
jum = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;
outer_pp(logic_far_pp) = jum(logic_far_pp);

%% + - D_matrix region 2

sum_1 = K_2_1_minus(SegitigaPlus,:,:)./(4*pi)-k.^2.*K_2_2_minus(SegitigaPlus,:,:)./(8*pi)+...
    repmat(greensmooth(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_minus(SegitigaPlus,:)./(4*pi)-k.^2.*K_1_2_minus(SegitigaPlus,:)./(8*pi)+...
    greensmooth(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
outer_pm = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;

int_smooth_div = greensmooth(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
int_smooth = repmat(greensmooth(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;
Cmn = bsxfun(@times,div_f_plus,div_f_minus.')./4;
int_arcioni_div = bsxfun(@times,div_f_plus,div_f_minus.').*I_1_pm./(4*pi);
int_arcioni = Cmn.*I_2_pm./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k.^2;
jum = outer_smooth + outer_arcioni;
outer_pm(logic_identical_pm) = jum(logic_identical_pm);

sum_1 = repmat(greenskalar(SegitigaPlus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar(SegitigaPlus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
jum = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi)...
    - sum_2.*repmat(div_f_plus.*luas(SegitigaPlus),1,jumSisi)./k.^2;
outer_pm(logic_far_pm) = jum(logic_far_pm);

%% - + D_matrix region 2

sum_1 = K_2_1_plus(SegitigaMinus,:,:)./(4*pi)-k.^2.*K_2_2_plus(SegitigaMinus,:,:)./(8*pi)+...
    repmat(greensmooth(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_plus(SegitigaMinus,:)./(4*pi)-k.^2.*K_1_2_plus(SegitigaMinus,:)./(8*pi)+...
    greensmooth(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
outer_mp = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;

int_smooth_div = greensmooth(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
int_smooth = repmat(greensmooth(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;
Cmn = bsxfun(@times,div_f_minus,div_f_plus.')./4;
int_arcioni_div = bsxfun(@times,div_f_minus,div_f_plus.').*I_1_mp./(4*pi);
int_arcioni = Cmn.*I_2_mp./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k.^2;
jum = outer_smooth + outer_arcioni;
outer_mp(logic_identical_mp) = jum(logic_identical_mp);

sum_1 = repmat(greenskalar(SegitigaMinus,SegitigaPlus).*repmat(luas(SegitigaPlus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_plus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar(SegitigaMinus,SegitigaPlus).*repmat(div_f_plus.'.*luas(SegitigaPlus).',jumSisi,1);
jum = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;
outer_mp(logic_far_mp) = jum(logic_far_mp);

%% - - D_matrix region 2

sum_1 = K_2_1_minus(SegitigaMinus,:,:)./(4*pi)-k.^2.*K_2_2_minus(SegitigaMinus,:,:)./(8*pi)+...
    repmat(greensmooth(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = K_1_1_minus(SegitigaMinus,:)./(4*pi)-k.^2.*K_1_2_minus(SegitigaMinus,:)./(8*pi)+...
    greensmooth(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
outer_mm = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;

int_smooth_div = greensmooth(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
int_smooth = repmat(greensmooth(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
outer_smooth = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),int_smooth,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - int_smooth_div.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;
Cmn = bsxfun(@times,div_f_minus,div_f_minus.')./4;
int_arcioni_div = bsxfun(@times,div_f_minus,div_f_minus.').*I_1_mm./(4*pi);
int_arcioni = Cmn.*I_2_mm./(4*pi);
outer_arcioni = int_arcioni - int_arcioni_div./k.^2;
jum = outer_smooth + outer_arcioni;
outer_mm(logic_identical_mm) = jum(logic_identical_mm);

sum_1 = repmat(greenskalar(SegitigaMinus,SegitigaMinus).*repmat(luas(SegitigaMinus).',jumSisi,1),1,1,3).*...
    permute(repmat(f_minus.',1,1,jumSisi),[3,2,1]);
sum_2 = greenskalar(SegitigaMinus,SegitigaMinus).*repmat(div_f_minus.'.*luas(SegitigaMinus).',jumSisi,1);
jum = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi)...
    - sum_2.*repmat(div_f_minus.*luas(SegitigaMinus),1,jumSisi)./k.^2;
outer_mm(logic_far_mm) = jum(logic_far_mm);

am_2 = omega.*mu_0.*(outer_pp+outer_pm+outer_mp+outer_mm)./(1i);

%% + + K_matrix region 1

sum_1 = K_4_1_plus(SegitigaPlus,:,:)./(4*pi)-k_bg.^2.*K_4_2_plus(SegitigaPlus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source_bg(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaPlus).',jumSisi,1,3);
outer_pp = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_pp(logic_identical_pp) = jum(logic_identical_pp);

% sum_1 = -k_bg.^2.*K_4_2_plus(SegitigaPlus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source_bg(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% inner_ssl_1 = pre_K_2_1_plus_SSL_1(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_plus_SSL_2(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_plus_SSL_3(SegitigaPlus,:,:)./(4*pi);
% Cmn = -bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaPlus),luas(SegitigaPlus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaPlus).',jumSisi,1)+....
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaPlus).',jumSisi,1)+....
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaPlus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_pp(logic_adjacent_pp) = jum(logic_adjacent_pp);

sum_1 = cross(gradientgreen_source_bg(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaPlus).',jumSisi,1,3);
jum2 = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
outer_pp(logic_far_pp) = jum2(logic_far_pp);

%% + - K_matrix region 1

sum_1 = K_4_1_minus(SegitigaPlus,:,:)./(4*pi)-k_bg.^2.*K_4_2_minus(SegitigaPlus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source_bg(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaMinus).',jumSisi,1,3);
outer_pm = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_pm(logic_identical_pm) = jum(logic_identical_pm);

% sum_1 = -k_bg.^2.*K_4_2_minus(SegitigaPlus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source_bg(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_plus_SSL_1(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_plus_SSL_2(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_plus_SSL_3(SegitigaMinus,:,:)./(4*pi);
% Cmn = bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaPlus),luas(SegitigaMinus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaMinus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_pm(logic_adjacent_pm) = jum(logic_adjacent_pm);

sum_1 = cross(gradientgreen_source_bg(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaMinus).',jumSisi,1,3);
jum2 = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
outer_pm(logic_far_pm) = jum2(logic_far_pm);

%% - + K_matrix region 1

sum_1 = K_4_1_plus(SegitigaMinus,:,:)./(4*pi)-k_bg.^2.*K_4_2_plus(SegitigaMinus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source_bg(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaPlus).',jumSisi,1,3);
outer_mp = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_mp(logic_identical_mp) = jum(logic_identical_mp);

% sum_1 = -k_bg.^2.*K_4_2_plus(SegitigaMinus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source_bg(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_minus_SSL_1(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_minus_SSL_2(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_minus_SSL_3(SegitigaPlus,:,:)./(4*pi);
% Cmn = bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaMinus),luas(SegitigaPlus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaPlus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaPlus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaPlus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_mp(logic_adjacent_mp) = jum(logic_adjacent_mp);

% sum_1 = cross(gradientgreen_source_bg(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% outer_mp(logic_far_mp) = jum2(logic_far_mp);

%% - - K_matrix region 1

sum_1 = K_4_1_minus(SegitigaMinus,:,:)./(4*pi)-k_bg.^2.*K_4_2_minus(SegitigaMinus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source_bg(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaMinus).',jumSisi,1,3);
outer_mm = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_mm(logic_identical_mm) = jum(logic_identical_mm);

% sum_1 = -k_bg.^2.*K_4_2_minus(SegitigaMinus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source_bg(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_minus_SSL_1(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_minus_SSL_2(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_minus_SSL_3(SegitigaMinus,:,:)./(4*pi);
% Cmn = -bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaMinus),luas(SegitigaMinus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 =cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaMinus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_mm(logic_adjacent_mm) = jum(logic_adjacent_mm);

% sum_1 = cross(gradientgreen_source_bg(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% outer_mm(logic_far_mm) = jum2(logic_far_mm);

bm_1 = outer_pp+outer_pm+outer_mp+outer_mm;

%% + + K_matrix region 2

sum_1 = K_4_1_plus(SegitigaPlus,:,:)./(4*pi)-k.^2.*K_4_2_plus(SegitigaPlus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaPlus).',jumSisi,1,3);
outer_pp = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_pp(logic_identical_pp) = jum(logic_identical_pp);

% sum_1 = -k.^2.*K_4_2_plus(SegitigaPlus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_plus_SSL_1(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_plus_SSL_2(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_plus_SSL_3(SegitigaPlus,:,:)./(4*pi);
% Cmn = -bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaPlus),luas(SegitigaPlus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaPlus).',jumSisi,1)+....
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaPlus).',jumSisi,1)+....
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaPlus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_pp(logic_adjacent_pp) = jum(logic_adjacent_pp);

% sum_1 = cross(gradientgreen_source(SegitigaPlus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% outer_pp(logic_far_pp) = jum2(logic_far_pp);

%% + - K_matrix region 2

sum_1 = K_4_1_minus(SegitigaPlus,:,:)./(4*pi)-k.^2.*K_4_2_minus(SegitigaPlus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaMinus).',jumSisi,1,3);
outer_pm = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_pm(logic_identical_pm) = jum(logic_identical_pm);

% sum_1 = -k.^2.*K_4_2_minus(SegitigaPlus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_plus_SSL_1(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_plus_SSL_2(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_plus_SSL_3(SegitigaMinus,:,:)./(4*pi);
% Cmn = bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaPlus),luas(SegitigaMinus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_plus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaMinus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_pm(logic_adjacent_pm) = jum(logic_adjacent_pm);

% sum_1 = cross(gradientgreen_source(SegitigaPlus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_plus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaPlus),1,jumSisi);
% outer_pm(logic_far_pm) = jum2(logic_far_pm);

%% - + K_matrix region 2

sum_1 = K_4_1_plus(SegitigaMinus,:,:)./(4*pi)-k.^2.*K_4_2_plus(SegitigaMinus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaPlus).',jumSisi,1,3);
outer_mp = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_mp(logic_identical_mp) = jum(logic_identical_mp);

% sum_1 = -k.^2.*K_4_2_plus(SegitigaMinus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_minus_SSL_1(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_minus_SSL_2(SegitigaPlus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_minus_SSL_3(SegitigaPlus,:,:)./(4*pi);
% Cmn = bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaMinus),luas(SegitigaPlus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_plus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaPlus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaPlus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaPlus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaPlus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_mp(logic_adjacent_mp) = jum(logic_adjacent_mp);

% sum_1 = cross(gradientgreen_source(SegitigaMinus,SegitigaPlus,:),permute(repmat(f_plus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaPlus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% outer_mp(logic_far_mp) = jum2(logic_far_mp);

%% - - K_matrix region 2

sum_1 = K_4_1_minus(SegitigaMinus,:,:)./(4*pi)-k.^2.*K_4_2_minus(SegitigaMinus,:,:)./(8*pi)+...
    cross(gradientgreensmooth_source(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
    repmat(luas(SegitigaMinus).',jumSisi,1,3);
outer_mm = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);

jum = zeros(jumSisi,jumSisi);
outer_mm(logic_identical_mm) = jum(logic_identical_mm);

% sum_1 = -k.^2.*K_4_2_minus(SegitigaMinus,:,:)./(8*pi)+...
%     cross(gradientgreensmooth_source(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% outer_std = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% %ssl part
% inner_ssl_1 = pre_K_2_1_minus_SSL_1(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_2 = pre_K_2_1_minus_SSL_2(SegitigaMinus,:,:)./(4*pi);
% inner_ssl_3 = pre_K_2_1_minus_SSL_3(SegitigaMinus,:,:)./(4*pi);
% Cmn = -bsxfun(@times,PanjangSisi,PanjangSisi.')./(4.*bsxfun(@times,luas(SegitigaMinus),luas(SegitigaMinus).'));
% ssl_1 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi1(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_2 =cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi2(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% ssl_3 = cross(reshape(repmat(FreeVertex_minus,jumSisi,1,1),jumSisi,jumSisi,3)-...
%     permute(repmat(FreeVertex_minus.',1,1,jumSisi),[3,2,1]),...
%     permute(reshape(repmat(unit_M_sisi3(SegitigaMinus).',3,jumSisi,1),3,jumSisi,jumSisi),[3,2,1]));
% outer_ssl = Cmn.*(dot(ssl_1,inner_ssl_1,3).*repmat(l1(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_2,inner_ssl_2,3).*repmat(l2(SegitigaMinus).',jumSisi,1)+...
%     dot(ssl_3,inner_ssl_3,3).*repmat(l3(SegitigaMinus).',jumSisi,1));
% jum = outer_std + outer_ssl;
% outer_mm(logic_adjacent_mm) = jum(logic_adjacent_mm);

% sum_1 = cross(gradientgreen_source(SegitigaMinus,SegitigaMinus,:),permute(repmat(f_minus.',1,1,jumSisi),[3,2,1])).*...
%     repmat(luas(SegitigaMinus).',jumSisi,1,3);
% jum2 = dot(reshape(repmat(f_minus,jumSisi,1,1),jumSisi,jumSisi,3),sum_1,3).*repmat(luas(SegitigaMinus),1,jumSisi);
% outer_mm(logic_far_mm) = jum2(logic_far_mm);

bm_2 = outer_pp+outer_pm+outer_mp+outer_mm;

%% Surface current (clear)

A = [(am_1+am_2) -(bm_1+bm_2) ; (bm_1+bm_2) (am_1/imp_bg.^2+am_2/imp.^2)];

% cur_mat = A\q.';

% LU fact
[L,U,P] = lu(A);
y = L\(P*q.');
cur_mat = U\y;

save surface_current_coef cur_mat

for a = 1 : jumSegitiga
    i=[0 0 0];
    q=[0 0 0];
    q2=0;
    i2=0;
    for m=1:jumSisi
        if(SegitigaPlus(m)==a)
            i=i+cur_mat(m).*f_plus(m,:);
            i2=i2+cur_mat(m).*div_f_plus(m,:);
        end
        if(SegitigaMinus(m)==a)
            i=i+cur_mat(m).*f_minus(m,:);
            i2=i2+cur_mat(m).*div_f_minus(m,:);
        end
        if(SegitigaPlus(m)==a)
            q=q+cur_mat(m+jumSisi).*f_plus(m,:);
            q2=q2+cur_mat(m+jumSisi).*div_f_plus(m,:);
        end
        if(SegitigaMinus(m)==a)
            q=q+cur_mat(m+jumSisi).*f_minus(m,:);
            q2=q2+cur_mat(m+jumSisi).*div_f_minus(m,:);
        end
    end
    J(a,:)=i;
    J_(a,:)=i2;
    M(a,:)=q;
    M_(a,:)=q2;
    mag_J(a) = norm(i);
    mag_M(a) = norm(q);
end

save Current J...
        J_...
        M...
        M_...
        cur_mat

for a = 1 : jumSegitiga
    E_surf(a,:) = cross(faceNorm(a,:),reshape(M(a,:),1,3))-(1i)*J_(a)*faceNorm(a,:)./(omega*eps_0*eps_bg);
    H_surf(a,:) = -cross(faceNorm(a,:),reshape(J(a,:),1,3))-(1i)*M_(a)*faceNorm(a,:)./(omega*mu_0);
    E_inc(a,:) = E_0.*exp((1i)*dot(titik_tengah(a,:),kv_bg));
    H_inc(a,:) = cross(k_vec,E_0).*exp((1i)*dot(titik_tengah(a,:),kv_bg))./(mu_0*c);
    mag_E(a) = norm(E_surf(a,:));
    mag_H(a) = norm(H_surf(a,:));
end

for m=1:jumSegitiga
    N=t(m,1:3);
    X(:,m) = p(N,1);
    Y(:,m) = p(N,2);
    Z(:,m) = p(N,3);
end

D = repmat(mag_E,3,1);
D2 = repmat(mag_H,3,1);

figure()
h=fill3(X, Y, Z, D);
title('Electric surface field');
shading interp;
colorbar;
colormap(jet);
rotate3d
xlabel('x')
ylabel('y')
zlabel('z')

figure()
h=fill3(X, Y, Z, D2);
title('Magnetic surface field');
shading interp;
colorbar;
colormap(jet);
rotate3d
xlabel('x')
ylabel('y')
zlabel('z')