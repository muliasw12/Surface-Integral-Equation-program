clearvars -except tStart
load parameter_segitiga

f_plus = repmat(PanjangSisi,1,3).*rho_plus./(2*repmat(luas(SegitigaPlus),1,3));
f_minus = -repmat(PanjangSisi,1,3).*rho_minus./(2*repmat(luas(SegitigaMinus),1,3));
div_f_plus = PanjangSisi./luas(SegitigaPlus);
div_f_minus = -PanjangSisi./luas(SegitigaMinus);

[x,y] = meshgrid(1:jumSegitiga,1:jumSegitiga);
source_tri_vect_1 = reshape(titik_tengah(x,1),jumSegitiga,jumSegitiga);
source_tri_vect_2 = reshape(titik_tengah(x,2),jumSegitiga,jumSegitiga);
source_tri_vect_3 = reshape(titik_tengah(x,3),jumSegitiga,jumSegitiga);
obs_point_vect_1 = reshape(titik_tengah(y,1),jumSegitiga,jumSegitiga);
obs_point_vect_2 = reshape(titik_tengah(y,2),jumSegitiga,jumSegitiga);
obs_point_vect_3 = reshape(titik_tengah(y,3),jumSegitiga,jumSegitiga);
obs_point = reshape([obs_point_vect_1 obs_point_vect_2 obs_point_vect_3],...
    jumSegitiga,jumSegitiga,3);
clear obs_point_vect_1 obs_point_vect_2 obs_point_vect_3
source_tri = t(1:jumSegitiga,:);
q1 = p(source_tri(:,1),:);
q2 = p(source_tri(:,2),:);
q3 = p(source_tri(:,3),:);
q1_1 = reshape(q1(x,1),jumSegitiga,jumSegitiga);
q1_2 = reshape(q1(x,2),jumSegitiga,jumSegitiga);
q1_3 = reshape(q1(x,3),jumSegitiga,jumSegitiga);
q2_1 = reshape(q2(x,1),jumSegitiga,jumSegitiga);
q2_2 = reshape(q2(x,2),jumSegitiga,jumSegitiga);
q2_3 = reshape(q2(x,3),jumSegitiga,jumSegitiga);
q3_1 = reshape(q3(x,1),jumSegitiga,jumSegitiga);
q3_2 = reshape(q3(x,2),jumSegitiga,jumSegitiga);
q3_3 = reshape(q3(x,3),jumSegitiga,jumSegitiga);
q1 = reshape([q1_1 q1_2 q1_3],jumSegitiga,jumSegitiga,3);
q2 = reshape([q2_1 q2_2 q2_3],jumSegitiga,jumSegitiga,3);
q3 = reshape([q3_1 q3_2 q3_3],jumSegitiga,jumSegitiga,3);
clear q1_1 q1_2 q1_3 q2_1 q2_2 q2_3 q3_1 q3_2 q3_3
vec_sisi1 = q2-q1;
vec_sisi2 = q3-q2;
vec_sisi3 = q1-q3;
source_point1 = (q2+q1)./2;
source_point2 = (q3+q2)./2;
source_point3 = (q1+q3)./2;
source_point1_2 = q1;
source_point2_2 = q2;
source_point3_2 = q3;
R_sisi1_plus = obs_point-q2;
R_sisi1_plus = sqrt(R_sisi1_plus(:,:,1).^2+R_sisi1_plus(:,:,2).^2+R_sisi1_plus(:,:,3).^2);
R_sisi2_plus = obs_point-q3;
R_sisi2_plus = sqrt(R_sisi2_plus(:,:,1).^2+R_sisi2_plus(:,:,2).^2+R_sisi2_plus(:,:,3).^2);
R_sisi3_plus = obs_point-q1;
R_sisi3_plus = sqrt(R_sisi3_plus(:,:,1).^2+R_sisi3_plus(:,:,2).^2+R_sisi3_plus(:,:,3).^2);
R_sisi1_minus = R_sisi3_plus;
R_sisi2_minus = R_sisi1_plus;
R_sisi3_minus = R_sisi2_plus;

unit_edge_sisi1 = vec_sisi1./reshape(repmat(sqrt(vec_sisi1(:,:,1).^2+...
    vec_sisi1(:,:,2).^2+vec_sisi1(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);
unit_edge_sisi2 = vec_sisi2./reshape(repmat(sqrt(vec_sisi2(:,:,1).^2+...
    vec_sisi2(:,:,2).^2+vec_sisi2(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);
unit_edge_sisi3 = vec_sisi3./reshape(repmat(sqrt(vec_sisi3(:,:,1).^2+...
    vec_sisi3(:,:,2).^2+vec_sisi3(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);
S_sisi1_plus = dot(q2-obs_point,unit_edge_sisi1,3);
S_sisi2_plus = dot(q3-obs_point,unit_edge_sisi2,3);
S_sisi3_plus = dot(q1-obs_point,unit_edge_sisi3,3);
S_sisi1_minus = dot(q1-obs_point,unit_edge_sisi1,3);
S_sisi2_minus = dot(q2-obs_point,unit_edge_sisi2,3);
S_sisi3_minus = dot(q3-obs_point,unit_edge_sisi3,3);
faceNorm_1 = reshape(faceNorm(x,1),jumSegitiga,jumSegitiga);
faceNorm_2 = reshape(faceNorm(x,2),jumSegitiga,jumSegitiga);
faceNorm_3 = reshape(faceNorm(x,3),jumSegitiga,jumSegitiga);
faceNorm = reshape([faceNorm_1 faceNorm_2 faceNorm_3],jumSegitiga,jumSegitiga,3);
clear faceNorm_1 faceNorm_2 faceNorm_3

unit_M_sisi1 = cross(unit_edge_sisi1,faceNorm);
unit_M_sisi1 = unit_M_sisi1./reshape(repmat(sqrt(unit_M_sisi1(:,:,1).^2+...
    unit_M_sisi1(:,:,2).^2+unit_M_sisi1(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);
unit_M_sisi2 = cross(unit_edge_sisi2,faceNorm);
unit_M_sisi2 = unit_M_sisi2./reshape(repmat(sqrt(unit_M_sisi2(:,:,1).^2+...
    unit_M_sisi2(:,:,2).^2+unit_M_sisi2(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);
unit_M_sisi3 = cross(unit_edge_sisi3,faceNorm);
unit_M_sisi3 = unit_M_sisi3./reshape(repmat(sqrt(unit_M_sisi3(:,:,1).^2+...
    unit_M_sisi3(:,:,2).^2+unit_M_sisi3(:,:,3).^2),1,3),jumSegitiga,jumSegitiga,3);

W_0 = dot(faceNorm,obs_point-source_point1,3);
temp = zeros(jumSegitiga,jumSegitiga);
% W_0(logic_identical) = temp(logic_identical);
rho = obs_point - reshape(repmat(W_0,1,1,3),jumSegitiga,jumSegitiga,3).*faceNorm;
R0_sisi1_kuad = R_sisi1_plus.^2-S_sisi1_plus.^2;
R0_sisi2_kuad = R_sisi2_plus.^2-S_sisi2_plus.^2;
R0_sisi3_kuad = R_sisi3_plus.^2-S_sisi3_plus.^2;
T_0_1 = dot(obs_point-source_point1,unit_M_sisi1,3);
T_0_2 = dot(obs_point-source_point2,unit_M_sisi2,3);
T_0_3 = dot(obs_point-source_point3,unit_M_sisi3,3);
a_1 = q1-obs_point;
a_1 = a_1./reshape(repmat(sqrt(a_1(:,:,1).^2+a_1(:,:,2).^2+a_1(:,:,3).^2),...
    1,3),jumSegitiga,jumSegitiga,3);
a_2 = q2-obs_point;
a_2 = a_2./reshape(repmat(sqrt(a_2(:,:,1).^2+a_2(:,:,2).^2+a_2(:,:,3).^2),...
    1,3),jumSegitiga,jumSegitiga,3);
a_3 = q3-obs_point;
a_3 = a_3./reshape(repmat(sqrt(a_3(:,:,1).^2+a_3(:,:,2).^2+a_3(:,:,3).^2),...
    1,3),jumSegitiga,jumSegitiga,3);
X = 1 + dot(a_1,a_2,3) + dot(a_1,a_3,3) + dot(a_2,a_3,3);
Y = abs(dot(a_1,cross(a_2,a_3),3));

AngleExcess = zeros(jumSegitiga,jumSegitiga);
temp = 2.*atan2(Y,X);
AngleExcess(W_0>0)=temp(W_0>0);
temp = -2.*atan2(Y,X);
AngleExcess(W_0<0)=temp(W_0<0);

basic_line_1 = log((R_sisi1_minus-S_sisi1_minus)./(R_sisi1_plus-S_sisi1_plus));
temp = log((R_sisi1_plus+S_sisi1_plus)./(R_sisi1_minus+S_sisi1_minus));
logic_RS_1 = abs(R_sisi1_minus+S_sisi1_minus)>abs(R_sisi1_plus-S_sisi1_plus);
basic_line_1(logic_RS_1) = temp(logic_RS_1);

basic_line_2 = log((R_sisi2_minus-S_sisi2_minus)./(R_sisi2_plus-S_sisi2_plus));
temp = log((R_sisi2_plus+S_sisi2_plus)./(R_sisi2_minus+S_sisi2_minus));
logic_RS_2 = abs(R_sisi2_minus+S_sisi2_minus)>abs(R_sisi2_plus-S_sisi2_plus);
basic_line_2(logic_RS_2) = temp(logic_RS_2);

basic_line_3 = log((R_sisi3_minus-S_sisi3_minus)./(R_sisi3_plus-S_sisi3_plus));
temp = log((R_sisi3_plus+S_sisi3_plus)./(R_sisi3_minus+S_sisi3_minus));
logic_RS_3 = abs(R_sisi3_minus+S_sisi3_minus)>abs(R_sisi3_plus-S_sisi3_plus);
basic_line_3(logic_RS_3) = temp(logic_RS_3);

adv1_line_1 = R0_sisi1_kuad.*basic_line_1./2+(R_sisi1_plus.*S_sisi1_plus-...
    R_sisi1_minus.*S_sisi1_minus)./2;
adv1_line_2 = R0_sisi2_kuad.*basic_line_2./2+(R_sisi2_plus.*S_sisi2_plus-...
    R_sisi2_minus.*S_sisi2_minus)./2;
adv1_line_3 = R0_sisi3_kuad.*basic_line_3./2+(R_sisi3_plus.*S_sisi3_plus-...
    R_sisi3_minus.*S_sisi3_minus)./2;
adv2_line_1 = R0_sisi1_kuad.*adv1_line_1.*3./4+(R_sisi1_plus.^3.*S_sisi1_plus-...
    R_sisi1_minus.^3.*S_sisi1_minus)./4;
adv2_line_2 = R0_sisi2_kuad.*adv1_line_2.*3./4+(R_sisi2_plus.^3.*S_sisi2_plus-...
    R_sisi2_minus.^3.*S_sisi2_minus)./4;
adv2_line_3 = R0_sisi3_kuad.*adv1_line_3.*3./4+(R_sisi3_plus.^3.*S_sisi3_plus-...
    R_sisi3_minus.^3.*S_sisi3_minus)./4;

basic_surf = AngleExcess./W_0;
adv1_surf = -W_0.^2.*basic_surf-(T_0_1.*basic_line_1+T_0_2.*basic_line_2+...
    T_0_3.*basic_line_3);
adv2_surf = W_0.^2.*adv1_surf./3-(T_0_1.*adv1_line_1+T_0_2.*adv1_line_2+...
    T_0_3.*adv1_line_3)./3;
temp = -(T_0_1.*basic_line_1+T_0_2.*basic_line_2+T_0_3.*basic_line_3);
adv1_surf(W_0 == 0) = temp(W_0 == 0);
temp = -(T_0_1.*adv1_line_1+T_0_2.*adv1_line_2+T_0_3.*adv1_line_3)./3;
adv2_surf(W_0 == 0) = temp(W_0 == 0);

basic_I_nm = unit_M_sisi1.*reshape(repmat(basic_line_1,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi2.*reshape(repmat(basic_line_2,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi3.*reshape(repmat(basic_line_3,1,1,3),jumSegitiga,jumSegitiga,3);
adv1_I_nm = unit_M_sisi1.*reshape(repmat(adv1_line_1,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi2.*reshape(repmat(adv1_line_2,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi3.*reshape(repmat(adv1_line_3,1,1,3),jumSegitiga,jumSegitiga,3);
adv2_I_nm = unit_M_sisi1.*reshape(repmat(adv2_line_1,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi2.*reshape(repmat(adv2_line_2,1,1,3),jumSegitiga,jumSegitiga,3)+...
    unit_M_sisi3.*reshape(repmat(adv2_line_3,1,1,3),jumSegitiga,jumSegitiga,3);

[x,y] = meshgrid(1:jumSisi,1:jumSegitiga);

K_1_1_plus = repmat(div_f_plus.',jumSegitiga,1).*adv1_surf(:,SegitigaPlus);
K_1_1_minus = repmat(div_f_minus.',jumSegitiga,1).*adv1_surf(:,SegitigaMinus);
K_1_2_plus = repmat(div_f_plus.',jumSegitiga,1).*adv2_surf(:,SegitigaPlus);
K_1_2_minus = repmat(div_f_minus.',jumSegitiga,1).*adv2_surf(:,SegitigaMinus);

FreeVertex_plus_1 = reshape(FreeVertex_plus(x,1),jumSegitiga,jumSisi);
FreeVertex_plus_2 = reshape(FreeVertex_plus(x,2),jumSegitiga,jumSisi);
FreeVertex_plus_3 = reshape(FreeVertex_plus(x,3),jumSegitiga,jumSisi);
FreeVertex_plus = reshape([FreeVertex_plus_1 FreeVertex_plus_2 FreeVertex_plus_3],...
    jumSegitiga,jumSisi,3);
clear FreeVertex_plus_1 FreeVertex_plus_2 FreeVertex_plus_3
FreeVertex_minus_1 = reshape(FreeVertex_minus(x,1),jumSegitiga,jumSisi);
FreeVertex_minus_2 = reshape(FreeVertex_minus(x,2),jumSegitiga,jumSisi);
FreeVertex_minus_3 = reshape(FreeVertex_minus(x,3),jumSegitiga,jumSisi);
FreeVertex_minus = reshape([FreeVertex_minus_1 FreeVertex_minus_2 FreeVertex_minus_3],...
    jumSegitiga,jumSisi,3);
clear FreeVertex_minus_1 FreeVertex_minus_2 FreeVertex_minus_3

K_2_1_plus = repmat(div_f_plus.',jumSegitiga,1,3).*((rho(:,SegitigaPlus,:)-FreeVertex_plus).*...
    reshape(repmat(adv1_surf(:,SegitigaPlus),1,1,3),jumSegitiga,jumSisi,3)+...
    adv1_I_nm(:,SegitigaPlus,:))./2;
K_2_1_minus = repmat(div_f_minus.',jumSegitiga,1,3).*((rho(:,SegitigaMinus,:)-FreeVertex_minus).*...
    reshape(repmat(adv1_surf(:,SegitigaMinus),1,1,3),jumSegitiga,jumSisi,3)+...
    adv1_I_nm(:,SegitigaMinus,:))./2;
K_2_2_plus = repmat(div_f_plus.',jumSegitiga,1,3).*((rho(:,SegitigaPlus,:)-FreeVertex_plus).*...
    reshape(repmat(adv2_surf(:,SegitigaPlus),1,1,3),jumSegitiga,jumSisi,3)+...
    adv2_I_nm(:,SegitigaPlus,:)./3)./2;
K_2_2_minus = repmat(div_f_minus.',jumSegitiga,1,3).*((rho(:,SegitigaMinus,:)-FreeVertex_minus).*...
    reshape(repmat(adv2_surf(:,SegitigaMinus),1,1,3),jumSegitiga,jumSisi,3)+...
    adv2_I_nm(:,SegitigaMinus,:)./3)./2;

logic_W_plus = repmat(W_0(:,SegitigaPlus),1,1,3) == 0;
logic_W_minus = repmat(W_0(:,SegitigaMinus),1,1,3) == 0;
K_3_1_plus = basic_I_nm(:,SegitigaPlus,:)-(repmat(-W_0(:,SegitigaPlus).*...
    basic_surf(:,SegitigaPlus),1,1,3).*faceNorm(:,SegitigaPlus,:));
K_3_2_plus = adv1_I_nm(:,SegitigaPlus,:)-(repmat(W_0(:,SegitigaPlus).*...
    adv1_surf(:,SegitigaPlus),1,1,3).*faceNorm(:,SegitigaPlus,:));
K_3_1_minus = basic_I_nm(:,SegitigaMinus,:)-(repmat(-W_0(:,SegitigaMinus).*...
    basic_surf(:,SegitigaMinus),1,1,3).*faceNorm(:,SegitigaMinus,:));
K_3_2_minus = adv1_I_nm(:,SegitigaMinus,:)-(repmat(W_0(:,SegitigaMinus).*...
    adv1_surf(:,SegitigaMinus),1,1,3).*faceNorm(:,SegitigaMinus,:));
temp = basic_I_nm(:,SegitigaPlus,:)+repmat(AngleExcess(:,SegitigaPlus),1,1,3).*...
    faceNorm(:,SegitigaPlus,:);
K_3_1_plus(logic_W_plus) = temp(logic_W_plus);
temp = adv1_I_nm(:,SegitigaPlus,:).*zeros(jumSegitiga,jumSisi,3);
K_3_2_plus(logic_W_plus) = temp(logic_W_plus);
temp = basic_I_nm(:,SegitigaMinus,:)+repmat(AngleExcess(:,SegitigaMinus),1,1,3).*...
    faceNorm(:,SegitigaMinus,:);
K_3_1_minus(logic_W_minus) = temp(logic_W_minus);
temp = adv1_I_nm(:,SegitigaMinus,:).*zeros(jumSegitiga,jumSisi,3);
K_3_2_minus(logic_W_minus) = temp(logic_W_minus);

logic_K_4_plus = y == SegitigaPlus(x);
logic_K_4_minus = y == SegitigaMinus(x);
obs_point_vect_1 = reshape(titik_tengah(y,1),jumSegitiga,jumSisi);
obs_point_vect_2 = reshape(titik_tengah(y,2),jumSegitiga,jumSisi);
obs_point_vect_3 = reshape(titik_tengah(y,3),jumSegitiga,jumSisi);
obs_point = reshape([obs_point_vect_1 obs_point_vect_2 obs_point_vect_3],...
    jumSegitiga,jumSisi,3);
clear obs_point_vect_1 obs_point_vect_2 obs_point_vect_3
K_4_1_plus = repmat(div_f_plus.',jumSegitiga,1,3).*...
    cross(K_3_1_plus,obs_point-FreeVertex_plus)./2;
K_4_2_plus = repmat(div_f_plus.',jumSegitiga,1,3).*...
    cross(K_3_2_plus,obs_point-FreeVertex_plus)./2;
K_4_1_minus = repmat(div_f_minus.',jumSegitiga,1,3).*...
    cross(K_3_1_minus,obs_point-FreeVertex_minus)./2;
K_4_2_minus = repmat(div_f_minus.',jumSegitiga,1,3).*...
    cross(K_3_2_minus,obs_point-FreeVertex_minus)./2;
temp = zeros(jumSegitiga,jumSisi,3);
K_4_1_plus(logic_K_4_plus) = temp(logic_K_4_plus);
K_4_2_plus(logic_K_4_plus) = temp(logic_K_4_plus);
K_4_1_minus(logic_K_4_minus) = temp(logic_K_4_minus);
K_4_2_minus(logic_K_4_minus) = temp(logic_K_4_minus);

save singular_integral K_4_1_plus...
   K_4_1_minus...
   K_4_2_plus...
   K_4_2_minus...
   K_2_1_plus...
   K_2_1_minus...
   K_2_2_plus...
   K_2_2_minus...
   K_1_1_plus...
   K_1_1_minus...
   K_1_2_plus...
   K_1_2_minus...
   unit_M_sisi1...
   unit_M_sisi2...
   unit_M_sisi3
