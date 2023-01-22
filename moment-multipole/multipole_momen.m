clearvars -except avg_bench
close all
load csca_1cube % csca from SIE

%%%%% input %%%%%
wl_interp = (300:2:800).'; % make sure that wl is consistent with the csca
N_multipole = 3;
%%%%%%%%%%%%%%%%%

%% main program
eV_interp = 1240./wl_interp;
h_bar_eV = 6.582119e-16;
omega = eV_interp./h_bar_eV;
lamda = 1240*10^-9./eV_interp;
k = 2*pi().*sqrt(eps_bg)./lamda;
eps_0 = 8.854e-12;
mu_0 = 1.257e-6;
imp_bg = sqrt(mu_0/(eps_0*eps_bg));
c  = 1/sqrt(eps_0*mu_0);
freq = c./lamda.*1e-12;
E_inc = 1;
titik_x = titik_tengah_ff(:,1);
titik_y = titik_tengah_ff(:,2);
titik_z = titik_tengah_ff(:,3);

r = sqrt(titik_x.^2+titik_y.^2+titik_z.^2);
r2 = sqrt(titik_x.^2+titik_y.^2);
sin_theta = r2./r;
cos_theta = titik_z./r;
tan_theta = r2./titik_z;
sin_phi = titik_y./r2;
cos_phi = titik_x./r2;
tan_phi = sin_phi./cos_phi;
phi = atan2(sin_phi,cos_phi);
theta = atan2(sin_theta,cos_theta);
temp = zeros(length(phi),1);
temp2 = zeros(length(phi),3);

% hai = 123;
for hai = 1 : size(E_farfield_sca,1)
    E_sca_rho = sin_theta.*cos_phi.*E_farfield_sca(hai,:,1).'+...
        sin_theta.*sin_phi.*E_farfield_sca(hai,:,2).'+...
        cos_theta.*E_farfield_sca(hai,:,3).';
    E_sca_theta = cos_theta.*cos_phi.*E_farfield_sca(hai,:,1).'+...
        cos_theta.*sin_phi.*E_farfield_sca(hai,:,2).'-...
        sin_theta.*E_farfield_sca(hai,:,3).';
    E_sca_phi = -sin_phi.*E_farfield_sca(hai,:,1).'+cos_phi.*E_farfield_sca(hai,:,2).';
    E_sca_sph = [E_sca_rho E_sca_theta E_sca_phi];
    
    H_sca_rho = sin_theta.*cos_phi.*H_farfield_sca(hai,:,1).'+...
        sin_theta.*sin_phi.*H_farfield_sca(hai,:,2).'+...
        cos_theta.*H_farfield_sca(hai,:,3).';
    H_sca_theta = cos_theta.*cos_phi.*H_farfield_sca(hai,:,1).'+...
        cos_theta.*sin_phi.*H_farfield_sca(hai,:,2).'-...
        sin_theta.*H_farfield_sca(hai,:,3).';
    H_sca_phi = -sin_phi.*H_farfield_sca(hai,:,1).'+cos_phi.*H_farfield_sca(hai,:,2).';
    H_sca_sph = [H_sca_rho H_sca_theta H_sca_phi];
    
    %% vec_Mo and vec_Ne
   
    m_ = -N_multipole:N_multipole;
    jum_m = length(m_);
    jum_n = (jum_m-1)/2;
    
    for n = 1 : jum_n
        for u = 1 : jum_m
            m = m_(u);
            for a = 1 : length(phi)
                var_r = r(a);
                var_cos_theta = cos_theta(a);
                var_sin_theta = sin_theta(a);
                var_phi = phi(a);
                var_rho = k(hai)*var_r;
                
                f_sphhankel_1 = sqrt(pi./(2*var_rho)).*besselh(n+0.5,var_rho);
                d_sphhankel_1 = f_sphhankel_1+var_rho.*(-f_sphhankel_1./(2.*var_rho)+...
                    (sqrt(pi./(2*var_rho)).*besselh(n-1+0.5,var_rho)-...
                    sqrt(pi./(2*var_rho)).*besselh(n+1+0.5,var_rho))./2);
                
                f_ass_Legendre = legendreP_mod(n,m,var_cos_theta); % m = 0,1,2,3,...
                d_Legendre = -1./var_sin_theta.*((1+n).*var_cos_theta.*f_ass_Legendre+...
                    (-1+m-n).*legendreP_mod(n+1,m,var_cos_theta));
                
                vec_M = [0 (1i).*m.*f_ass_Legendre./var_sin_theta.*...
                    f_sphhankel_1.*exp((1i).*m.*var_phi) ...
                    -d_Legendre.*f_sphhankel_1.*exp((1i).*m.*var_phi)];
                vec_N = [n.*(n+1).*f_ass_Legendre.*f_sphhankel_1.*...
                    exp((1i).*m.*var_phi)./var_rho d_Legendre.*...
                    d_sphhankel_1.*exp((1i).*m.*var_phi)./var_rho ...
                    (1i).*m.*f_ass_Legendre./var_sin_theta.*...
                    d_sphhankel_1.*exp((1i).*m.*var_phi)./var_rho;];
                atas_M(n,u,a) = dot(E_sca_sph(a,:),conj(vec_M));
                atas_N(n,u,a) = dot(E_sca_sph(a,:),conj(vec_N));
                bawah_M(n,u,a) = norm(vec_M).^2;
                bawah_N(n,u,a) = norm(vec_N).^2;
            end
        end
    end
    sum_atas_M = sum(atas_M,3);
    sum_atas_N = sum(atas_N,3);
    sum_bawah_M = sum(bawah_M,3);
    sum_bawah_N = sum(bawah_N,3);
    
    for n = 1 : jum_n
        for u = 1 : jum_m
            m =  m_(u);
            if (n+m)<0||(n-m)<0
                continue;
            else
                En= (1i).^(n+2.*m-1).*sqrt((2.*n+1).*factorial(n-m)./...
                    factorial(n+m))./(2.*sqrt(pi()));
                b_mn(hai,n,u) = sum_atas_M(n,u)./(sum_bawah_M(n,u).*En);
                a_mn(hai,n,u) = sum_atas_N(n,u)./(sum_bawah_N(n,u).*En);
            end
        end
    end 
    prog_2 = hai*100/size(E_farfield_sca,1)
end

%% C_sca benchmarking (SIE-inversed)
lamda2 = lamda*10^9;
for a = 1 : length(lamda) % c_sca SIE
    for b = 1 : jumSegitiga_ff
        sum_sca(b,:) = luas_ff(b).*dot(faceNorm_ff(b,:),...
            real(reshape(cross(E_farfield_sca(a,b,:),conj(H_farfield_sca(a,b,:))),1,3)/2));
    end
    P_sca = sum(sum_sca);
    illumination = 1/(2*imp_bg);
    c_sca_num(a,:) = P_sca/illumination;
    a
end

for i = 1 : jum_n % c_sca num
    for j = 1 : jum_m
        lala2(:,i,j) = i*(i+1).*(abs(a_mn(:,i,j)).^2+abs(b_mn(:,i,j)).^2);
    end
end
c_sca_bohren_total = sum(sum(lala2,3),2)./k.^2;

%% c_sca contributions (12)

for i = 1 : jum_n
    for j = 1 : jum_m
        lala_E(:,i,j) = i*(i+1).*abs(a_mn(:,i,j)).^2;
        lala_M(:,i,j) = i*(i+1).*abs(b_mn(:,i,j)).^2;
    end
end

col = hsv(jum_n);
c_sca_cont_E =  sum(lala_E,3)./k.^2;
c_sca_cont_M =  sum(lala_M,3)./k.^2;

%% bench c_sca komponen dipol dan kuadrupol

b_n_dipol = squeeze(b_mn(:,1,:));
a_n_dipol = squeeze(a_mn(:,1,:));
b_n_kuadrupol = squeeze(b_mn(:,2,:));
a_n_kuadrupol = squeeze(a_mn(:,2,:));

for i = 1 : jum_m
    lala_dipol(:,i) = j*(j+1).*(abs(a_n_dipol(:,i)).^2+abs(b_n_dipol(:,i)).^2);
    lala_kuadrupol(:,i) = j*(j+1).*(abs(a_n_kuadrupol(:,i)).^2+abs(b_n_kuadrupol(:,i)).^2);
end
c_sca_dipol = sum(lala_dipol,2)./k.^2;
c_sca_kuadrupol = sum(lala_kuadrupol,2)./k.^2;

ind_m2 = (jum_m+1)/2-2;
ind_m1 = (jum_m+1)/2-1;
ind_0  = (jum_m+1)/2;
ind_p1 = (jum_m+1)/2+1;
ind_p2 = (jum_m+1)/2+2;

C_0 = sqrt(6*pi()).*(1i)./(c.*imp_bg.*k);
D_0 = 6.*sqrt(30*pi())./((1i).*imp_bg.*c.*k.^2);

p_x = C_0.*(a_n_dipol(:,ind_p1)-a_n_dipol(:,ind_m1));
p_y = C_0.*(1i).*(a_n_dipol(:,ind_p1)+a_n_dipol(:,ind_m1));
p_z = -C_0.*sqrt(2).*a_n_dipol(:,ind_0);
p = p_x.^2+p_y.^2+p_z.^2;

m_x = c.*C_0.*(b_n_dipol(:,ind_p1)-b_n_dipol(:,ind_m1));
m_y = c.*C_0.*(1i).*(b_n_dipol(:,ind_p1)+b_n_dipol(:,ind_m1));
m_z = -c.*C_0.*sqrt(2).*b_n_dipol(:,ind_0);
m = m_x.^2+m_y.^2+m_z.^2;

qE_xx = D_0.*((1i).*(a_n_kuadrupol(:,ind_p2)+a_n_kuadrupol(:,ind_m2))-...
    (1i).*sqrt(6).*a_n_kuadrupol(:,ind_0)./2);
qE_yy = D_0.*(-(1i).*(a_n_kuadrupol(:,ind_p2)+a_n_kuadrupol(:,ind_m2))-...
    (1i).*sqrt(6).*a_n_kuadrupol(:,ind_0)./2);
qE_zz = D_0.*(1i).*sqrt(6).*a_n_kuadrupol(:,ind_0);
qE_xy = D_0.*(a_n_kuadrupol(:,ind_m2)-a_n_kuadrupol(:,ind_p2));
qE_xz = D_0.*(1i).*(a_n_kuadrupol(:,ind_m1)-a_n_kuadrupol(:,ind_p1));
qE_yz = D_0.*(a_n_kuadrupol(:,ind_m1)+a_n_kuadrupol(:,ind_p1));
qE_yx = qE_xy;
qE_zx = qE_xz;
qE_zy = qE_yz;

qM_xx = D_0.*c.*((1i).*(b_n_kuadrupol(:,ind_p2)+b_n_kuadrupol(:,ind_m2))-...
    (1i).*sqrt(6).*b_n_kuadrupol(:,ind_0)./2);
qM_yy = D_0.*c.*(-(1i).*(b_n_kuadrupol(:,ind_p2)+b_n_kuadrupol(:,ind_m2))-...
    (1i).*sqrt(6).*b_n_kuadrupol(:,ind_0)./2);
qM_zz = D_0.*c.*(1i).*sqrt(6).*b_n_kuadrupol(:,ind_0);
qM_xy = D_0.*c.*(b_n_kuadrupol(:,ind_m2)-b_n_kuadrupol(:,ind_p2));
qM_xz = D_0.*c.*(1i).*(b_n_kuadrupol(:,ind_m1)-b_n_kuadrupol(:,ind_p1));
qM_yz = D_0.*c.*(b_n_kuadrupol(:,ind_m1)+b_n_kuadrupol(:,ind_p1));
qM_yx = qM_xy;
qM_zx = qM_xz;
qM_zy = qM_yz;

c_sca_p = (1./(k.^2.*C_0.^2)).*(abs((1i).*p_x+p_y).^2./2+abs(p_z).^2+...
    abs((1i).*p_x-p_y).^2);
c_sca_m = (1./(k.^2.*(c^2*C_0.^2))).*(abs((1i).*m_x+m_y).^2./2+abs(m_z).^2+...
    abs((1i).*m_x-m_y).^2);
c_sca_qE = (1./(3.*k.^2.*(2*D_0.^2))).*(abs(2.*qE_xx+2.*(1i).*qE_xy+qE_zz).^2./4+...
    abs(qE_xz+(1i).*qE_yz).^2+abs(qE_xz-(1i).*qE_yz).^2+...
    abs(2.*(1i).*qE_xx+2.*qE_xy+qE_zz).^2./4+2/3.*abs(qE_zz).^2);
c_sca_qM = (1./(3.*k.^2.*(2*c.^2.*D_0.^2))).*(abs(2.*qM_xx+2.*(1i).*qM_xy+qM_zz).^2./4+...
    abs(qM_xz+(1i).*qM_yz).^2+abs((1i).*qM_yz-qM_xz).^2+...
    abs(2.*qM_xx-2.*(1i).*qM_xy+qM_zz).^2./4+2/3.*abs(qM_zz).^2);
c_sca_multipole = c_sca_p + c_sca_m + c_sca_qE + c_sca_qM;

% csca_tot = (k.^4./(6*pi*(eps_0)^2*(E_inc)^2))...
%     .*(C_0.*(p.^2)+((c.*C_0.*m.^2)./c));

figure()
plot(lamda2,abs(c_sca_p),'LineWidth',2)
hold on
plot(lamda2,abs(c_sca_m),'LineWidth',2)
hold on
plot(lamda2,abs(c_sca_qE),'LineWidth',2)
hold on
plot(lamda2,abs(c_sca_qM),'LineWidth',2)
hold on
plot(lamda2,abs(c_sca_multipole),'LineWidth',2)
hold on
plot(lamda2,abs(c_sca_num),'LineWidth',2)
ylabel('C_{sca} (m^2)')
xlabel('\lambda (nm)')
legend('c_{sca} p','c_{sca} m','c_{sca} qE_{xy}','c_{sca} qM_{xy}','c_{sca} multipole','c_{sca} SIE')

% figure()
% plot(lamda2,abs(c_sca_dipol),'LineWidth',2)
% hold on
% plot(lamda2,abs(c_sca_kuadrupol),'LineWidth',2)
% ylabel('C_{sca} (m^2)')
% xlabel('\lambda (nm)')
% legend('c_{sca} dipol','c_{sca} kuadrupol')

% figure()
% plot(lamda2,abs(p_x),'LineWidth',2)
% ylabel('p_x')
% figure()
% plot(lamda2,abs(p_y),'LineWidth',2)
% ylabel('p_y')
% figure()
% plot(lamda2,abs(p_z),'LineWidth',2)
% ylabel('p_z')
% figure()
% plot(lamda2,abs(m_x),'LineWidth',2)
% ylabel('m_x')
% figure()
% plot(lamda2,abs(m_y),'LineWidth',2)
% ylabel('m_y')
% figure()
% plot(lamda2,abs(m_z),'LineWidth',2)
% ylabel('m_z')
xlabel('\lambda (nm)')