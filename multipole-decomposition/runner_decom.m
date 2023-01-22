clearvars -except avg_bench
close all
load csca_square_5nm.mat % csca from SIE

%%%%% input %%%%%
wl_interp = (350:2:800).'; % make sure that wl is consistent with the csca
N_multipole = 2;
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
figure()
plot(lamda.*1e9,c_sca_bohren_total,'-','LineWidth',2)
hold on
plot(lamda.*1e9,c_sca_num,'--','LineWidth',2)
for i = 1 : jum_n
    plot(lamda.*1e9,c_sca_cont_E(:,i),'color',col(i,:),'LineWidth',2)
    hold on
    plot(lamda.*1e9,c_sca_cont_M(:,i),'--','color',col(i,:),'LineWidth',2)
    hold on
end
xlabel('\lambda (nm)')
%xlim([min(lamda_interp) max(lamda_interp)])
legend('c_{sca} decomp','c_{sca} SIE','E_1','M_1','E_2','M_2')
ylim([0 1.5*10^(-21)])

%save res_aggregate_5_overlap_Ex.mat...
    %c_sca_bohren_total...
    %c_sca_cont_E...
    %c_sca_cont_M...
    %c_sca_num...
    %lamda