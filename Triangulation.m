%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Program Triangulasi - SIE Method 														%
%% Source code ini mencari parameter2 semua elemen mesh segitiga pada permukaan scatterer   %
%% Nanda Perdana - 2016																		%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = p.'.*10^-9;
t = t.';

%% Jumlah segitiga, vector sisi segitiga, luas segitiga, keliling segitiga, dan vektor normal face segitiga
jumSegitiga = length(t);
for i = 1:jumSegitiga
    N = t(:,i);
    vec_sisi1 = p(:,N(3))-p(:,N(2));
    vec_sisi2 = p(:,N(1))-p(:,N(3));
    vec_sisi3 = p(:,N(2))-p(:,N(1));
    luas(i) = norm(cross(vec_sisi1,vec_sisi2))/2;
    keliling(i) = norm(vec_sisi1)+norm(vec_sisi2)+norm(vec_sisi3);
    faceNorm(i,:) = cross(vec_sisi1,vec_sisi2)./norm(cross(vec_sisi1,vec_sisi2));
end

%% titik tengah
for a = 1 : jumSegitiga
    tri = t(:,a);
    point1_x = p(1,tri(1));
    point1_y = p(2,tri(1));
    point1_z = p(3,tri(1));
    point2_x = p(1,tri(2));
    point2_y = p(2,tri(2));
    point2_z = p(3,tri(2));
    point3_x = p(1,tri(3));
    point3_y = p(2,tri(3));
    point3_z = p(3,tri(3));
    
    x = point1_x*1/3+point2_x*1/3+point3_x*1/3;
    y = point1_y*1/3+point2_y*1/3+point3_y*1/3;
    z = point1_z*1/3+point2_z*1/3+point3_z*1/3;
    titik_tengah(:,a) = [x y z];
end

%% Definisi indeks sisi, segi3 + dan -
sisi = [];
n = 0;
for i = 1:jumSegitiga
    N = t(1:3,i);
    for j = i+1:jumSegitiga
        M = t(1:3,j);
        a=1-all([N-M(1) N-M(2) N-M(3)]);
        if(sum(a)==2)
            n=n+1;
            sisi=[sisi M(find(a))];
            SegitigaPlus(n)=i;
            SegitigaMinus(n)=j;
        end
    end
end
jumSisi = length(sisi);

for i = 1:jumSisi
    PanjangSisi(i) = norm(p(:,sisi(1,i))-p(:,sisi(2,i)));
end

%% rho+- dan freevertex +-

for i = 1 : jumSisi
    plus = SegitigaPlus(i);
    plus1 = t(1,plus);
    plus2 = t(2,plus);
    plus3 = t(3,plus);
    if ((plus1~=sisi(1,i))&&(plus1~=sisi(2,i)))
        titik_sisi = plus1;
    end
    if ((plus2~=sisi(1,i))&&(plus2~=sisi(2,i)))
        titik_sisi = plus2;
    end
    if ((plus3~=sisi(1,i))&&(plus3~=sisi(2,i)))
        titik_sisi = plus3;
    end
    FreeVertex_plus(:,i) = p(:,titik_sisi);
    rho_plus(:,i) = titik_tengah(:,plus)-FreeVertex_plus(:,i);
end

for i = 1 : jumSisi
    mn = SegitigaMinus(i);
    minus1 = t(1,mn);
    minus2 = t(2,mn);
    minus3 = t(3,mn);
    if ((minus1~=sisi(1,i))&&(minus1~=sisi(2,i)))
        titik_sisi = minus1;
    end
    if ((minus2~=sisi(1,i))&&(minus2~=sisi(2,i)))
        titik_sisi = minus2;
    end
    if ((minus3~=sisi(1,i))&&(minus3~=sisi(2,i)))
        titik_sisi = minus3;
    end
    FreeVertex_minus(:,i) = p(:,titik_sisi);
    rho_minus(:,i) = titik_tengah(:,mn)-FreeVertex_minus(:,i);
end

%% logical hubungan antar segi3

p = p.';
t = t.';

for i = 1 : jumSegitiga
    obs = t(i,:);
    for j = 1 : jumSegitiga
        temp = 0;
        source = t(j,:);
        jar_center(i,j) = norm(titik_tengah(:,i)-titik_tengah(:,j));
        if i == j
            logic_identical(i,j) = 1;
        end
        if jar_center(i,j) > (keliling(i)+keliling(j))/2
            logic_far(i,j) = 1;
        end
        if p(obs(1),:) == p(source(1),:)
            temp = temp + 1;
        elseif p(obs(1),:) == p(source(2),:)
            temp = temp + 1;
        elseif p(obs(1),:) == p(source(3),:)
            temp = temp + 1;
        end
        if p(obs(2),:) == p(source(1),:)
            temp = temp + 1;
        elseif p(obs(2),:) == p(source(2),:)
            temp = temp + 1;
        elseif p(obs(2),:) == p(source(3),:)
            temp = temp + 1;
        end
        if p(obs(3),:) == p(source(1),:)
            temp = temp + 1;
        elseif p(obs(3),:) == p(source(2),:)
            temp = temp + 1;
        elseif p(obs(3),:) == p(source(3),:)
            temp = temp + 1;
        end
        
        if temp == 2
            logic_adjacent(i,j) = 1;
        else
            logic_adjacent(i,j) = 0;
        end
        if temp == 1
            logic_touch(i,j) = 1;
        else
            logic_touch(i,j) = 0;
        end
        
    end
end

%% Transpose correction

PanjangSisi = PanjangSisi.';
SegitigaPlus = SegitigaPlus.';
SegitigaMinus = SegitigaMinus.';
luas = luas.';
sisi = sisi.';
titik_tengah = titik_tengah.';
rho_plus = rho_plus.';
rho_minus = rho_minus.';
FreeVertex_plus = FreeVertex_plus.';
FreeVertex_minus = FreeVertex_minus.';
keliling = keliling.';

%% file saving

save parameter_segitiga p ...
    t ...
    PanjangSisi ...
    SegitigaMinus ...
    SegitigaPlus ...
    jumSegitiga ...
    jumSisi ...
    titik_tengah...
    rho_plus...
    rho_minus...
    FreeVertex_plus...
    FreeVertex_minus...
    luas ...
    sisi ...
    keliling ...
    jar_center ...
    faceNorm ...
    logic_touch ...
    logic_identical ...
    logic_far ...
    logic_adjacent
   
clearvars -except tStart

%% Mesh far-field dan parameter2nya 
%% anggap mesh kotak / meshgrid
 
load sph_farfield_20micro
p = p.';
t = t.';
jumSegitiga_ff = length(t);

for i = 1 : jumSegitiga_ff
    N = t(:,i);
    vec_sisi1 = p(:,N(3))-p(:,N(2));
    vec_sisi2 = p(:,N(1))-p(:,N(3));
    vec_sisi3 = p(:,N(2))-p(:,N(1));
    luas_ff(i,:) = norm(cross(vec_sisi1,vec_sisi2))/2;
    faceNorm_ff(i,:) = cross(vec_sisi3,-vec_sisi2)./norm(cross(vec_sisi3,-vec_sisi2));
    titik_tengah_ff(i,:)=sum(p(:,N),2)/3;
end

theta = 0:0.025:2*pi;
r = 20*10^-6;
clear x...
    y...
    z
x = r.*cos(theta);
y = r.*sin(theta);
titik_obs_2d_xz = [x; zeros(1,length(x)); y];
titik_obs_2d_yz = [zeros(1,length(x)); x; y];

%% File saving
save parameter_segitiga_farfield luas_ff...
    faceNorm_ff...
    titik_tengah_ff...
    jumSegitiga_ff

save parameter_segitiga_farfield_2d titik_obs_2d_xz...
    titik_obs_2d_yz
