clearvars -except tStart
load parameter_segitiga

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Singularity Arcioni - SIE Method 												
% Source code ini mencari elemen matriks singular pada elemen segitiga kasus khusus
% (segitiga identik dan bertetangga satu sisi) dengan metode pada publikasi Arcioni    

% P. Arcioni, M. Bressan, L. Perregrini. "On the evaluation of the double surface 
% integrals arising in the application of the boundary integral method to 3-D problems". 
% IEEE Trans. Microw. Theory Tech., 45(3), 436-439 (1997)

%% Nanda Perdana - 2016																		
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_1_pm = zeros(jumSisi,jumSisi);
I_1_mp = zeros(jumSisi,jumSisi);
I_2_pm = zeros(jumSisi,jumSisi);
I_2_mp = zeros(jumSisi,jumSisi);

% Penentuan logika 2 segitiga yang identik

logic_identical_pp = logical(logic_identical(SegitigaPlus,SegitigaPlus));
logic_identical_pm = logical(logic_identical(SegitigaPlus,SegitigaMinus));
logic_identical_mp = logical(logic_identical(SegitigaMinus,SegitigaPlus));
logic_identical_mm = logical(logic_identical(SegitigaMinus,SegitigaMinus));

for a = 1 : jumSisi
    for b = 1 : jumSisi        
        
        % Plus plus
        if logic_identical_pp(a,b) 
            if a==b % identik
                vec_sisi1 = PanjangSisi(a);
                vec_sisi2 = FreeVertex_plus(a,:)-p(sisi(a,1),:);
                vec_sisi3 = FreeVertex_plus(a,:)-p(sisi(a,2),:);
                edge1 = vec_sisi1;
                edge2 = norm(vec_sisi2);
                edge3 = norm(vec_sisi3);
                half_kel = keliling(SegitigaPlus(a))/2;
                
                I_1_pp(a,b) = -4/3*luas(SegitigaPlus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(10+3*(edge3^2-edge1^2)/edge2^2-3*(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5-3*(edge1^2-edge2^2)/edge3^2-2*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5+3*(edge3^2-edge1^2)/edge2^2+2*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 2*log(1-edge1/half_kel)/edge1*(edge1^2-3*edge2^2-3*edge3^2-8*luas(SegitigaPlus(a))^2/edge1^2);
                part5_1 = 4*log(1-edge2/half_kel)/edge2*(edge1^2-2*edge2^2-4*edge3^2+6*luas(SegitigaPlus(a))^2/edge2^2);
                part6_1 = 4*log(1-edge3/half_kel)/edge3*(edge1^2-4*edge2^2-2*edge3^2+6*luas(SegitigaPlus(a))^2/edge3^2);
                
                I_2_pp(a,b) = (part1_1-part2_1-part3_1+part4_1+part5_1+part6_1)*luas(SegitigaPlus(a))^2/30;
            else % tetangga satu sisi
                edge1 = PanjangSisi(a,:);
                edge2 = PanjangSisi(b,:);
                edge3 = keliling(SegitigaPlus(a))-edge1-edge2;
                temp = edge1;
                edge1 = edge3;
                edge3 = temp;
                half_kel = keliling(SegitigaPlus(a))/2;
                
                I_1_pp(a,b) = -4/3*luas(SegitigaPlus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(-10+(edge3^2-edge1^2)/edge2^2-(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5+(edge1^2-edge2^2)/edge3^2-6*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5-(edge3^2-edge1^2)/edge2^2+6*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 12*log(1-edge1/half_kel)/edge1*(2*edge1^2-edge2^2-edge3^2+4*luas(SegitigaPlus(a))^2/edge1^2);
                part5_1 = 2*log(1-edge2/half_kel)/edge2*(9*edge1^2-3*edge2^2-edge3^2+4*luas(SegitigaPlus(a))^2/edge2^2);
                part6_1 = 2*log(1-edge3/half_kel)/edge3*(9*edge1^2-edge2^2-3*edge3^2+4*luas(SegitigaPlus(a))^2/edge3^2);
                
                I_2_pp(a,b) = (part1_1+part2_1+part3_1+part4_1+part5_1+part6_1)*luas(SegitigaPlus(a))^2/60;          
            end
        end
        
        % Plus minus
        if logic_identical_pm(a,b)
            if a==b% identik
                vec_sisi1 = PanjangSisi(a);
                vec_sisi2 = FreeVertex_plus(a,:)-p(sisi(a,1),:);
                vec_sisi3 = FreeVertex_plus(a,:)-p(sisi(a,2),:);
                edge1 = vec_sisi1;
                edge2 = norm(vec_sisi2);
                edge3 = norm(vec_sisi3);
                half_kel = keliling(SegitigaPlus(a))/2;
                
                I_1_pm(a,b) = -4/3*luas(SegitigaPlus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(10+3*(edge3^2-edge1^2)/edge2^2-3*(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5-3*(edge1^2-edge2^2)/edge3^2-2*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5+3*(edge3^2-edge1^2)/edge2^2+2*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 2*log(1-edge1/half_kel)/edge1*(edge1^2-3*edge2^2-3*edge3^2-8*luas(SegitigaPlus(a))^2/edge1^2);
                part5_1 = 4*log(1-edge2/half_kel)/edge2*(edge1^2-2*edge2^2-4*edge3^2+6*luas(SegitigaPlus(a))^2/edge2^2);
                part6_1 = 4*log(1-edge3/half_kel)/edge3*(edge1^2-4*edge2^2-2*edge3^2+6*luas(SegitigaPlus(a))^2/edge3^2);
                
                I_2_pm(a,b) = (part1_1-part2_1-part3_1+part4_1+part5_1+part6_1)*luas(SegitigaPlus(a))^2/30;
                
            else % tetangga satu sisi
                edge1 = PanjangSisi(a,:);
                edge2 = PanjangSisi(b,:);
                edge3 = keliling(SegitigaPlus(a))-edge1-edge2;
                temp = edge1;
                edge1 = edge3;
                edge3 = temp;
                half_kel = keliling(SegitigaPlus(a))/2;
                
                I_1_pm(a,b) = -4/3*luas(SegitigaPlus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(-10+(edge3^2-edge1^2)/edge2^2-(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5+(edge1^2-edge2^2)/edge3^2-6*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5-(edge3^2-edge1^2)/edge2^2+6*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 12*log(1-edge1/half_kel)/edge1*(2*edge1^2-edge2^2-edge3^2+4*luas(SegitigaPlus(a))^2/edge1^2);
                part5_1 = 2*log(1-edge2/half_kel)/edge2*(9*edge1^2-3*edge2^2-edge3^2+4*luas(SegitigaPlus(a))^2/edge2^2);
                part6_1 = 2*log(1-edge3/half_kel)/edge3*(9*edge1^2-edge2^2-3*edge3^2+4*luas(SegitigaPlus(a))^2/edge3^2);
                
                I_2_pm(a,b) = (part1_1+part2_1+part3_1+part4_1+part5_1+part6_1)*luas(SegitigaPlus(a))^2/60;          
            end
        end
        
        % Minus plus
        if logic_identical_mp(a,b)
            if a==b % identik
                vec_sisi1 = PanjangSisi(a);
                vec_sisi2 = FreeVertex_minus(a,:)-p(sisi(a,1),:);
                vec_sisi3 = FreeVertex_minus(a,:)-p(sisi(a,2),:);
                edge1 = vec_sisi1;
                edge2 = norm(vec_sisi2);
                edge3 = norm(vec_sisi3);
                half_kel = keliling(SegitigaMinus(a))/2;
                
                I_1_mp(a,b) = -4/3*luas(SegitigaMinus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(10+3*(edge3^2-edge1^2)/edge2^2-3*(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5-3*(edge1^2-edge2^2)/edge3^2-2*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5+3*(edge3^2-edge1^2)/edge2^2+2*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 2*log(1-edge1/half_kel)/edge1*(edge1^2-3*edge2^2-3*edge3^2-8*luas(SegitigaMinus(a))^2/edge1^2);
                part5_1 = 4*log(1-edge2/half_kel)/edge2*(edge1^2-2*edge2^2-4*edge3^2+6*luas(SegitigaMinus(a))^2/edge2^2);
                part6_1 = 4*log(1-edge3/half_kel)/edge3*(edge1^2-4*edge2^2-2*edge3^2+6*luas(SegitigaMinus(a))^2/edge3^2);
                
                I_2_mp(a,b) = (part1_1-part2_1-part3_1+part4_1+part5_1+part6_1)*luas(SegitigaMinus(a))^2/30;
            else % tetangga satu sisi
                edge1 = PanjangSisi(a,:);
                edge2 = PanjangSisi(b,:);
                edge3 = keliling(SegitigaMinus(a))-edge1-edge2;
                temp = edge1;
                edge1 = edge3;
                edge3 = temp;
                half_kel = keliling(SegitigaMinus(a))/2;
                
                I_1_mp(a,b) = -4/3*luas(SegitigaMinus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(-10+(edge3^2-edge1^2)/edge2^2-(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5+(edge1^2-edge2^2)/edge3^2-6*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5-(edge3^2-edge1^2)/edge2^2+6*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 12*log(1-edge1/half_kel)/edge1*(2*edge1^2-edge2^2-edge3^2+4*luas(SegitigaMinus(a))^2/edge1^2);
                part5_1 = 2*log(1-edge2/half_kel)/edge2*(9*edge1^2-3*edge2^2-edge3^2+4*luas(SegitigaMinus(a))^2/edge2^2);
                part6_1 = 2*log(1-edge3/half_kel)/edge3*(9*edge1^2-edge2^2-3*edge3^2+4*luas(SegitigaMinus(a))^2/edge3^2);
                
                I_2_mp(a,b) = (part1_1+part2_1+part3_1+part4_1+part5_1+part6_1)*luas(SegitigaMinus(a))^2/60;          
            end
        end
        
        % Minus minus
        if logic_identical_mm(a,b)
            if a==b % identik
                vec_sisi1 = PanjangSisi(a);
                vec_sisi2 = FreeVertex_minus(a,:)-p(sisi(a,1),:);
                vec_sisi3 = FreeVertex_minus(a,:)-p(sisi(a,2),:);
                edge1 = vec_sisi1;
                edge2 = norm(vec_sisi2);
                edge3 = norm(vec_sisi3);
                half_kel = keliling(SegitigaMinus(a))/2;
                
                I_1_mm(a,b) = -4/3*luas(SegitigaMinus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(10+3*(edge3^2-edge1^2)/edge2^2-3*(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5-3*(edge1^2-edge2^2)/edge3^2-2*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5+3*(edge3^2-edge1^2)/edge2^2+2*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 2*log(1-edge1/half_kel)/edge1*(edge1^2-3*edge2^2-3*edge3^2-8*luas(SegitigaMinus(a))^2/edge1^2);
                part5_1 = 4*log(1-edge2/half_kel)/edge2*(edge1^2-2*edge2^2-4*edge3^2+6*luas(SegitigaMinus(a))^2/edge2^2);
                part6_1 = 4*log(1-edge3/half_kel)/edge3*(edge1^2-4*edge2^2-2*edge3^2+6*luas(SegitigaMinus(a))^2/edge3^2);
                
                I_2_mm(a,b) = (part1_1-part2_1-part3_1+part4_1+part5_1+part6_1)*luas(SegitigaMinus(a))^2/30;
            else % tetangga satu sisi
                edge1 = PanjangSisi(a,:);
                edge2 = PanjangSisi(b,:);
                edge3 = keliling(SegitigaMinus(a))-edge1-edge2;
                temp = edge1;
                edge1 = edge3;
                edge3 = temp;
                half_kel = keliling(SegitigaMinus(a))/2;
                
                I_1_mm(a,b) = -4/3*luas(SegitigaMinus(a)).^2.*(log(1-edge1/half_kel)/edge1+...
                    log(1-edge2/half_kel)/edge2+...
                    log(1-edge3/half_kel)/edge3);
                
                part1_1 = edge1*(-10+(edge3^2-edge1^2)/edge2^2-(edge1^2-edge2^2)/edge3^2);
                part2_1 = edge2*(5+(edge1^2-edge2^2)/edge3^2-6*(edge2^2-edge3^2)/edge1^2);
                part3_1 = edge3*(5-(edge3^2-edge1^2)/edge2^2+6*(edge2^2-edge3^2)/edge1^2);
                part4_1 = 12*log(1-edge1/half_kel)/edge1*(2*edge1^2-edge2^2-edge3^2+4*luas(SegitigaMinus(a))^2/edge1^2);
                part5_1 = 2*log(1-edge2/half_kel)/edge2*(9*edge1^2-3*edge2^2-edge3^2+4*luas(SegitigaMinus(a))^2/edge2^2);
                part6_1 = 2*log(1-edge3/half_kel)/edge3*(9*edge1^2-edge2^2-3*edge3^2+4*luas(SegitigaMinus(a))^2/edge3^2);
                
                I_2_mm(a,b) = (part1_1+part2_1+part3_1+part4_1+part5_1+part6_1)*luas(SegitigaMinus(a))^2/60;          
            end
        end
    end
end

% save file

save singular_arcioni_tp I_1_pp...
    I_1_pm...
    I_1_mp...
    I_1_mm...
    I_2_pp...
    I_2_pm...
    I_2_mp...
    I_2_mm...