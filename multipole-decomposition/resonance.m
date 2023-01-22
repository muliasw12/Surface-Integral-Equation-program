close all
load res_aggregate_5_overlap_Ex  %from multipole decomposition

csca_E = c_sca_cont_E(106,:); %sesuaikan dengan titik resonance
csca_M = c_sca_cont_M(106,:); %sesuaikan dengan titik resonance
csca_tot = c_sca_bohren_total(106,:); %sesuaikan dengan titik resonance
total = 0;

%ratio contribution
for i = 1 : 10 %make sure N multipole in multipole decomposition
    res_E = csca_E./csca_tot;
    res_M = csca_M./csca_tot;
end

res_E = res_E.';
res_M = res_M.';
res = res_E + res_M;
sum_res = sum(res);
total = cat(2, res_E, res_M, res);
%result
disp(res)
disp(sum_res)




    

