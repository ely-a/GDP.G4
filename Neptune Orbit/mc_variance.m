
count = 1;
percentiles = [];
for i = 10:N_MC
    percentiles(count) = prctile(max_m_TPS(1:i),99);
    count = count+1;
end

figure
plot(10:N_MC, percentiles)