function  barsem(to_bar,err_bar,p_val)
[nbars,ngroups] = size(to_bar);
bar(to_bar');
hold all
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = [];
for i = 1:nbars
    x(i,:) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x(i,:), to_bar(i,:), err_bar(i,:), '.g');
end

sig_data = cell(1,ngroups);
for i=1:ngroups
    sig_data{i} = [x(1,i),x(2,i)];
end
sigstar(sig_data,p_val);

end

