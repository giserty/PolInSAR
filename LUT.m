function [hvmap,extmap] = LUT(fai0,num_hv,num_ext,ywv,m,n,theta,kz)

hv_vector = linspace(0.01,60.01,num_hv);
ext_vector = linspace(0.001,0.115,num_ext);

mindist = ones(m,n) .* 1e9;
mindist = mindist(:);
limit2piclip = ones(m,n);
limit2piclip = limit2piclip(:);
kzclip = kz(:);

hvfit = zeros(m,n);
hvfit = hvfit(:);
extfit = zeros(m,n);
extfit = extfit(:);


for i = 1:length(hv_vector)
%     disp(i);
    hv_val = hv_vector(i);
    for j = 1:length(ext_vector)
        ext_val = ext_vector(j);
        y_v = ywv ./ exp(fai0*1i);
        dist = abs(y_v-2.*ext_val.*(exp(2.*ext_val.*hv_val./cos(theta)+kz.*hv_val.*1i)-1)./((2.*ext_val+kz.*cos(theta).*1i).*(exp(2.*ext_val.*hv_val./cos(theta))-1)));
        dist = dist(:);
        ind_limit = limit2piclip & (hv_val > abs(2*pi./kzclip));
        if any(ind_limit(:))
            dist(ind_limit) = 1e10;
        end
        ind = dist < mindist;
        if any(ind(:))
            mindist(ind) = dist(ind);
            hvfit(ind) = hv_val;  % hv_val是数值  hvfit是列表
            extfit(ind) = ext_val;
        end
    end
end

hvmap = reshape(hvfit,m,n);
extmap = reshape(extfit,m,n);

end
