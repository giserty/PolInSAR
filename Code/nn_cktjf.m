function [yd_lidar,yd_sg,yd_num ] = nn_cktjf(winh,winl,sg,lidar,bhs)
%将矩阵划分为winh行 winl列进行样地统计  有效的样地数>bhs才计入
yd_num=0;
t=size(sg);
h=t(1);
l=t(2);
for i=1:winh:h
    for j=1:winl:l
        num=0;
        sgnum=0;
        lidarnum=0;
        if i+winh-1<=h
            h1=i+winh-1;
        else
            h1=h;
        end
        if j+winl-1<=l
            l1=j+winl-1;
        else
            l1=l;
        end
        for k1=i:h1
            for k2=j:l1
                % if sgmask(k1,k2)&& sg(k1,k2)>0 &&lidar(k1,k2)>0
                if sg(k1,k2)>0 &&lidar(k1,k2)>0
                    num=num+1;
                   sgnum=sgnum+sg(k1,k2);
                   lidarnum=lidarnum+lidar(k1,k2);
                end
            end
        end
        if num/winh/winl>=bhs
            yd_num=yd_num+1;
            yd_sg(yd_num)=sgnum/num;
            yd_lidar(yd_num)=lidarnum/num;
        end
    end
end


end

