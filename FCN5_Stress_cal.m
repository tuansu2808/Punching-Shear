function [stress_tem,downstand_sum] = FCN5_Stress_cal(hcx, hcy,...
    downstand_sum, Perimeter_num, Rz, ex, ey,...
    w_x, w_y, distributed_load, inner_load,...
    inner_x, inner_y, Col_position, d_eff,...
    Design_column_or_wall, drop_or_slab)


% downstand_sum = Drop_sum;
% Perimeter_num = Drop_perimeter_num;
% w_x = 0;
% w_y = 0;
% distributed_load = distributed_load_in_drop;
% inner_load = zeros(5,1);
% inner_x = hcx/1000;
% inner_y = hcy/1000;
% d_eff = d_drop_eff;
% drop_or_slab = 0;


% % % 
% hcx = Drop_x;
% hcy = Drop_y;
% downstand_sum = Slab_sum;
% Perimeter_num = Slab_perimeter_num;
% w_x = w_x_drop;
% w_y = w_y_drop;
% distributed_load = distributed_load_in_slab;
% inner_load = inner_force;
% inner_x = Drop_x_original/1000;
% inner_y = Drop_y_original/1000;

stress_tem = zeros(Perimeter_num,1);
downstand_sum(1).area_inner = 0;

for i0 = 1:Perimeter_num
    ele_tem = downstand_sum(i0).All_element;
    num_ele = length(ele_tem(:,3));
    area_section = sum(ele_tem(:,3).*(ones(num_ele,1)-ele_tem(:,4)).*ele_tem(:,5));

    downstand_sum(i0).total_length = sum(ele_tem(:,3).*(ones(num_ele,1)-ele_tem(:,4)));
    downstand_sum(i0).area = area_section;

    distance_to_edge = downstand_sum(i0).distance_to_edge;
    if str2double(Design_column_or_wall) == 1 && drop_or_slab == 1
        bx = hcx + 2*distance_to_edge;
    else
        bx = distance_to_edge + 1.5*d_eff + distance_to_edge/(3^0.5);
    end

    by = hcy + 2*distance_to_edge;
    
    num_load = length(ex);
    
    beta_tem = zeros(num_load,1);
    if Col_position == 1
        for i1 = 1:num_load
            beta_tem(i1,1) = 1 + 1.8*((ex(i1,1)/by)^2+(ey(i1,1)/bx)^2)^0.5;
        end
    else
        beta_tem = ones(num_load,1);
    end
    
    Rz_tem = Rz(:,1);

    area_section = downstand_sum(i0).area;
    area_inner = downstand_sum(i0).area_inner;

    if str2double(Design_column_or_wall) == 0 && drop_or_slab == 0
        Rz_tem = Rz(:,1);
    else
        if Col_position == 1
            for i10 = 1:5
                Rz_tem(i10,1) = max([Rz_tem(i10,1) - 2*(w_x+w_y)/1000 - ...
                    (area_inner/10^6)*distributed_load(i10,1) ... 
                    - inner_load(i10,1),0]);
            end
        elseif Col_position == 2
            for i10 = 1:5
                if i0 == 1
                    Rz_tem(i10,1) = max([Rz_tem(i10,1) - (w_x+2*w_y)/1000 - ...
                        inner_load(i10,1),0]);
                else
                    Rz_tem(i10,1) = max([Rz_tem(i10,1) - (w_x+2*w_y)/1000 - ...
                        (area_inner/10^6)*distributed_load(i10,1) ... 
                        - inner_load(i10,1),0]);
                end
            end
        elseif Col_position == 3
            for i10 = 1:5
                if i0 == 1
                    Rz_tem(i10,1) = max([Rz_tem(i10,1) - (w_x+w_y)/1000 - ...
                        inner_load(i10,1),0]);
                else
                    Rz_tem(i10,1) = max([Rz_tem(i10,1) - (w_x+w_y)/1000 - ...
                        (area_inner/10^6)*distributed_load(i10,1) ... 
                        - inner_load(i10,1),0]);
                end
            end
        end
    end
        
    stress_tem(i0,1) = 1000*max(beta_tem.*Rz_tem)/area_section;

    downstand_sum(i0).beta_tem = beta_tem;
    downstand_sum(i0).Rz_tem = Rz_tem;
    downstand_sum(i0).Max_force = max(beta_tem.*Rz_tem);
    downstand_sum(i0).bx = bx;
    downstand_sum(i0).by = by;
end