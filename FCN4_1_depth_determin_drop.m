function [Drop_sum] = FCN4_1_depth_determin_drop(Drop_sum,...
    S_thickness, Drop_thickness, Drop_x, Drop_y, Drop_perimeter_num)

x_drop_max = Drop_x/2;
x_drop_min = -Drop_x/2;
y_drop_max = Drop_y/2;
y_drop_min = -Drop_y/2;


peri_num = Drop_perimeter_num;
downstand_tem = Drop_sum;

for i1 = 1:peri_num
    downstand_tem(i1).All_element = [];

    for i2 = 1:4
        x_mid = downstand_tem(i1).L_gross(i2).x_mid;
        y_mid = downstand_tem(i1).L_gross(i2).y_mid;

        num_point = length(x_mid);
        
        for i3 = 1:num_point
            x_tem = x_mid(i3);
            y_tem = y_mid(i3);

            if x_tem < x_drop_max && x_tem > x_drop_min ...
                    && y_tem < y_drop_max && y_tem > y_drop_min
                downstand_tem(i1).L_gross(i2).depth(i3) = Drop_thickness;
            else
                downstand_tem(i1).L_gross(i2).depth(i3) = S_thickness;
            end
        end

        ele_tem = [x_mid', y_mid'...
            ,downstand_tem(i1).L_gross(i2).ele_size*ones(length(x_mid),1)...
            ,downstand_tem(i1).L_gross(i2).intersect' ...
            ,downstand_tem(i1).L_gross(i2).depth'];

        downstand_tem(i1).All_element = [downstand_tem(i1).All_element;ele_tem];
    end
    downstand_tem(i1).All_element = unique(downstand_tem(i1).All_element,"rows");
end

Drop_sum = downstand_tem;