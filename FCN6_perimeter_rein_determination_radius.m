function Downstand_sum = FCN6_perimeter_rein_determination_radius...
    (d_effective, hcx, hcy, ...
    Hole, Perimeter_num, ...
    number_of_hole, Col_position,...
    hcx_col, hcy_col)


% Gross length of each each drop perimeters
Downstand_sum = struct;
ele_size = 50;

for i1 = 1:Perimeter_num

    distance_to_edge = 0.5*d_effective + (i1-1)*0.75*d_effective;

    Downstand_sum(i1).distance_to_edge = distance_to_edge;
    Downstand_sum(i1).area_inside = hcx*hcy+...
        2*distance_to_edge*(hcx+hcy)+...
        pi*distance_to_edge*distance_to_edge;

        Downstand_sum(i1).L_gross(1).x(1) = -hcx/2-distance_to_edge;
        Downstand_sum(i1).L_gross(1).x(2) = hcx/2+distance_to_edge;
        Downstand_sum(i1).L_gross(1).y(1) = hcy/2+distance_to_edge;
        Downstand_sum(i1).L_gross(1).y(2) = hcy/2+distance_to_edge;
        
        Downstand_sum(i1).L_gross(1).x_mid = -hcx/2-distance_to_edge+ele_size/2:...
                                        ele_size:hcx/2+distance_to_edge-ele_size/2;
        num_ele = length(Downstand_sum(i1).L_gross(1).x_mid);
        Downstand_sum(i1).L_gross(1).y_mid = ...
                                        (hcy/2+distance_to_edge)*ones(1,num_ele);
        Downstand_sum(i1).L_gross(1).ele_size = ele_size;
    



    Downstand_sum(i1).L_gross(2).x(1) = hcx/2+distance_to_edge;
    Downstand_sum(i1).L_gross(2).x(2) = hcx/2+distance_to_edge;
    Downstand_sum(i1).L_gross(2).y(1) = hcy/2+distance_to_edge;
    Downstand_sum(i1).L_gross(2).y(2) = -hcy/2-distance_to_edge;

    Downstand_sum(i1).L_gross(2).y_mid = -hcy/2-distance_to_edge+ele_size/2:...
                                    ele_size:hcy/2+distance_to_edge-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(2).y_mid);
    Downstand_sum(i1).L_gross(2).x_mid = ...
                                    (hcx/2+distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(2).ele_size = ele_size;




    Downstand_sum(i1).L_gross(3).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum(i1).L_gross(3).x(2) = hcx/2+distance_to_edge;
    Downstand_sum(i1).L_gross(3).y(1) = -hcy/2-distance_to_edge;
    Downstand_sum(i1).L_gross(3).y(2) = -hcy/2-distance_to_edge;

    Downstand_sum(i1).L_gross(3).x_mid = -hcx/2-distance_to_edge+ele_size/2:...
                                    ele_size:hcx/2+distance_to_edge-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(3).x_mid);
    Downstand_sum(i1).L_gross(3).y_mid = ...
                                    (-hcy/2-distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(3).ele_size = ele_size;




    Downstand_sum(i1).L_gross(4).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum(i1).L_gross(4).x(2) = -hcx/2-distance_to_edge;
    Downstand_sum(i1).L_gross(4).y(1) = hcy/2+distance_to_edge;
    Downstand_sum(i1).L_gross(4).y(2) = -hcy/2-distance_to_edge;

    Downstand_sum(i1).L_gross(4).y_mid = -hcy/2-distance_to_edge+ele_size/2:...
                                    ele_size:hcy/2+distance_to_edge-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(4).y_mid);
    Downstand_sum(i1).L_gross(4).x_mid = ...
                                    (-hcx/2-distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(4).ele_size = ele_size;

end

% Intersection

for i1 = 1:Perimeter_num
    for i2 = 1:4
        num_ele = length(Downstand_sum(i1).L_gross(i2).x_mid);
        intersec_tem = zeros(num_ele,number_of_hole);
        element_all = [Downstand_sum(i1).L_gross(i2).x_mid',...
            Downstand_sum(i1).L_gross(i2).y_mid'];

        for i3 = 1:number_of_hole
            shape_tem = Hole(i3).shape;
            in_tem = inpolygon(element_all(:,1),element_all(:,2)...
                ,shape_tem.Vertices(:,1),shape_tem.Vertices(:,2));
            intersec_tem(:,i3) = in_tem;
        end
        
        for i3 = 1:num_ele
            inter_sum = sum(intersec_tem(i3,:));
            if inter_sum == 0
                Downstand_sum(i1).L_gross(i2).intersect(i3) = 0;
            else
                Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
            end

            if Col_position == 2
                if Downstand_sum(i1).L_gross(i2).x_mid(i3) <=  -hcx_col/2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end
            elseif Col_position == 3
                if Downstand_sum(i1).L_gross(i2).x_mid(i3) <=  -hcx_col/2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end

                if Downstand_sum(i1).L_gross(i2).y_mid(i3) <=  -hcy_col/2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end
            end

        end
    end
end
