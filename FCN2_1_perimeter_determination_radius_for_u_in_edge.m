function Downstand_sum = FCN2_1_perimeter_determination_radius_for_u_in_edge...
    (d_effective, hcx, hcy, ...
    Hole, Perimeter_num, number_of_hole, jump_step, Col_position)

% Gross length of each each drop perimeters
Downstand_sum = struct;
ele_size = 10;
ele_size_angular = 0.01;

for i1 = 1:Perimeter_num

    if i1 == 1
        distance_to_edge = 0;
    else
        distance_to_edge = 2*d_effective + (i1-2)*jump_step*d_effective;
    end

    Downstand_sum(i1).distance_to_edge = distance_to_edge;
    Downstand_sum(i1).area_inside = hcx*hcy+...
        2*distance_to_edge*(hcx+hcy)+...
        pi*distance_to_edge*distance_to_edge;


    Downstand_sum(i1).L_gross(1).x(1) = -hcx/2;
    Downstand_sum(i1).L_gross(1).x(2) = hcx/2;
    Downstand_sum(i1).L_gross(1).y(1) = hcy/2+distance_to_edge;
    Downstand_sum(i1).L_gross(1).y(2) = hcy/2+distance_to_edge;
    
    Downstand_sum(i1).L_gross(1).x_mid = -hcx/2+ele_size/2:ele_size:hcx/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(1).x_mid);
    Downstand_sum(i1).L_gross(1).y_mid = (hcy/2+distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(1).ele_size = ele_size;
    

    Downstand_sum(i1).L_gross(2).x(1) = hcx/2+distance_to_edge;
    Downstand_sum(i1).L_gross(2).x(2) = hcx/2+distance_to_edge;
    Downstand_sum(i1).L_gross(2).y(1) = hcy/2;
    Downstand_sum(i1).L_gross(2).y(2) = -hcy/2;

    Downstand_sum(i1).L_gross(2).y_mid = -hcy/2+ele_size/2:ele_size:hcy/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(2).y_mid);
    Downstand_sum(i1).L_gross(2).x_mid = (hcx/2+distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(2).ele_size = ele_size;


    Downstand_sum(i1).L_gross(3).x(1) = -hcx/2;
    Downstand_sum(i1).L_gross(3).x(2) = hcx/2;
    Downstand_sum(i1).L_gross(3).y(1) = -hcy/2-distance_to_edge;
    Downstand_sum(i1).L_gross(3).y(2) = -hcy/2-distance_to_edge;

    Downstand_sum(i1).L_gross(3).x_mid = -hcx/2+ele_size/2:ele_size:hcx/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(3).x_mid);
    Downstand_sum(i1).L_gross(3).y_mid = (-hcy/2-distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(3).ele_size = ele_size;


    Downstand_sum(i1).L_gross(4).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum(i1).L_gross(4).x(2) = -hcx/2-distance_to_edge;
    Downstand_sum(i1).L_gross(4).y(1) = hcy/2;
    Downstand_sum(i1).L_gross(4).y(2) = -hcy/2;

    Downstand_sum(i1).L_gross(4).y_mid = -hcy/2+ele_size/2:ele_size:hcy/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross(4).y_mid);
    Downstand_sum(i1).L_gross(4).x_mid = (-hcx/2-distance_to_edge)*ones(1,num_ele);
    Downstand_sum(i1).L_gross(4).ele_size = ele_size;


    Downstand_sum(i1).C_gross(1).cx = -hcx/2;
    Downstand_sum(i1).C_gross(1).cy = hcy/2;
    Downstand_sum(i1).C_gross(1).r = distance_to_edge;
    Downstand_sum(i1).C_gross(1).angle(1) = pi/2;
    Downstand_sum(i1).C_gross(1).angle(2) = pi;
    Downstand_sum(i1).C_gross(1).ang = pi/2+ele_size_angular/2:...
                            ele_size_angular:pi-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross(1).ang;
    Downstand_sum(i1).C_gross(1).x_mid = distance_to_edge*cos(ang)-hcx/2;
    Downstand_sum(i1).C_gross(1).y_mid = distance_to_edge*sin(ang)+hcy/2;
    Downstand_sum(i1).C_gross(1).ele_size = distance_to_edge*ele_size_angular;


    Downstand_sum(i1).C_gross(2).cx = hcx/2;
    Downstand_sum(i1).C_gross(2).cy = hcy/2;
    Downstand_sum(i1).C_gross(2).r = distance_to_edge;
    Downstand_sum(i1).C_gross(2).angle(1) = 0;
    Downstand_sum(i1).C_gross(2).angle(2) = pi/2;
    Downstand_sum(i1).C_gross(2).ang = 0+ele_size_angular/2:...
                        ele_size_angular:pi/2-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross(2).ang;
    Downstand_sum(i1).C_gross(2).x_mid = distance_to_edge*cos(ang)+hcx/2;
    Downstand_sum(i1).C_gross(2).y_mid = distance_to_edge*sin(ang)+hcy/2;
    Downstand_sum(i1).C_gross(2).ele_size = distance_to_edge*ele_size_angular;


    Downstand_sum(i1).C_gross(3).cx = hcx/2;
    Downstand_sum(i1).C_gross(3).cy = -hcy/2;
    Downstand_sum(i1).C_gross(3).r = distance_to_edge;
    Downstand_sum(i1).C_gross(3).angle(1) = 3/2*pi;
    Downstand_sum(i1).C_gross(3).angle(2) = 2*pi;
    Downstand_sum(i1).C_gross(3).ang = 3/2*pi+ele_size_angular/2:...
                            ele_size_angular:2*pi-ele_size_angular/2;
    ang = Downstand_sum(i1).C_gross(3).ang;
    Downstand_sum(i1).C_gross(3).x_mid = distance_to_edge*cos(ang)+hcx/2;
    Downstand_sum(i1).C_gross(3).y_mid = distance_to_edge*sin(ang)-hcy/2;
    Downstand_sum(i1).C_gross(3).ele_size = distance_to_edge*ele_size_angular;


    Downstand_sum(i1).C_gross(4).cx = -hcx/2;
    Downstand_sum(i1).C_gross(4).cy = -hcy/2;
    Downstand_sum(i1).C_gross(4).r = distance_to_edge;
    Downstand_sum(i1).C_gross(4).angle(1) = pi;
    Downstand_sum(i1).C_gross(4).angle(2) = 3/2*pi;
    Downstand_sum(i1).C_gross(4).ang = pi+ele_size_angular/2:...
                        ele_size_angular:3/2*pi-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross(4).ang;
    Downstand_sum(i1).C_gross(4).x_mid = distance_to_edge*cos(ang)-hcx/2;
    Downstand_sum(i1).C_gross(4).y_mid = distance_to_edge*sin(ang)-hcy/2;
    Downstand_sum(i1).C_gross(4).ele_size = distance_to_edge*ele_size_angular;













    if i1 == 1
        distance_to_edge_2 = 0;
    else
        distance_to_edge_2 = distance_to_edge - 1.5*d_effective;
    end

    Downstand_sum(i1).distance_to_edge_inner = distance_to_edge_2;


    Downstand_sum(i1).L_gross_inner(1).x(1) = -hcx/2;
    Downstand_sum(i1).L_gross_inner(1).x(2) = hcx/2;
    Downstand_sum(i1).L_gross_inner(1).y(1) = hcy/2+distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(1).y(2) = hcy/2+distance_to_edge_2;
    
    Downstand_sum(i1).L_gross_inner(1).x_mid = -hcx/2+ele_size/2:ele_size:hcx/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross_inner(1).x_mid);
    Downstand_sum(i1).L_gross_inner(1).y_mid = (hcy/2+distance_to_edge_2)*ones(1,num_ele);
    Downstand_sum(i1).L_gross_inner(1).ele_size = ele_size;
    

    Downstand_sum(i1).L_gross_inner(2).x(1) = hcx/2+distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(2).x(2) = hcx/2+distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(2).y(1) = hcy/2;
    Downstand_sum(i1).L_gross_inner(2).y(2) = -hcy/2;

    Downstand_sum(i1).L_gross_inner(2).y_mid = -hcy/2+ele_size/2:ele_size:hcy/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross_inner(2).y_mid);
    Downstand_sum(i1).L_gross_inner(2).x_mid = (hcx/2+distance_to_edge_2)*ones(1,num_ele);
    Downstand_sum(i1).L_gross_inner(2).ele_size = ele_size;


    Downstand_sum(i1).L_gross_inner(3).x(1) = -hcx/2;
    Downstand_sum(i1).L_gross_inner(3).x(2) = hcx/2;
    Downstand_sum(i1).L_gross_inner(3).y(1) = -hcy/2-distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(3).y(2) = -hcy/2-distance_to_edge_2;

    Downstand_sum(i1).L_gross_inner(3).x_mid = -hcx/2+ele_size/2:ele_size:hcx/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross_inner(3).x_mid);
    Downstand_sum(i1).L_gross_inner(3).y_mid = (-hcy/2-distance_to_edge_2)*ones(1,num_ele);
    Downstand_sum(i1).L_gross_inner(3).ele_size = ele_size;


    Downstand_sum(i1).L_gross_inner(4).x(1) = -hcx/2-distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(4).x(2) = -hcx/2-distance_to_edge_2;
    Downstand_sum(i1).L_gross_inner(4).y(1) = hcy/2;
    Downstand_sum(i1).L_gross_inner(4).y(2) = -hcy/2;

    Downstand_sum(i1).L_gross_inner(4).y_mid = -hcy/2+ele_size/2:ele_size:hcy/2-ele_size/2;
    num_ele = length(Downstand_sum(i1).L_gross_inner(4).y_mid);
    Downstand_sum(i1).L_gross_inner(4).x_mid = (-hcx/2-distance_to_edge_2)*ones(1,num_ele);
    Downstand_sum(i1).L_gross_inner(4).ele_size = ele_size;


    Downstand_sum(i1).C_gross_inner(1).cx = -hcx/2;
    Downstand_sum(i1).C_gross_inner(1).cy = hcy/2;
    Downstand_sum(i1).C_gross_inner(1).r = distance_to_edge_2;
    Downstand_sum(i1).C_gross_inner(1).angle(1) = pi/2;
    Downstand_sum(i1).C_gross_inner(1).angle(2) = pi;
    Downstand_sum(i1).C_gross_inner(1).ang = pi/2+ele_size_angular/2:...
                            ele_size_angular:pi-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross_inner(1).ang;
    Downstand_sum(i1).C_gross_inner(1).x_mid = distance_to_edge_2*cos(ang)-hcx/2;
    Downstand_sum(i1).C_gross_inner(1).y_mid = distance_to_edge_2*sin(ang)+hcy/2;
    Downstand_sum(i1).C_gross_inner(1).ele_size = distance_to_edge_2*ele_size_angular;


    Downstand_sum(i1).C_gross_inner(2).cx = hcx/2;
    Downstand_sum(i1).C_gross_inner(2).cy = hcy/2;
    Downstand_sum(i1).C_gross_inner(2).r = distance_to_edge_2;
    Downstand_sum(i1).C_gross_inner(2).angle(1) = 0;
    Downstand_sum(i1).C_gross_inner(2).angle(2) = pi/2;
    Downstand_sum(i1).C_gross_inner(2).ang = 0+ele_size_angular/2:...
                        ele_size_angular:pi/2-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross_inner(2).ang;
    Downstand_sum(i1).C_gross_inner(2).x_mid = distance_to_edge_2*cos(ang)+hcx/2;
    Downstand_sum(i1).C_gross_inner(2).y_mid = distance_to_edge_2*sin(ang)+hcy/2;
    Downstand_sum(i1).C_gross_inner(2).ele_size = distance_to_edge_2*ele_size_angular;


    Downstand_sum(i1).C_gross_inner(3).cx = hcx/2;
    Downstand_sum(i1).C_gross_inner(3).cy = -hcy/2;
    Downstand_sum(i1).C_gross_inner(3).r = distance_to_edge_2;
    Downstand_sum(i1).C_gross_inner(3).angle(1) = 3/2*pi;
    Downstand_sum(i1).C_gross_inner(3).angle(2) = 2*pi;
    Downstand_sum(i1).C_gross_inner(3).ang = 3/2*pi+ele_size_angular/2:...
                            ele_size_angular:2*pi-ele_size_angular/2;
    ang = Downstand_sum(i1).C_gross_inner(3).ang;
    Downstand_sum(i1).C_gross_inner(3).x_mid = distance_to_edge_2*cos(ang)+hcx/2;
    Downstand_sum(i1).C_gross_inner(3).y_mid = distance_to_edge_2*sin(ang)-hcy/2;
    Downstand_sum(i1).C_gross_inner(3).ele_size = distance_to_edge_2*ele_size_angular;


    Downstand_sum(i1).C_gross_inner(4).cx = -hcx/2;
    Downstand_sum(i1).C_gross_inner(4).cy = -hcy/2;
    Downstand_sum(i1).C_gross_inner(4).r = distance_to_edge_2;
    Downstand_sum(i1).C_gross_inner(4).angle(1) = pi;
    Downstand_sum(i1).C_gross_inner(4).angle(2) = 3/2*pi;
    Downstand_sum(i1).C_gross_inner(4).ang = pi+ele_size_angular/2:...
                        ele_size_angular:3/2*pi-ele_size_angular/2; 
    ang = Downstand_sum(i1).C_gross_inner(4).ang;
    Downstand_sum(i1).C_gross_inner(4).x_mid = distance_to_edge_2*cos(ang)-hcx/2;
    Downstand_sum(i1).C_gross_inner(4).y_mid = distance_to_edge_2*sin(ang)-hcy/2;
    Downstand_sum(i1).C_gross_inner(4).ele_size = distance_to_edge_2*ele_size_angular;

end

% Intersection

x_2 = -hcx/2;
y_2 = -hcy/2;

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
                if Downstand_sum(i1).L_gross(i2).x_mid(i3) <=  x_2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end
            elseif Col_position == 3
                if Downstand_sum(i1).L_gross(i2).x_mid(i3) <=  x_2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end

                if Downstand_sum(i1).L_gross(i2).y_mid(i3) <=  y_2
                    Downstand_sum(i1).L_gross(i2).intersect(i3) = 1;
                end
            end

        end
    end

    for i2 = 1:4
        num_ele = length(Downstand_sum(i1).C_gross(i2).x_mid);
        intersec_tem = zeros(num_ele,number_of_hole);
        element_all = [Downstand_sum(i1).C_gross(i2).x_mid',...
            Downstand_sum(i1).C_gross(i2).y_mid'];

        for i3 = 1:number_of_hole
            shape_tem = Hole(i3).shape;
            in_tem = inpolygon(element_all(:,1),element_all(:,2)...
                ,shape_tem.Vertices(:,1),shape_tem.Vertices(:,2));
            intersec_tem(:,i3) = in_tem;
        end
        
        for i3 = 1:num_ele
            inter_sum = sum(intersec_tem(i3,:));
            if inter_sum == 0
                Downstand_sum(i1).C_gross(i2).intersect(i3) = 0;
            else
                Downstand_sum(i1).C_gross(i2).intersect(i3) = 1;
            end

            if Col_position == 2
                if Downstand_sum(i1).C_gross(i2).x_mid(i3) <=  x_2
                    Downstand_sum(i1).C_gross(i2).intersect(i3) = 1;
                end
            elseif Col_position == 3
                if Downstand_sum(i1).C_gross(i2).x_mid(i3) <=  x_2
                    Downstand_sum(i1).C_gross(i2).intersect(i3) = 1;
                end

                if Downstand_sum(i1).C_gross(i2).y_mid(i3) <=  y_2
                    Downstand_sum(i1).C_gross(i2).intersect(i3) = 1;
                end
            end
        end
    end
end
