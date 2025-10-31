function Downstand_sum_rein = FCN8_Arrangement_of_shear...
    (Asw_rein,Downstand_sum_rein, d_eff, hcx, hcy, ...
    Hole, Perimeter_rein_num, number_of_hole, Area_link,...
    Downstand_sum, ~, Last_peri, outter_peri_density,...
    Drop_x, Drop_y, Col_position, hcx_col, hcy_col)

% Gross length of each each drop perimeters
stand_peri_num = length(Downstand_sum);

for i1 = 1:Perimeter_rein_num
    distance_to_edge = Downstand_sum_rein(i1).distance_to_edge;

    ele_size = Asw_rein(i1,3);

    distance_tip = min(d_eff,ele_size);

    all_length = 2*(hcx + hcy + 4*distance_to_edge);
    all_rein_num = ceil(all_length/ele_size);

    Num_x = max([ceil(all_rein_num/all_length*(hcx+2*distance_to_edge)),...
        ceil((hcx+2*distance_to_edge-distance_tip)/(1.5*d_eff)+1)]);

    Num_y = max([ceil(all_rein_num/all_length*(hcy+2*distance_to_edge)),...
        ceil((hcy+2*distance_to_edge-distance_tip)/(1.5*d_eff)+1)]);

    Downstand_sum_rein(i1).L_rein(1).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(1).x(2) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(1).y(1) = hcy/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(1).y(2) = hcy/2+distance_to_edge;
    
    if Num_x ~= 1
        tem = linspace...
            (-hcx/2 - distance_to_edge + distance_tip/2,...
            hcx/2 + distance_to_edge - distance_tip/2,Num_x);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    Downstand_sum_rein(i1).L_rein(1).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(1).num = Num_x;

    Downstand_sum_rein(i1).L_rein(1).x_mid = tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(1).x_mid);
    Downstand_sum_rein(i1).L_rein(1).y_mid = ...
                                    (hcy/2+distance_to_edge)*ones(1,num_ele);   

    Downstand_sum_rein(i1).L_rein(2).x(1) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(2).x(2) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(2).y(1) = hcy/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(2).y(2) = -hcy/2-distance_to_edge;

    if Num_y ~= 1
        tem = linspace...
            (-hcy/2 - distance_to_edge + distance_tip/2,...
            hcy/2 + distance_to_edge - distance_tip/2,Num_y);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    Downstand_sum_rein(i1).L_rein(2).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(2).num = Num_y;

    Downstand_sum_rein(i1).L_rein(2).y_mid = tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(2).y_mid);
    Downstand_sum_rein(i1).L_rein(2).x_mid = ...
                                    (hcx/2+distance_to_edge)*ones(1,num_ele);


    Downstand_sum_rein(i1).L_rein(3).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(3).x(2) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(3).y(1) = -hcy/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(3).y(2) = -hcy/2-distance_to_edge;

    if Num_x ~= 1
        tem = linspace...
            (-hcx/2 - distance_to_edge + distance_tip/2,...
            hcx/2 + distance_to_edge - distance_tip/2,Num_x);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    Downstand_sum_rein(i1).L_rein(3).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(3).num = Num_x;

    Downstand_sum_rein(i1).L_rein(3).x_mid = tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(3).x_mid);
    Downstand_sum_rein(i1).L_rein(3).y_mid = ...
                                    (-hcy/2-distance_to_edge)*ones(1,num_ele);


    Downstand_sum_rein(i1).L_rein(4).x(1) = -hcx/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(4).x(2) = -hcx/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(4).y(1) = hcy/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(4).y(2) = -hcy/2-distance_to_edge;

    if Num_y ~= 1
        tem = linspace...
            (-hcy/2 - distance_to_edge + distance_tip/2,...
            hcy/2 + distance_to_edge - distance_tip/2,Num_y);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    Downstand_sum_rein(i1).L_rein(4).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(4).num = Num_y;

    Downstand_sum_rein(i1).L_rein(4).y_mid = tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(4).y_mid);
    Downstand_sum_rein(i1).L_rein(4).x_mid = ...
                                    (-hcx/2-distance_to_edge)*ones(1,num_ele);




    % Intersection
    for i2 = 1:4
        num_ele = length(Downstand_sum_rein(i1).L_rein(i2).x_mid);
        intersec_tem = zeros(num_ele,number_of_hole);
        element_all = [Downstand_sum_rein(i1).L_rein(i2).x_mid',...
            Downstand_sum_rein(i1).L_rein(i2).y_mid'];

        for i3 = 1:number_of_hole
            shape_tem = Hole(i3).shape;
            in_tem = inpolygon(element_all(:,1),element_all(:,2)...
                ,shape_tem.Vertices(:,1),shape_tem.Vertices(:,2));
            intersec_tem(:,i3) = in_tem;
        end
        
        for i3 = 1:num_ele
            inter_sum = sum(intersec_tem(i3,:));
            if inter_sum == 0
                Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 0;
            else
                Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
            end

            if Col_position == 2
                if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                end
            elseif Col_position == 3
                if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                end

                if Downstand_sum_rein(i1).L_rein(i2).y_mid(i3) <=  -hcy_col/2
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                end
            end

        end
    end
end


% Adding more shear due to uneven drop

last_edge_distance = Downstand_sum_rein(Perimeter_rein_num).distance_to_edge;

if Last_peri == stand_peri_num && Drop_x/2 - hcx/2 - last_edge_distance >= 0.75*d_eff
    distance_4_add = Drop_x/2 - hcx/2 - last_edge_distance;
    Num_peri_add = floor(distance_4_add/(0.75*d_eff));
    for i1 = (Perimeter_rein_num+1):(Perimeter_rein_num+Num_peri_add)
        distance_to_edge = last_edge_distance + (i1-Perimeter_rein_num)*0.75*d_eff;
        Downstand_sum_rein(i1).distance_to_edge = distance_to_edge;

        ele_size = Asw_rein(Perimeter_rein_num,3);
        distance_tip = min(d_eff,ele_size);

        all_length = hcy + 2*last_edge_distance;
        all_rein_num = ceil(all_length*outter_peri_density/Area_link);
        
        Num_y = max([all_rein_num,...
        ceil((all_length-distance_tip)/(1.5*d_eff)+1)]);
    
        Downstand_sum_rein(i1).L_rein(1).x(1) = -hcx/2-distance_to_edge;
        Downstand_sum_rein(i1).L_rein(1).x(2) = hcx/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(1).y(1) = hcy/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(1).y(2) = hcy/2+last_edge_distance;
            
        Downstand_sum_rein(i1).L_rein(1).ele_size = 0;
        Downstand_sum_rein(i1).L_rein(1).num = 0;
    
        Downstand_sum_rein(i1).L_rein(1).x_mid = [-hcx/2-distance_to_edge,...
                                                   hcx/2+distance_to_edge];
        Downstand_sum_rein(i1).L_rein(1).y_mid = [hcy/2+last_edge_distance,...
                                                   hcy/2+last_edge_distance];  
    

        Downstand_sum_rein(i1).L_rein(2).x(1) = hcx/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(2).x(2) = hcx/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(2).y(1) = hcy/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(2).y(2) = -hcy/2-last_edge_distance;

        if Num_y ~= 1
            tem = linspace...
                (-hcy/2 - last_edge_distance + distance_tip/2,...
                hcy/2 + last_edge_distance - distance_tip/2,Num_y);
            ele_size_tem = abs(tem(1) - tem(2));
        else
            tem = 0;
            ele_size_tem = 0;
        end
    
        Downstand_sum_rein(i1).L_rein(2).ele_size = ele_size_tem;
        Downstand_sum_rein(i1).L_rein(2).num = Num_y;
    
        Downstand_sum_rein(i1).L_rein(2).y_mid = tem;
        num_ele = length(Downstand_sum_rein(i1).L_rein(2).y_mid);
        Downstand_sum_rein(i1).L_rein(2).x_mid = ...
                                        (hcx/2+distance_to_edge)*ones(1,num_ele);
    
    
        Downstand_sum_rein(i1).L_rein(3).x(1) = -hcx/2-distance_to_edge;
        Downstand_sum_rein(i1).L_rein(3).x(2) = hcx/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(3).y(1) = -hcy/2-last_edge_distance;
        Downstand_sum_rein(i1).L_rein(3).y(2) = -hcy/2-last_edge_distance;
        
        Downstand_sum_rein(i1).L_rein(3).ele_size = 0;
        Downstand_sum_rein(i1).L_rein(3).num = 0;
    
        Downstand_sum_rein(i1).L_rein(3).x_mid = [-hcx/2-distance_to_edge,...
                                                   hcx/2+distance_to_edge];
        Downstand_sum_rein(i1).L_rein(3).y_mid = [-hcy/2-last_edge_distance,...
                                                   -hcy/2-last_edge_distance];
    
    
        Downstand_sum_rein(i1).L_rein(4).x(1) = -hcx/2-distance_to_edge;
        Downstand_sum_rein(i1).L_rein(4).x(2) = -hcx/2-distance_to_edge;
        Downstand_sum_rein(i1).L_rein(4).y(1) = hcy/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(4).y(2) = -hcy/2-last_edge_distance;
    
        if Num_y ~= 1
            tem = linspace...
                (-hcy/2 - last_edge_distance + distance_tip/2,...
                hcy/2 + last_edge_distance - distance_tip/2,Num_y);
            ele_size_tem = abs(tem(1) - tem(2));
        else
            tem = 0;
            ele_size_tem = 0;
        end
    
        Downstand_sum_rein(i1).L_rein(4).ele_size = ele_size_tem;
        Downstand_sum_rein(i1).L_rein(4).num = Num_y;
    
        Downstand_sum_rein(i1).L_rein(4).y_mid = tem;
        num_ele = length(Downstand_sum_rein(i1).L_rein(4).y_mid);
        Downstand_sum_rein(i1).L_rein(4).x_mid = ...
                                        (-hcx/2-distance_to_edge)*ones(1,num_ele);
    
    
    
    
        % Intersection
        for i2 = 2:2:4
            num_ele = length(Downstand_sum_rein(i1).L_rein(i2).x_mid);
            intersec_tem = zeros(num_ele,number_of_hole);
            element_all = [Downstand_sum_rein(i1).L_rein(i2).x_mid',...
                Downstand_sum_rein(i1).L_rein(i2).y_mid'];
    
            for i3 = 1:number_of_hole
                shape_tem = Hole(i3).shape;
                in_tem = inpolygon(element_all(:,1),element_all(:,2)...
                    ,shape_tem.Vertices(:,1),shape_tem.Vertices(:,2));
                intersec_tem(:,i3) = in_tem;
            end
            
            for i3 = 1:num_ele
                inter_sum = sum(intersec_tem(i3,:));
                if inter_sum == 0
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 0;
                else
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                end


                if Col_position == 2
                    if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
                elseif Col_position == 3
                    if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
    
                    if Downstand_sum_rein(i1).L_rein(i2).y_mid(i3) <=  -hcy_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
                end

            end
        end

        for i2 = 1:2:3
            num_ele = length(Downstand_sum_rein(i1).L_rein(i2).x_mid);
            for i3 = 1:num_ele
                Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
            end
        end

    end
end




if Last_peri == stand_peri_num && Drop_y/2 - hcy/2 - last_edge_distance >= 0.75*d_eff
    distance_4_add = Drop_y/2 - hcy/2 - last_edge_distance;
    Num_peri_add = floor(distance_4_add/(0.75*d_eff));
    for i1 = (Perimeter_rein_num+1):(Perimeter_rein_num+Num_peri_add)
        distance_to_edge = last_edge_distance + (i1-Perimeter_rein_num)*0.75*d_eff;
        Downstand_sum_rein(i1).distance_to_edge = distance_to_edge;
    
        all_length = hcx + 2*last_edge_distance;
        all_rein_num = ceil(all_length*outter_peri_density/Area_link);

        ele_size = Asw_rein(Perimeter_rein_num,3);
        distance_tip = min(d_eff,ele_size);

        Num_x = max([all_rein_num,...
        ceil((all_length-distance_tip)/(1.5*d_eff)+1)]);
    
        Downstand_sum_rein(i1).L_rein(1).x(1) = -hcx/2-last_edge_distance;
        Downstand_sum_rein(i1).L_rein(1).x(2) = hcx/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(1).y(1) = hcy/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(1).y(2) = hcy/2+distance_to_edge;
        
        if Num_x ~= 1
            tem = linspace...
                (-hcx/2 - last_edge_distance + distance_tip/2,...
                hcx/2 + last_edge_distance - distance_tip/2,Num_x);
            ele_size_tem = abs(tem(1) - tem(2));
        else
            tem = 0;
            ele_size_tem = 0;
        end
    
        Downstand_sum_rein(i1).L_rein(1).ele_size = ele_size_tem;
        Downstand_sum_rein(i1).L_rein(1).num = Num_x;
    
        Downstand_sum_rein(i1).L_rein(1).x_mid = tem;
        num_ele = length(Downstand_sum_rein(i1).L_rein(1).x_mid);
        Downstand_sum_rein(i1).L_rein(1).y_mid = ...
                                        (hcy/2+distance_to_edge)*ones(1,num_ele);   
    

        Downstand_sum_rein(i1).L_rein(2).x(1) = hcx/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(2).x(2) = hcx/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(2).y(1) = hcy/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(2).y(2) = -hcy/2-distance_to_edge;
        
        Downstand_sum_rein(i1).L_rein(2).ele_size = 0;
        Downstand_sum_rein(i1).L_rein(2).num = 0;
    
        Downstand_sum_rein(i1).L_rein(2).y_mid = [hcy/2+distance_to_edge, -hcy/2-distance_to_edge];
        Downstand_sum_rein(i1).L_rein(2).x_mid = [hcx/2+last_edge_distance, hcx/2+last_edge_distance];
    
    
        Downstand_sum_rein(i1).L_rein(3).x(1) = -hcx/2-last_edge_distance;
        Downstand_sum_rein(i1).L_rein(3).x(2) = hcx/2+last_edge_distance;
        Downstand_sum_rein(i1).L_rein(3).y(1) = -hcy/2-distance_to_edge;
        Downstand_sum_rein(i1).L_rein(3).y(2) = -hcy/2-distance_to_edge;
    
        if Num_x ~= 1
            tem = linspace...
                (-hcx/2 - last_edge_distance + distance_tip/2,...
                hcx/2 + last_edge_distance - distance_tip/2,Num_x);
            ele_size_tem = abs(tem(1) - tem(2));
        else
            tem = 0;
            ele_size_tem = 0;
        end
    
        Downstand_sum_rein(i1).L_rein(3).ele_size = ele_size_tem;
        Downstand_sum_rein(i1).L_rein(3).num = Num_x;
    
        Downstand_sum_rein(i1).L_rein(3).x_mid = tem;
        num_ele = length(Downstand_sum_rein(i1).L_rein(3).x_mid);
        Downstand_sum_rein(i1).L_rein(3).y_mid = ...
                                        (-hcy/2-distance_to_edge)*ones(1,num_ele);
    
    

        Downstand_sum_rein(i1).L_rein(4).x(1) = -hcx/2-last_edge_distance;
        Downstand_sum_rein(i1).L_rein(4).x(2) = -hcx/2-last_edge_distance;
        Downstand_sum_rein(i1).L_rein(4).y(1) = hcy/2+distance_to_edge;
        Downstand_sum_rein(i1).L_rein(4).y(2) = -hcy/2-distance_to_edge;
       
        Downstand_sum_rein(i1).L_rein(4).ele_size = 0;
        Downstand_sum_rein(i1).L_rein(4).num = 0;
    
        Downstand_sum_rein(i1).L_rein(4).y_mid = [-hcy/2-distance_to_edge, hcy/2+distance_to_edge];
        Downstand_sum_rein(i1).L_rein(4).x_mid = [-hcx/2-last_edge_distance, -hcx/2-last_edge_distance];
    
    
    
    
        % Intersection
        for i2 = 1:2:3
            num_ele = length(Downstand_sum_rein(i1).L_rein(i2).x_mid);
            intersec_tem = zeros(num_ele,number_of_hole);
            element_all = [Downstand_sum_rein(i1).L_rein(i2).x_mid',...
                Downstand_sum_rein(i1).L_rein(i2).y_mid'];
    
            for i3 = 1:number_of_hole
                shape_tem = Hole(i3).shape;
                in_tem = inpolygon(element_all(:,1),element_all(:,2)...
                    ,shape_tem.Vertices(:,1),shape_tem.Vertices(:,2));
                intersec_tem(:,i3) = in_tem;
            end
            
            for i3 = 1:num_ele
                inter_sum = sum(intersec_tem(i3,:));
                if inter_sum == 0
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 0;
                else
                    Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                end

                if Col_position == 2
                    if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
                elseif Col_position == 3
                    if Downstand_sum_rein(i1).L_rein(i2).x_mid(i3) <=  -hcx_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
    
                    if Downstand_sum_rein(i1).L_rein(i2).y_mid(i3) <=  -hcy_col/2
                        Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
                    end
                end
                
            end
        end

        for i2 = 2:2:4
            num_ele = length(Downstand_sum_rein(i1).L_rein(i2).x_mid);
            
            for i3 = 1:num_ele
                Downstand_sum_rein(i1).L_rein(i2).intersect(i3) = 1;
            end
        end
    end
end