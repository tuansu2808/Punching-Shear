function Downstand_sum_rein = FCN13_Arrangement_of_shear_wall...
    (Asw_rein,Downstand_sum_rein, d_eff, hcx, hcy, ...
    Hole, Perimeter_rein_num, number_of_hole, Area_link,...
    Downstand_sum, ~, Last_peri, outter_peri_density,...
    Drop_x, Drop_y)


% Asw_rein = Asw_drop_rein;
% Downstand_sum_rein = Drop_sum_rein;
% d_eff = d_tem;
% hcx = hcx;
% hcy = hcy;
% Hole = Hole;
% Perimeter_rein_num = Drop_perimeter_rein_num;
% number_of_hole = number_of_hole_input;
% Area_link = Area_link;
% Downstand_sum = Drop_sum;
% Last_peri = Last_drop;
% outter_peri_density = outter_peri_density;
% Drop_x = Drop_x;
% Drop_y = Drop_y;



% Gross length of each each drop perimeters
stand_peri_num = length(Downstand_sum);

for i1 = 1:Perimeter_rein_num

    distance_to_edge = Downstand_sum_rein(i1).distance_to_edge;

    ele_size = Asw_rein(i1,3);

    distance_tip = min(d_eff,ele_size);

    length_x = distance_to_edge + 1.5*d_eff + distance_to_edge/(3^0.5);
    length_y = hcy + 2*distance_to_edge;

    all_length = length_y + 2*length_x;

    all_rein_num = ceil(all_length/ele_size);

    Num_x = max([ceil(all_rein_num/all_length*length_x),...
                ceil((length_x-distance_tip)/(1.5*d_eff)+1)]);

    Num_y = max([ceil(all_rein_num/all_length*length_y),...
                ceil((length_y-distance_tip)/(1.5*d_eff)+1)]);

    %%
    Downstand_sum_rein(i1).L_rein(1).x(1) = hcx/2-1.5*d_eff-distance_to_edge/(3^0.5);
    Downstand_sum_rein(i1).L_rein(1).x(2) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(1).y(1) = hcy/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(1).y(2) = hcy/2+distance_to_edge;
    
    Length_tem = length_x;

    if Num_x ~= 1
        tem = linspace...
            (distance_tip/2,Length_tem-distance_tip/2,Num_x);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    Downstand_sum_rein(i1).L_rein(1).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(1).num = Num_x;

    first_point = Downstand_sum_rein(i1).L_rein(1).x(1);
    num_tem_2 = length(tem);
    first_point_vector = first_point*ones(1,num_tem_2);

    Downstand_sum_rein(i1).L_rein(1).x_mid = first_point_vector+tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(1).x_mid);
    Downstand_sum_rein(i1).L_rein(1).y_mid = ...
                                    (hcy/2+distance_to_edge)*ones(1,num_ele);   

    %%
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

    Downstand_sum_rein(i1).L_rein(2).y_mid = unique([-tem,tem]);
    num_ele = length(Downstand_sum_rein(i1).L_rein(2).y_mid);
    Downstand_sum_rein(i1).L_rein(2).x_mid = ...
                                    (hcx/2+distance_to_edge)*ones(1,num_ele);


    %%
    Downstand_sum_rein(i1).L_rein(3).x(1) = hcx/2-1.5*d_eff-distance_to_edge/(3^0.5);
    Downstand_sum_rein(i1).L_rein(3).x(2) = hcx/2+distance_to_edge;
    Downstand_sum_rein(i1).L_rein(3).y(1) = -hcy/2-distance_to_edge;
    Downstand_sum_rein(i1).L_rein(3).y(2) = -hcy/2-distance_to_edge;

    Length_tem = Downstand_sum_rein(i1).L_rein(3).x(2) - ...
        Downstand_sum_rein(i1).L_rein(3).x(1);

    if Num_x ~= 1
        tem = linspace...
            (distance_tip/2,Length_tem-distance_tip/2,Num_x);
        ele_size_tem = abs(tem(1) - tem(2));
    else
        tem = 0;
        ele_size_tem = 0;
    end

    first_point = Downstand_sum_rein(i1).L_rein(3).x(1);
    num_tem_2 = length(tem);
    first_point_vector = first_point*ones(1,num_tem_2);

    Downstand_sum_rein(i1).L_rein(3).ele_size = ele_size_tem;
    Downstand_sum_rein(i1).L_rein(3).num = Num_x;

    Downstand_sum_rein(i1).L_rein(3).x_mid = first_point_vector+tem;
    num_ele = length(Downstand_sum_rein(i1).L_rein(3).x_mid);
    Downstand_sum_rein(i1).L_rein(3).y_mid = ...
                                    (-hcy/2-distance_to_edge)*ones(1,num_ele);


    %%
    Downstand_sum_rein(i1).L_rein(4).x(1) = Downstand_sum_rein(i1).L_rein(3).x(1);
    Downstand_sum_rein(i1).L_rein(4).x(2) = Downstand_sum_rein(i1).L_rein(3).x(2);
    Downstand_sum_rein(i1).L_rein(4).y(1) = Downstand_sum_rein(i1).L_rein(3).y(1);
    Downstand_sum_rein(i1).L_rein(4).y(2) = Downstand_sum_rein(i1).L_rein(3).y(2);

    Downstand_sum_rein(i1).L_rein(4).ele_size = 0;
    Downstand_sum_rein(i1).L_rein(4).num = 1;

    Downstand_sum_rein(i1).L_rein(4).y_mid = Downstand_sum_rein(i1).L_rein(3).y_mid(end);
    Downstand_sum_rein(i1).L_rein(4).x_mid = Downstand_sum_rein(i1).L_rein(3).x_mid(end);

end

% Intersection

for i1 = 1:Perimeter_rein_num
    for i2 = 1:3
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
        end
    end

    Downstand_sum_rein(i1).L_rein(4).intersect(1) = 1;

end




if outter_peri_density ~= 0

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
            
            Num_y = all_rein_num;
        
            %%
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
        
    
            %%
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
        
        
            %%
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
        
        
            %%
            Downstand_sum_rein(i1).L_rein(4).x(1) = -hcx/2-distance_to_edge;
            Downstand_sum_rein(i1).L_rein(4).x(2) = -hcx/2-distance_to_edge;
            Downstand_sum_rein(i1).L_rein(4).y(1) = hcy/2+last_edge_distance;
            Downstand_sum_rein(i1).L_rein(4).y(2) = -hcy/2-last_edge_distance;
        
            Downstand_sum_rein(i1).L_rein(4).ele_size = 0;
            Downstand_sum_rein(i1).L_rein(4).num = 0;
        
            Downstand_sum_rein(i1).L_rein(4).y_mid = [-hcx/2-distance_to_edge,...
                                                      -hcx/2-distance_to_edge];
            Downstand_sum_rein(i1).L_rein(4).x_mid = [-hcy/2-distance_to_edge,...
                                                       hcy/2+distance_to_edge];
        
        
        
        
            % Intersection
            for i2 = 2
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
                end
            end
    
            for i2 = [1,3,4]
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
        
            all_length = last_edge_distance + 1.5*d_eff + last_edge_distance/(3^0.5);
            all_rein_num = ceil(all_length*outter_peri_density/Area_link);
    
            ele_size = Asw_rein(Perimeter_rein_num,3);
            distance_tip = min(d_eff,ele_size);
    
            Num_x = all_rein_num;
        
            %%
            Downstand_sum_rein(i1).L_rein(1).x(1) = hcx/2-1.5*d_eff-last_edge_distance/(3^0.5);
            Downstand_sum_rein(i1).L_rein(1).x(2) = hcx/2+last_edge_distance;
            Downstand_sum_rein(i1).L_rein(1).y(1) = hcy/2+distance_to_edge;
            Downstand_sum_rein(i1).L_rein(1).y(2) = hcy/2+distance_to_edge;
            
            Length_tem = Downstand_sum_rein(i1).L_rein(1).x(2) - ...
                Downstand_sum_rein(i1).L_rein(1).x(1);
        
            if Num_x ~= 1
                tem = linspace...
                    (distance_tip/2,Length_tem-distance_tip/2,Num_x);
                ele_size_tem = abs(tem(1) - tem(2));
            else
                tem = 0;
                ele_size_tem = 0;
            end
        
            Downstand_sum_rein(i1).L_rein(1).ele_size = ele_size_tem;
            Downstand_sum_rein(i1).L_rein(1).num = Num_x;
        
            first_point = Downstand_sum_rein(i1).L_rein(1).x(1);
            num_tem_2 = length(tem);
            first_point_vector = first_point*ones(1,num_tem_2);
        
            Downstand_sum_rein(i1).L_rein(1).x_mid = first_point_vector+tem;
            num_ele = length(Downstand_sum_rein(i1).L_rein(1).x_mid);
            Downstand_sum_rein(i1).L_rein(1).y_mid = ...
                                            (hcy/2+distance_to_edge)*ones(1,num_ele);   
            
            %%
    
            Downstand_sum_rein(i1).L_rein(2).x(1) = hcx/2+last_edge_distance;
            Downstand_sum_rein(i1).L_rein(2).x(2) = hcx/2+last_edge_distance;
            Downstand_sum_rein(i1).L_rein(2).y(1) = hcy/2+distance_to_edge;
            Downstand_sum_rein(i1).L_rein(2).y(2) = -hcy/2-distance_to_edge;
            
            Downstand_sum_rein(i1).L_rein(2).ele_size = 0;
            Downstand_sum_rein(i1).L_rein(2).num = 0;
        
            Downstand_sum_rein(i1).L_rein(2).y_mid = [hcy/2+distance_to_edge, -hcy/2-distance_to_edge];
            Downstand_sum_rein(i1).L_rein(2).x_mid = [hcx/2+last_edge_distance, hcx/2+last_edge_distance];
        
            %%
            Downstand_sum_rein(i1).L_rein(3).x(1) = hcx/2-1.5*d_eff-last_edge_distance/(3^0.5);
            Downstand_sum_rein(i1).L_rein(3).x(2) = hcx/2+last_edge_distance;
            Downstand_sum_rein(i1).L_rein(3).y(1) = -hcy/2-distance_to_edge;
            Downstand_sum_rein(i1).L_rein(3).y(2) = -hcy/2-distance_to_edge;
        
            Length_tem = Downstand_sum_rein(i1).L_rein(3).x(2) - ...
                Downstand_sum_rein(i1).L_rein(3).x(1);
        
            if Num_x ~= 1
                tem = linspace...
                    (distance_tip/2,Length_tem-distance_tip/2,Num_x);
                ele_size_tem = abs(tem(1) - tem(2));
            else
                tem = 0;
                ele_size_tem = 0;
            end
        
            first_point = Downstand_sum_rein(i1).L_rein(3).x(1);
            num_tem_2 = length(tem);
            first_point_vector = first_point*ones(1,num_tem_2);
        
            Downstand_sum_rein(i1).L_rein(3).ele_size = ele_size_tem;
            Downstand_sum_rein(i1).L_rein(3).num = num_tem_2;
        
            Downstand_sum_rein(i1).L_rein(3).x_mid = first_point_vector+tem;
            num_ele = length(Downstand_sum_rein(i1).L_rein(3).x_mid);
            Downstand_sum_rein(i1).L_rein(3).y_mid = ...
                                            (-hcy/2-distance_to_edge)*ones(1,num_ele);    
        
            %%
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
end