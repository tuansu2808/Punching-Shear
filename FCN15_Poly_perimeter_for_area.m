function [down_stand_sum,drop_area] = FCN15_Poly_perimeter_for_area...
    (down_stand_sum, number_of_hole, Hole, Col_position,...
    hcx, hcy, Drop_x_original, Drop_y_original, drop_or_slab)                         

% down_stand_sum = Slab_sum;
% number_of_hole = number_of_hole_input;



%%_column or drop area:

if drop_or_slab == 0
    Column_drop_area = hcx*hcy;
    down_stand_sum(1).area_inner = 0;
else
    if Col_position == 1
        col_drop = polyshape([-Drop_x_original/2,-Drop_x_original/2,Drop_x_original/2,Drop_x_original/2,-Drop_x_original/2],...
                             [-Drop_y_original/2,Drop_y_original/2,Drop_y_original/2,-Drop_y_original/2,-Drop_y_original/2,]);
    elseif Col_position == 2
        col_drop = polyshape([-hcx/2,-hcx/2,Drop_x_original-hcx/2,Drop_x_original-hcx/2,-hcx/2],...
                             [-Drop_y_original/2,Drop_y_original/2,Drop_y_original/2,-Drop_y_original/2,-Drop_y_original/2,]);
    elseif Col_position == 3
        col_drop = polyshape([-hcx/2,-hcx/2,Drop_x_original-hcx/2,Drop_x_original-hcx/2,-hcx/2],...
                             [-hcy/2,Drop_y_original-hcy/2,Drop_y_original-hcy/2,-hcy/2,-hcy/2,]);
    end
    
    Column_drop_area_full = polyarea(col_drop.Vertices(:,1),col_drop.Vertices(:,2));
    
    area_hole_col_drop = zeros(number_of_hole,1);
    for i2 = 1:number_of_hole
        shape = Hole(i2).shape;
        tem = subtract(col_drop,shape);
        area_hole_col_drop(i2,1) = Column_drop_area_full - polyarea(tem.Vertices(:,1),tem.Vertices(:,2));
    end

    Column_drop_area = Column_drop_area_full - sum(area_hole_col_drop);
    down_stand_sum(1).area_inner = Column_drop_area;
end




peri_num = length(down_stand_sum);

for i1 = 2:peri_num
    ele_tem_f = [down_stand_sum(i1).L_gross(1).x_mid', down_stand_sum(i1).L_gross(1).y_mid'];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(2).x_mid', down_stand_sum(i1).C_gross(2).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).L_gross(2).x_mid', down_stand_sum(i1).L_gross(2).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(3).x_mid', down_stand_sum(i1).C_gross(3).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).L_gross(3).x_mid', down_stand_sum(i1).L_gross(3).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(4).x_mid', down_stand_sum(i1).C_gross(4).y_mid'])];
    ele_tem_f = [ele_tem_f;[down_stand_sum(i1).L_gross(4).x_mid', down_stand_sum(i1).L_gross(4).y_mid']];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(1).x_mid', down_stand_sum(i1).C_gross(1).y_mid'])];
    
    tem = polyshape(ele_tem_f(:,1),ele_tem_f(:,2));


    
    % if Col_position == 1
    %     area_cut = 0;
    % else
        
    if Col_position == 2
        tem = polyshape(ele_tem_f(:,1),ele_tem_f(:,2));
        shape = polyshape([-hcx/2,-hcx/2,-100000,-100000],...
                          [-100000,100000,100000,-100000]);
        tem = subtract(tem,shape);
        % area_cut = area_full - polyarea(tem.Vertices(:,1),tem.Vertices(:,2));
    elseif Col_position == 3
        tem = polyshape(ele_tem_f(:,1),ele_tem_f(:,2));
        shape = polyshape([100000,100000,-hcx/2,-hcx/2,-100000,-100000],...
                          [-100000,-hcy/2,-hcy/2,100000,100000,-100000]);
        tem = subtract(tem,shape);
        % area_cut = area_full - polyarea(tem.Vertices(:,1),tem.Vertices(:,2));
    end


    tem_full = tem;

    area_full = polyarea(tem_full.Vertices(:,1),tem_full.Vertices(:,2));

    area_hole = zeros(number_of_hole,1);
    for i2 = 1:number_of_hole
        shape = Hole(i2).shape;
        tem = subtract(tem_full,shape);
        area_hole(i2,1) = area_full - polyarea(tem.Vertices(:,1),tem.Vertices(:,2));
    end



    % area_inner = area_full - sum(area_hole) - area_cut - Column_drop_area;
    area_inner = area_full - sum(area_hole) - Column_drop_area;


    down_stand_sum(i1).area_inner = area_inner;
    down_stand_sum(i1).poly_full = tem;
end


drop_area = Column_drop_area;