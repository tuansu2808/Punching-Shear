function down_stand_sum = FCN9_Poly_perimeter(down_stand_sum, peri_num, number_of_hole, Hole)                         

for i1 = 2:peri_num
    ele_tem_f = [down_stand_sum(i1).L_gross(1).x_mid', down_stand_sum(i1).L_gross(1).y_mid'];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(2).x_mid', down_stand_sum(i1).C_gross(2).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).L_gross(2).x_mid', down_stand_sum(i1).L_gross(2).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(3).x_mid', down_stand_sum(i1).C_gross(3).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).L_gross(3).x_mid', down_stand_sum(i1).L_gross(3).y_mid'])];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(4).x_mid', down_stand_sum(i1).C_gross(4).y_mid'])];
    ele_tem_f = [ele_tem_f;[down_stand_sum(i1).L_gross(4).x_mid', down_stand_sum(i1).L_gross(4).y_mid']];
    ele_tem_f = [ele_tem_f;flipud([down_stand_sum(i1).C_gross(1).x_mid', down_stand_sum(i1).C_gross(1).y_mid'])];
    
    if i1 == 2
        ele_tem_f_2 = [down_stand_sum(1).L_gross(1).x_mid', down_stand_sum(1).L_gross(1).y_mid'];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(1).L_gross(2).x_mid', down_stand_sum(1).L_gross(2).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(1).L_gross(3).x_mid', down_stand_sum(1).L_gross(3).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;[down_stand_sum(1).L_gross(4).x_mid', down_stand_sum(1).L_gross(4).y_mid']];
    else
        ele_tem_f_2 = [down_stand_sum(i1).L_gross_inner(1).x_mid', down_stand_sum(i1).L_gross_inner(1).y_mid'];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).C_gross_inner(2).x_mid', down_stand_sum(i1).C_gross_inner(2).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).L_gross_inner(2).x_mid', down_stand_sum(i1).L_gross_inner(2).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).C_gross_inner(3).x_mid', down_stand_sum(i1).C_gross_inner(3).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).L_gross_inner(3).x_mid', down_stand_sum(i1).L_gross_inner(3).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).C_gross_inner(4).x_mid', down_stand_sum(i1).C_gross_inner(4).y_mid'])];
        ele_tem_f_2 = [ele_tem_f_2;[down_stand_sum(i1).L_gross_inner(4).x_mid', down_stand_sum(i1).L_gross_inner(4).y_mid']];
        ele_tem_f_2 = [ele_tem_f_2;flipud([down_stand_sum(i1).C_gross_inner(1).x_mid', down_stand_sum(i1).C_gross_inner(1).y_mid'])];
    end

    tem = polyshape(ele_tem_f(:,1),ele_tem_f(:,2));

    tem_2 = polyshape(ele_tem_f_2(:,1),ele_tem_f_2(:,2));
    
    s_tem = subtract(tem,tem_2);
    
% hold on    
% plot(tem);
% plot(tem_2);
% plot(s_tem);


    for i2 = 1:number_of_hole
        shape = Hole(i2).shape;
        s_tem = subtract(s_tem,shape);
    end

    down_stand_sum(i1).poly_gon = s_tem;
end
