function Report_data = FCN10_report_column(Completed_input,Completed_output)

Report_data = cell(200,20);

i = 1;
Report_data(i,1) = {Completed_input.fck};
Report_data(i,2) = {Completed_input.gamma_c};
Report_data(i,3) = {Completed_input.fyk};
Report_data(i,4) = {Completed_input.gamma_s};
Report_data(i,5) = {Completed_input.Dia_link};
Report_data(i,6) = {Completed_input.One_stran_area};
Report_data(i,7) = {Completed_input.f_pe*Completed_input.fpu};
Report_data(i,8) = {Completed_input.gamma_s};

i = i + 1;
Report_data(i,1) = {Completed_input.Span_x/1000};
Report_data(i,2) = {Completed_input.Span_y/1000};
Report_data(i,3) = {Completed_input.hcx/1000};
Report_data(i,4) = {Completed_input.hcy/1000};
Report_data(i,5) = {Completed_input.Drop_x/1000};
Report_data(i,6) = {Completed_input.Drop_y/1000};
Report_data(i,7) = {Completed_input.S_thickness/1000};
Report_data(i,8) = {Completed_input.Drop_thickness/1000};
Report_data(i,9) = {(Completed_input.Cx + Completed_input.Cy)/2};
Report_data(i,10) = {(Completed_input.Cpx + Completed_input.Cpy)/2};

i = i + 1;
Report_data(i,1) = {Completed_input.Total_p_x};
Report_data(i,2) = {Completed_input.Num_p_x_col};
Report_data(i,3) = {Completed_input.Num_p_x_drop};
Report_data(i,4) = {Completed_input.Drape_x};

Report_data(i,5) = {Completed_input.Total_p_y};
Report_data(i,6) = {Completed_input.Num_p_y_col};
Report_data(i,7) = {Completed_input.Num_p_y_drop};
Report_data(i,8) = {Completed_input.Drape_y};


i = i + 1;
Report_data(i,1) = {Completed_input.Dia_x};
Report_data(i,2) = {Completed_input.Dia_x_slab};
Report_data(i,3) = {Completed_input.spacing_x_drop};
Report_data(i,4) = {Completed_input.Spacing_slab_x};

Report_data(i,5) = {Completed_input.Dia_y};
Report_data(i,6) = {Completed_input.Dia_y_slab};
Report_data(i,7) = {Completed_input.spacing_y_drop};
Report_data(i,8) = {Completed_input.Spacing_slab_y};

for i2 = 1:8
    i = i + 1;
    if  i2 <= size(Completed_input.hole_table,1)
        Report_data(i,1) = {Completed_input.hole_table(i2,2)};
        Report_data(i,2) = {Completed_input.hole_table(i2,3)}; 
        Report_data(i,3) = {Completed_input.hole_table(i2,4)};
        Report_data(i,4) = {Completed_input.hole_table(i2,5)};
    end
end

for i2 = 1:5
    i = i + 1;
    Report_data(i,1) = {Completed_input.Rz(i2,1)};
    Report_data(i,2) = {Completed_input.Rz(i2,2)};
    Report_data(i,3) = {Completed_input.Rz(i2,3)};
    Report_data(i,4) = {Completed_input.factor_SDL_LL(i2,1)};
    Report_data(i,5) = {Completed_input.factor_SDL_LL(i2,2)};
    Report_data(i,6) = {Completed_input.ex(i2,1)};
    Report_data(i,7) = {Completed_input.ey(i2,1)};
end

i = i + 1;
Report_data(i,1) = {Completed_output.stress_tem_drop(1,1)};
Report_data(i,2) = {Completed_output.V_Rd_max};

for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_drop

        Rz_max = max(Completed_output.Drop_sum(i2).Rz_tem);
        f_tem = find(Completed_output.Drop_sum(i2).Rz_tem == Rz_max);
        beta = Completed_output.Drop_sum(i2).beta_tem(f_tem);

        Report_data(i,1) = {Completed_output.drop_table{i2,3}};
        Report_data(i,2) = {Completed_output.drop_table{i2,6}};
        Report_data(i,3) = {Completed_output.drop_table{i2,7}};
        Report_data(i,4) = {Rz_max/beta};
        Report_data(i,5) = {Completed_output.Drop_sum(i2).bx};
        Report_data(i,6) = {Completed_output.Drop_sum(i2).by};
        Report_data(i,7) = {beta};
        Report_data(i,8) = {Rz_max};
        Report_data(i,9) = {Completed_output.drop_table{i2,9}};
        
    end
end


for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_drop

        Report_data(i,1) = {Completed_output.drop_table{i2,2}};
        Report_data(i,2) = Completed_output.drop_table{i2,4};
        Report_data(i,3) = {Completed_output.drop_table{i2,10}};

    end
end


for i2 = 1:14
    i = i + 1;
    try
        if  i2 <= size(Completed_output.drop_rein_table{:,:},1)
            Report_data(i,1) = {Completed_output.drop_rein_table{i2,2}};
            Report_data(i,2) = {Completed_output.drop_rein_table{i2,3}};
            Report_data(i,3) = {Completed_output.drop_rein_table{i2,4}};
            Report_data(i,4) = {Completed_output.drop_rein_table{i2,5}};
            Report_data(i,5) = {Completed_output.drop_rein_table{i2,6}};
            Report_data(i,6) = {Completed_output.drop_rein_table{i2,7}};
            Report_data(i,7) = {Completed_output.drop_rein_table{i2,8}};
            Report_data(i,8) = {Completed_output.drop_rein_table{i2,9}};
        end
    end
end





i = i + 1;
Report_data(i,1) = {Completed_output.stress_tem_slab(1,1)};
Report_data(i,2) = {Completed_output.V_Rd_max};
Report_data(i,3) = {2*(Completed_output.w_x_drop + Completed_output.w_y_drop)/1000};

for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_slab

        Rz_max = max(Completed_output.Slab_sum(i2).Rz_tem);
        f_tem = find(Completed_output.Slab_sum(i2).Rz_tem == Rz_max);
        beta = Completed_output.Slab_sum(i2).beta_tem(f_tem);

        Report_data(i,1) = {Completed_output.slab_table{i2,3}};
        Report_data(i,2) = {Completed_output.slab_table{i2,6}};
        Report_data(i,3) = {Completed_output.slab_table{i2,7}};
        Report_data(i,4) = {Rz_max/beta};
        Report_data(i,5) = {Completed_output.Slab_sum(i2).bx};
        Report_data(i,6) = {Completed_output.Slab_sum(i2).by};
        Report_data(i,7) = {beta};
        Report_data(i,8) = {Rz_max};
        Report_data(i,9) = {Completed_output.slab_table{i2,9}};
        
    end
end


for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_slab

        Report_data(i,1) = {Completed_output.slab_table{i2,2}};
        Report_data(i,2) = Completed_output.slab_table{i2,4};
        Report_data(i,3) = {Completed_output.slab_table{i2,10}};

    end
end


for i2 = 1:14
    i = i + 1;
    try
        if  i2 <= size(Completed_output.slab_rein_table{:,:},1)
            Report_data(i,1) = {Completed_output.slab_rein_table{i2,2}};
            Report_data(i,2) = {Completed_output.slab_rein_table{i2,3}};
            Report_data(i,3) = {Completed_output.slab_rein_table{i2,4}};
            Report_data(i,4) = {Completed_output.slab_rein_table{i2,5}};
            Report_data(i,5) = {Completed_output.slab_rein_table{i2,6}};
            Report_data(i,6) = {Completed_output.slab_rein_table{i2,7}};
            Report_data(i,7) = {Completed_output.slab_rein_table{i2,8}};
            Report_data(i,8) = {Completed_output.slab_rein_table{i2,9}};
        end
    end
end

i = i + 1;
Report_data(i,1) = {Completed_output.Data_name};
Report_data(i,2) = {Completed_input.Col_position};
Report_data(i,3) = {Completed_input.Design_column_or_wall};


for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_drop
        if Completed_input.Col_position == 1
            Report_data(i,1) = {Completed_output.Drop_sum(i2).area_inner};
        elseif Completed_input.Col_position == 2
            Report_data(i,2) = {Completed_output.Drop_sum(i2).area_inner/2};
        elseif Completed_input.Col_position == 3
            Report_data(i,3) = {Completed_output.Drop_sum(i2).area_inner/4};
        end        
    end
end

for i2 = 2:10
    i = i + 1;
    if  i2 <= Completed_output.Last_drop
        if Completed_input.Col_position == 1
            Report_data(i,1) = {Completed_output.Slab_sum(i2).area_inner};
        elseif Completed_input.Col_position == 2
            Report_data(i,2) = {Completed_output.Slab_sum(i2).area_inner/2};
        elseif Completed_input.Col_position == 3
            Report_data(i,3) = {Completed_output.Slab_sum(i2).area_inner/4};
        end        
    end
end

i = i + 1;
Report_data(i,1) = {Completed_input.LL};
Report_data(i,2) = {Completed_input.SDL};
