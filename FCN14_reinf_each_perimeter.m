function [Downstand_sum,Downrein_sum] = FCN14_reinf_each_perimeter...
        (Downstand_sum, Downrein_sum)


stand_peri_num = length(Downstand_sum);
rein_peri_num = length(Downrein_sum);

for i1 = 1:rein_peri_num

    tem_1 = ...
        [Downrein_sum(i1).L_rein(1).x_mid';...
        Downrein_sum(i1).L_rein(2).x_mid';...
        Downrein_sum(i1).L_rein(3).x_mid';...
        Downrein_sum(i1).L_rein(4).x_mid'];

    tem_2 = ...
        [Downrein_sum(i1).L_rein(1).y_mid';...
        Downrein_sum(i1).L_rein(2).y_mid';...
        Downrein_sum(i1).L_rein(3).y_mid';...
        Downrein_sum(i1).L_rein(4).y_mid'];

    tem_3 = ...
        [Downrein_sum(i1).L_rein(1).intersect';...
        Downrein_sum(i1).L_rein(2).intersect';...
        Downrein_sum(i1).L_rein(3).intersect';...
        Downrein_sum(i1).L_rein(4).intersect'];

    Downrein_sum(i1).all_reinf = unique([tem_1, tem_2, tem_3],'row');
end

for i0 = 2:stand_peri_num
    poly_gon = Downstand_sum(i0).poly_gon;
    poly_full = Downstand_sum(i0).poly_full;

    Downstand_sum(i0).all_rein_length = 0;

    Downstand_sum(i0).shear_peri = [];

    for i1 = 1:rein_peri_num
        num_ele = length(Downrein_sum(i1).all_reinf);
        element_all = Downrein_sum(i1).all_reinf(:,1:2);

        in_tem = inpolygon(element_all(:,1),element_all(:,2)...
            ,poly_gon.Vertices(:,1),poly_gon.Vertices(:,2));

        [~,on_tem] = inpolygon(element_all(:,1),element_all(:,2)...
            ,poly_full.Vertices(:,1),poly_full.Vertices(:,2));

        in_tem = in_tem - on_tem;

        Downstand_sum(i0).rein_consider(i1).intersect = zeros(num_ele,2);
        Downstand_sum(i0).rein_consider(i1).intersect(:,2) = Downrein_sum(i1).all_reinf(:,3);

        for i3 = 1:num_ele
            if in_tem(i3,1) == 1 && ...
                    Downstand_sum(i0).rein_consider(i1).intersect(i3,2) == 0
                Downstand_sum(i0).rein_consider(i1).intersect(i3,1) = 1;
            else
                Downstand_sum(i0).rein_consider(i1).intersect(i3,1) = 0;
            end
        end

        Downstand_sum(i0).rein_consider(i1).total = ...
            sum(Downstand_sum(i0).rein_consider(i1).intersect(:,1));

        if Downstand_sum(i0).rein_consider(i1).total ~= 0
            Downstand_sum(i0).shear_peri = [Downstand_sum(i0).shear_peri,i1];
        end
    end
end
