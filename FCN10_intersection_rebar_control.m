function Downstand_sum = FCN10_intersection_rebar_control...
        (Downstand_sum, Downrein_sum)

stand_peri_num = length(Downstand_sum);
rein_peri_num = length(Downrein_sum);

for i0 = 2:stand_peri_num
    poly_gon = Downstand_sum(i0).poly_gon;
    poly_full = Downstand_sum(i0).poly_full;

    Downstand_sum(i0).all_rein_length = 0;

    Downstand_sum(i0).shear_peri = [];

    for i1 = 1:rein_peri_num
        num_ele = length(Downrein_sum(i1).All_element);
        element_all = Downrein_sum(i1).All_element(:,1:2);

        in_tem = inpolygon(element_all(:,1),element_all(:,2)...
            ,poly_gon.Vertices(:,1),poly_gon.Vertices(:,2));

        [~,on_tem] = inpolygon(element_all(:,1),element_all(:,2)...
            ,poly_full.Vertices(:,1),poly_full.Vertices(:,2));

        in_tem = in_tem-on_tem;

        Downstand_sum(i0).rein_in_peri(i1).intersect = zeros(num_ele,2);
        Downstand_sum(i0).rein_in_peri(i1).intersect(:,2) = Downrein_sum(i1).All_element(:,3);

        for i3 = 1:num_ele
            if in_tem(i3,1) == 1 && Downrein_sum(i1).All_element(i3,4) == 0
                Downstand_sum(i0).rein_in_peri(i1).intersect(i3,1) = 1;
            else
                Downstand_sum(i0).rein_in_peri(i1).intersect(i3,1) = 0;
            end
        end

        Downstand_sum(i0).rein_in_peri(i1).length_each = ...
            sum(Downstand_sum(i0).rein_in_peri(i1).intersect(:,2).*...
                Downstand_sum(i0).rein_in_peri(i1).intersect(:,1));

        if Downstand_sum(i0).rein_in_peri(i1).length_each ~= 0
            Downstand_sum(i0).shear_peri = [Downstand_sum(i0).shear_peri,i1];
        end

        Downstand_sum(i0).all_rein_length = ...
            Downstand_sum(i0).all_rein_length ...
            + Downstand_sum(i0).rein_in_peri(i1).length_each;
    end
end
