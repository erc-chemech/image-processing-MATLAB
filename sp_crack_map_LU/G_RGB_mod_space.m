R_space=linspace(-0.1,0.05,1e3);
G_space=linspace(-0.1,0.05,1e3);

RGB_mod=@(R,G) sqrt(R.^2+G.^2+(-R-G).^2);

RGB=nan(length(G_space),length(R_space));

for dum=1:length(R_space)
    for dum_1=1:length(G_space)
        RGB(dum_1,dum)=RGB_mod(R_space(dum),G_space(dum_1));
    end
end

RGB2=[];
for dum=1:length(G_space)
    RGB2=[RGB2;min(RGB(dum,:))];
end
clear RGB;

