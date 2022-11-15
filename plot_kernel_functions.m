close all
t=linspace(0,40,40*100);
sig_f = 8;
sig_l = 2;
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); 
%%
figure
subplot(211)
hold all
i_vec = 1:1:40;
j_vec=i_vec;
for i=1:40
    k_i = @(jj) k(i,jj);
    k_i_val = k_i(t);
    plot(t, k_i_val,'b')
    area(t, k_i_val,'FaceAlpha',0.25)
end

%%
sig_f = 8;
sig_l = 1/2;

k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); 
%%
subplot(212)
hold all
for i=1:40
    k_i2 = @(jj) k(i,jj);
    k_i2_val = k_i2(t);
    plot(t, k_i2_val,'r')
    area(t, k_i2_val,'FaceAlpha',0.25)
end