


h_s = @(s) 1/(s+1)*(0.3*s+1);
s_test = 1:100;
h_value = zeros(1,length(s_test));
for i = 1:length(s_test)
    s = s_test(i);
    h_value(i) = h_s(s);
end

%%
H = tf(1,[0.3 1.3 1],'InputDelay', 0); 
hd = c2d(H,0.01,'zoh');
