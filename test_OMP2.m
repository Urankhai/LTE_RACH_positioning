% x = OMP( A, b, k )

clear all
close all

A = randn(839,100);
k = 5;
x_vect = zeros(100,1);
x_supp = sort(randi(100,[1,k]));
x_vect(x_supp) = 100*rand(k,1);
disp(['x org = ',num2str(nonzeros(x_vect)')])

b_org = A*x_vect;
b = b_org;% + 10*randn(839,1);
figure
hold on
plot(b)
plot(b_org)

x = Mathwork_OMP( A, b, k );
disp(['x est = ',num2str(nonzeros(x)')])