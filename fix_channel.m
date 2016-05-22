function [ c] = fix_channel(w_t_1,w_t_2,X1,g)
D=8;%width of channel
w_t_1=w_t_1;
w_t_2=w_t_2;
X1=X1;
g=g;
teta_1_founded=atan(-w_t_1(2,1)/w_t_1(3,1));
% teta_1_founded_in_degree=(teta_1_founded*180)/pi
teta_2_founded=atan(-w_t_2(2,1)/w_t_2(3,1));
% teta_2_founded_in_degree=(teta_2_founded*180)/pi
mean_teta_1_and_teta_2=(teta_1_founded+teta_2_founded)/2;
% mean_teta_1_and_teta_2_in_degree=(mean_teta_1_and_teta_2*180)/pi
w_t_1_new(3,1)=w_t_1(3,1);
w_t_1_new(2,1)=-tan(mean_teta_1_and_teta_2)*w_t_1_new(3,1);
w_t_2_new(3,1)=w_t_2(3,1);
w_t_2_new(2,1)=-tan(mean_teta_1_and_teta_2)*w_t_2_new(3,1);
if (-w_t_1(1,1)/w_t_1(3,1))>= (-w_t_2(1,1)/w_t_2(3,1))
    x_old=(-w_t_1(1,1)/w_t_1(3,1))-(-w_t_2(1,1)/w_t_2(3,1));
else
    x_old=(-w_t_2(1,1)/w_t_2(3,1))-(-w_t_1(1,1)/w_t_1(3,1));
end
x_new=abs(D/cos(mean_teta_1_and_teta_2));
e=x_new-x_old;
if (-w_t_1(1,1)/w_t_1(3,1))>=(-w_t_2(1,1)/w_t_2(3,1))
    w_t_1_new(1,1)=-(-w_t_1(1,1)+w_t_1_new(3,1)*(e/2));
    w_t_2_new(1,1)=-(-w_t_2(1,1)-w_t_2_new(3,1)*(e/2));
else
    w_t_1_new(1,1)=-(-w_t_1(1,1)-w_t_1_new(3,1)*(e/2));
    w_t_2_new(1,1)=-(-w_t_2(1,1)+w_t_2_new(3,1)*(e/2));
end
w_t_1=w_t_1_new;
w_t_2=w_t_2_new;
c=g(w_t_1'*X1).*g(w_t_2'*X1);
% c=c>.5;
% c=g(w(1:3)'*X1).*g(w(4:6)'*X1)
end

