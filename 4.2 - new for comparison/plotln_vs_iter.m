value_diff2= optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);
value_diff1= optimal_value1(1:optimal_value_mark1) - optimal_value2(optimal_value_mark2);
value_diff3= optimal_value3(1:optimal_value_mark3) - optimal_value2(optimal_value_mark2);
value_diff_1= optimal_value_1(1:optimal_value_mark_1) - optimal_value2(optimal_value_mark2);
value_diff_2= optimal_value_2(1:optimal_value_mark_2) - optimal_value2(optimal_value_mark2);
value_diff_3= optimal_value_3(1:optimal_value_mark_3) - optimal_value2(optimal_value_mark2);
value_diff_4= optimal_value_4(1:optimal_value_mark_4) - optimal_value2(optimal_value_mark2);
figure
semilogx(1:(optimal_value_mark2-3), log(value_diff2(1:optimal_value_mark2-3)));
hold on
semilogx(1:(optimal_value_mark1-1), log(value_diff1(1:optimal_value_mark1-1)));
hold on
semilogx(1:(optimal_value_mark_1), log(value_diff_1(1:optimal_value_mark_1)));
hold on
semilogx(1:(optimal_value_mark_2), log(value_diff_2(1:optimal_value_mark_2)));
hold on
semilogx(1:(optimal_value_mark_3), log(value_diff_3(1:optimal_value_mark_3)));
hold on
semilogx(1:(optimal_value_mark_4), log(value_diff_4(1:optimal_value_mark_4)));
hold on
semilogx(1:(optimal_value_mark3), log(value_diff3(1:optimal_value_mark3)));
hold on
  xlabel('Iterations');
  ylabel('ln(F(x) - F(x^*))');
  title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');
 hold off