function [t_ut, t_ut2] = test(t_inn)
  t_ut = t_inn;
  t_inn = t_inn + 3;
  t_ut2 = t_inn;
end

