function [f_out] = diff2(f_in,ableiten)
syms substitute;
f_zwischen=subs(f_in,ableiten,substitute);
f_zwischen_diff=diff(f_zwischen,substitute);
f_out=subs(f_zwischen_diff,substitute,ableiten);
end