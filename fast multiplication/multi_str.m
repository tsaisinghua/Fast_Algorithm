function s = multi_str(s1,s2)
% s= s1*s2 in string form
N = length(s1)+length(s2);
v = zeros(1,N);
v1 = zeros(1,N);
v2 = zeros(1,N);
v1(1:length(s1)) = fliplr(double(s1)-48);  % fliplr:低位元在前面 
v2(1:length(s2)) = fliplr(double(s2)-48);  % double(s)-48:字串改數字 
% Convolution of v1 and v2
for i=1:N
    for j=1:i
        v(i) = v(i)+v1(j)*v2(i-j+1);
    end
end
% 進位
for i=1:N-1
    v(i+1) = v(i+1)+floor(v(i)/10);
    v(i) = mod(v(i),10);
end
% 砍掉0
v = fliplr(v);
for i=1:N
    if v(i) ~= 0
        break;
    end
end
s = char(v(i:end)+48);