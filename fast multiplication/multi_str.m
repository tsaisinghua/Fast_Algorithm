function s = multi_str(s1,s2)

N = length(s1)+length(s2);
v = zeros(1,N);
v1 = zeros(1,N);
v2 = zeros(1,N);
v1(1:length(s1)) = fliplr(double(s1)-48);
v2(1:length(s2)) = fliplr(double(s2)-48);
% s='123': s=123 => double(s): 49 50 51(ASKII code, 0=>48;1=>49;...) => double(s)-48【字串改數字〕: 1 2 3
% fliplr(翻轉位數：低位原擺前面) => ex.double('0123')= 48 49 50 51 => fliplr(double('0123'))= 51 50 49 48

% convolution of x1 and v2  => 若要快速乘法，需修改此段。
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
% 砍掉 0 
v = fliplr(v);
for i=1:N
    if v(i) ~= 0
        break;
    end
end
s = char(v(i:end)+48);