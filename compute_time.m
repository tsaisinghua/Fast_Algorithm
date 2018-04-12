a = 1.234; b=2.345;
tic;
for i=1:100000000
    a = a+b;
    b = a-b;
    b = a-b;
end
T1 = toc;
%fprintf('(+,-) time:%f\n',T1);
tic;
for i=1:100000000
    a = a+b;
    b = a-b;
    b = a-b;
    a = a+b;
    b = a-b;
    b = a-b;
end
T2 = toc;
%fprintf('(+,-) time:%f\n',T2);
fprintf('(+,-) time:%f\n',(T2-T1)/3.0);

tic;
for i=1:100000000
    a = a*b;
    b = a/b;
    b = a/b;
end
T1 = toc;
%fprintf('(*,/) time:%f\n',T1);
tic;
for i=1:100000000
    a = a*b;
    b = a/b;
    b = a/b;
    a = a*b;
    b = a/b;
    b = a/b;

end
T2 = toc;
%fprintf('(*,/) time:%f\n',T2);
fprintf('(*,/) time:%f\n',(T2-T1)/3.0);
tic;
for i=1:100000000
    a = sin(a);
end
T1 = toc;
tic;
for i=1:100000000
    a = sin(b);
    b = sin(a);
end
T2 = toc;
fprintf('(sin) time:%f\n', T2-T1);

tic;
for i=1:10000
    a = randi(32768,100);
end
fprintf('randi time:%f\n',toc);