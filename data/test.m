load wave_[2018_4_17-14_12_52].csv
p=fft(wave__2018_4_17_14_12_52_(:,2));
p(2);
plot(abs(p(2:400)));