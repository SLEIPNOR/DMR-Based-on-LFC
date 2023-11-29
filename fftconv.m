function Conv2 = fftconv(a,b)



Conv2_size = size(a);

a_FT = fft2(a,Conv2_size(1),Conv2_size(2));
b_FT = fft2(b,Conv2_size(1),Conv2_size(2));

Conv2_FT = a_FT.*b_FT;
Conv2 = ifft2(Conv2_FT);
