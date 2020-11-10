function [MSE, MAE] = Indicator(x_clean, x_denoised)


x_clean = x_clean(:);
x_denoised = x_denoised(:);

MSE = mean((x_clean - x_denoised).^2);
MAE = mean(abs(x_clean - x_denoised));

end

