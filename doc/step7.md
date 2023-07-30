# Step 7:

a. Coordinates are received from the GUI as input, including the center and the two box sizes

b. Images are generated. They are simply a cropped version of the input image around $V_1$ along each direction.

c. Between the coordinates that are inputted from Step a, iterate through each box and sum each pixel's intensity. While adding each pixel's intensity, cound the number of pixels iterated to find a volume. Do this for both $V_1$ and $V_2$.

d. Calculate $V_0$ with $\frac{S_{1,SE} \cdot V_{2,SE} - S_{2,SE} \cdot V_{1,SE}}{S_{1,SE} - S_{2,SE}}$

e. Calculate $a$ (object radius) with $(\frac{3V_0}{4\pi})^{1/3}$

f. Calculate $\rho_{0,SE}$ with $\frac{S_{2,SE}}{V_{2,SE} - V_0}$ or $\frac{S_{1,SE} - S_{2,SE}}{V_{1,SE} - V_{2,SE}}$

g. Calculate ${SNR}_{SE}$ with $\frac{\rho_{0,SE}}{\sigma_{SE}}$. $\sigma_{SE}$ (or dsnr in code) is the input standard deviation on the GUI.

h. Calculate the uncertainty of $V_0$ ($\delta V_0$) with $\frac{\sqrt{\delta V}}{{SNR}_{SE}} \cdot \sqrt{V_{2,SE} + \frac{(V_{2,SE} - V_0)^2}{V_{1,SE} - V_{2,SE}}}$. $\delta V$ is set to be 1 as a constant.

i. Calculate the uncertainty of the radius $a$ ($\delta a$) with $a \cdot \frac{\delta V_0}{3V_0}$

j. Calculate the susceptibility ($\Delta \chi$ in ppm) with $2p / (\gamma_{\bar{R}} \cdot B_0 \cdot {TE}_{last} \cdot V_0)$. Some of these values are from Step 5 in the GUI. $\gamma_{\bar{R}}$ is a constant of 42.58.

k. Calculate the uncertainty of the susceptibility ($\delta \Delta \chi$) with $\delta \chi \cdot \sqrt{ (\frac{\delta p}{p})^2 + (\frac{\delta V_0}{V_0})^2}$

