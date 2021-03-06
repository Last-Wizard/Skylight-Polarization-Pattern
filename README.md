Skylight-Polarization-Pattern
=============================
>version: 1.0

>last_update_time: 2014-10-20

>reference:

>1. Pomozi I, Gál J, Horváth G, et al. Fine structure of the celestial polarization pattern and its temporal change during the total solar eclipse of 11 August 1999\[J\]. Remote sensing of Environment, 2001, 76(2): 181-201.
>2. 吴良海, 高隽, 范之国, 等. 基于复球面映射的大气偏振模式表征与分析\[J\]. 仪器仪表学报, 2011, 32(4): 870-876. (in Chinese)
>3. Radiative transfer theory [link](http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/the_vector_radiative_transfer_equation)
>4. Cornet C. Three-dimensional polarized Monte Carlo atmospheric radiative transfer model (3DMCPOL): 3D effects on polarized visible reflectances of a cirrus cloud[J]. Journal of Quantitative Spectroscopy and Radiative Transfer, 2010, 111(1): 174-186.

###理想Rayleigh单次散射条件下的大气偏振模式建模与仿真

>skylight polarization pattern based on Rayleigh model

>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img1.jpg)

>skylight polarization pattern based on vector radiative transfer equation [link](http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/the_vector_radiative_transfer_equation)

>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img2.jpg)


>偏振度dop: degree of polarization

>偏振化方向角aop: angle of polarization

### Rayleigh 3D
>三维空间上的大气偏振模式分布 (3-dimensional), 需要使用mayavi库 (mayavi library)

>直接投影

>复平面投影

### Rayleigh 2D
>二维平面中的大气偏振模式分布 (2-dimensional)

>直接投影 (Rayleigh_2D_directly.py)

>>偏振度 dop (太阳高度角30°, 方位角45°)

>>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img3.png)

>>偏振化方向角 aop (太阳高度角30°, 方位角45°)

>>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img4.png)

>复平面投影 (Rayleigh_2D_complex.py)

>>偏振度 aop (太阳高度角30°, 方位角45°)

>>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img5.png)

>>偏振化方向角 aop (太阳高度角30°, 方位角45°)

>>![img](https://github.com/ConanGit/gallery/blob/master/Skylight-Polarization-Pattern/img6.png)

###Environment
>python 2.7.8, numpy 1.8.1, matplotlib 1.3.1, spyder 2.3.0

