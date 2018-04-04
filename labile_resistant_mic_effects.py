from pylab import *

x=linspace(0,0.1,100)
plot(x,x/(x+k))
plot(0.005,0.005/(0.005+k),'bo')
text(0.01,0.005/(0.005+k),'Labile C')

plot(0.1,0.1/(0.1+k),'ro')
text(0.1,0.83,'Labile C',ha='right')

xlabel('Microbial biomass/substrate C ratio')
ylabel('Microbial effect on decomposition rate')
title('Microbial effects')
