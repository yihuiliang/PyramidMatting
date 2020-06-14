# Pyramid matting
Pyramid matting is a resource-adaptive multi-scale pixel pair optimization framework for image matting

If you use this source code please cite:  
```Yihui, Liang, Feng Fujian, and Cai Zhaoquan. "Pyramid Matting: A Resource-Adaptive Multi-Scale Pixel Pair Optimization Framework for Image Matting." IEEE Access 8 (2020): 93487-93498.```

# Dependance
The demo uses an evolutionary algorithm named competitive swarm optimization, which is proposed by Prof. Yaochu Jin and Prof. Ran Cheng.
```Cheng, Ran, and Yaochu Jin. "A competitive swarm optimizer for large scale optimization." IEEE transactions on cybernetics 45.2 (2014): 191-204.```

Please notice that pyramid matting framework can work with any other general evolutionary algorithms.
# Example usage:
```matlab
cd path     %go to the path of the source files
PyramidMattingDemo_CSO %run the demo
```