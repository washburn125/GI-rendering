# Global illumination using photon mapping

This project partially implements using C++ the 2-pass global illumination approach presented in [this article](http://graphics.stanford.edu/~henrik/papers/ewr7/egwr96.pdf).
The scene is rendered using recursive Monte-Carlo path tracing and photon mapping and (for comparison) using basic ray tracing and modulated light attenuation approach that takes into account only direct lightning.
You can find the core of both algorithms in RayTracer.h .

## Project report
Please read the .pdf file with finnal report for more details.

## To compile the project use:

```
mkdir build
cd build
cmake ..
```

## To run the project use:

```
./MyRayTracer -width WIDTH -height HEIGTH -output filename.pmm
```
where WIDTH and HEIGHT correspond to desired width and height of resulting image that will be saved as filename.pmm file.
