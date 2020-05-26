# MSP
A Python 3 model for the creation of two-body simulations with Aerobraking and basic Manoeuvres.

Originally created by Tristan Dijkstra, Andrea Battegazzore and Hugo Chassagnette for the creation of Mars Trajectories for Team Tumbleweed.

![orbit](http://tristandijkstra.com/d/orbit.gif)

## Prerequisites
`python -m pip install --upgrade pip`

`pip install tqdm`

`pip install numpy`

`pip install matplotlib`

`pip install pandas`


For creating matplotlib animations using the animator you need [ffmpeg](https://www.wikihow.com/Install-FFmpeg-on-Windows)

## Install
Currently the package is not available on pypi. To install it directly via github use:

`pip install git+https://github.com/IrTrez/msp.git`

## Usage
To use the model:

`from msp import msp`

To visualise results use the quickAnimate function from simtools or the fully fledged animate function from animator:

`from msp import simtools, animator`

Examples are available [here](https://github.com/IrTrez/msp-examples).
