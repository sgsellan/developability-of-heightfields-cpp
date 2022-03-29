# Developability of Heightfields via Rank Minimization

:warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning:  ***DISCLAIMER PLEASE READ THIS*** :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning:

**This is a bad, unoptimized, restricted implementation of the core functionality of the SIGGRAPH 2020 ["Developability of Heightfields via Rank Minimization"](http://dgp.toronto.edu/~sgsellan/pdf/compressed-developables.pdf) by [Silvia Sellán](http://dgp.toronto.edu/~sgsellan/), [Noam Aigerman](https://noamaig.github.io/) and [Alec Jacobson](http://www.cs.toronto.edu/~jacobson/).** **The proper code release for the paper, including scripts to recreate most paper figures, is [here](https://github.com/sgsellan/developability-of-heightfields)** **<- That hyperlink is most likely what you are looking for.**

:warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: 

Still, if running our algorithm would help with your own research or you want to be able to compare to our algorithm and you don't have access to Matlab, I provide this C++ implementation **without any guarantees at all**. I wrote this after the paper's publication as an exercise in learning C++ and may occassionally add functionality. Please do not hesitate to contact [sgsellan@cs.toronto.edu](mailto:sgsellan@cs.toronto.edu) if you find any issues or bugs in this code, or you need any additional functionality.


Please also note that while this C++ implementation of our method is hereby released under MIT License, the method itself is pending a US patent filed in 2020 by Adobe Inc.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ../ 
    make

This will create a `main` binary. A pre-compiled MacOS binary is also provided here.

## Run

From within the main directory just issue:

    ./main

A viewer app should launch showing a 3D mesh of a mountain range, and instructions should appear on the console.

![](img/range.png)

You can begin by hitting H in your keyboard to convert this 3D mesh into a heightfield from above, which will appear on screen.

![](img/range-h.png)

Now, you can further hit S to get the closest developable heightfield (this can take around a minute). 

![](img/range-d.png)

You are not limited to our sample. In general, you can run our code as

    ./main -i input_shape.obj -n grid_size -o omega

and follow the same instructions.



