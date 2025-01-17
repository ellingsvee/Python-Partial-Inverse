# Sparse partial inverses in Python. 

Code is a wrap of the C++ provided in this blog-post https://dansblog.netlify.app/posts/2024-09-05-partial-inverse/partial-inverse (https://github.com/dpsimpson/blog/tree/master/posts/2024-09-05-partial-inverse). Read https://www.sciencedirect.com/science/article/pii/S0378375807000845 to learn more about the technique.

To use the function in your Python project, build the ```partial_inverse_module.cpp``` file using ```pybind11```. Place the ```partial_inverse.cpython-...-darwin.so``` file in the same directory as the ```partial_inverse.py``` file. The only function is the ```pinv``` function, which takes in a ```scipy.sparse``` CSC-matrix and outputs its inverse. Make sure the matrix is positive definite!

If you want to modify the code, you can change the ```partial_inverse_module.cpp``` file and wrap it in the standard way using ```pybind11```.
