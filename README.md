# Sparse partial inverses in Python. 

Code is a wrap of the C++ provided in this blog-post https://dansblog.netlify.app/posts/2024-09-05-partial-inverse/partial-inverse (https://github.com/dpsimpson/blog/tree/master/posts/2024-09-05-partial-inverse). Read https://www.sciencedirect.com/science/article/pii/S0378375807000845 to learn more about the technique.

To use in your Python project, include the ```partial_inverse.cpython-312-darwin.so``` file and import the module by ```import partial_inverse```. The only function is the ```pinv``` function, which takes in a ```scipy.sparse``` CSC-matrix and outputs its inverse. Make sure the matrix is symmetric positive definite!
