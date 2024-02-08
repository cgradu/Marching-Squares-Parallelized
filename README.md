# Marching Squares Parallelized
## Author: Constantin Radu-Giorgian, 331CC

The "tema1_par.c" file contains a program designed to detect contours in images using the Marching Squares algorithm parallelized with threads.

### Structures and Definitions

- `thread_data`: Contains all necessary data for a thread to process a portion of the image.
- `#define` statements: Define various constants used in the program, such as scaling dimensions, steps, and limits.

### Main Functions

- `rescale_image()`: Resizes the image to specified dimensions using bicubic interpolation to maintain quality.
- `sample_grid()`: Samples a grid from the image, determining the value of each point (0 or 1) based on comparisons with a threshold (sigma).
- `march()`: Identifies the contour type for each subgrid of the sampled grid and replaces pixels in the original image with those of the corresponding contour image.
- `thread_function()`: Executed by each thread, including scaling, sampling, and contour detection steps. Uses `pthread_barrier_wait` for thread synchronization.
- `main()`: Checks input arguments and initializes necessary structures and variables. Creates a specified number of threads (P), each processing a portion of the image. After all threads finish, the processed image is written to the output file. Releases resources and terminates execution.

### Auxiliary Functions

- `init_contour_map()`: Initializes a contour map by reading predefined images representing different contour configurations.
- `update_image()`: Updates a section of an image with corresponding pixels from the contour image.
- `free_resources()`: Frees memory allocated for contour images, grid, and processed images.

### Parallelization Process

To adapt the Marching Squares algorithm for parallelization, functions for scaling, sampling, and contour detection were modified. Start and end variables are used to set the limit of each image portion processed by a thread.

In the rescaling function, start and end are calculated based on the new image dimensions, ensuring each thread processes a distinct image section, avoiding overlap and errors.

For sampling and contour detection, start and end are set using the step variable to evenly distribute the workload among threads.

The thread function runs these modified functions. A barrier is included between execution stages to synchronize threads and prevent conflicts in accessing common resources.

Other functions, not mentioned here, remain unchanged from the sequential version of the Marching Squares algorithm as presented in the "tema1.c" file. They retain their roles and functionality in the parallelized program version.