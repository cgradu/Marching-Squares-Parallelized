// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include "pthread_barrier_mac.h"
#include <string.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Helper function to find the minimum of two integers
int min(int a, int b){
    return (a > b) ? b : a;
}

// Structure to hold thread data, including image and processing information
typedef struct thread_data {
    int thread_id;
    int step_x;
    int step_y;
    int P;
    unsigned char **grid;
    ppm_image *image;
    ppm_image *new_image;
    ppm_image **contour_map;
    pthread_barrier_t *barrier;
} thread_data;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid(thread_data *data) {
    int step_x = data->step_x;
    int step_y = data->step_y;
    int sigma = SIGMA;
    ppm_image *image = data->new_image;
    unsigned char **grid = data->grid;

    int p = image->x / step_x;
    int q = image->y / step_y;
    int start_i = data->thread_id * (double)p / data->P;
    int end_i = min((data->thread_id + 1) * (double)p / data->P, p);

    for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
        for (int i = start_i; i < end_i; i++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {   
                grid[i][q] = 0;
            } else {
                grid[i][q] = 1;
            }
        }

        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[p][j] = 0;
            } else {
                grid[p][j] = 1;
            }
        }

    return;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(thread_data *data) {
    int step_x = data->step_x;
    int step_y = data->step_y;
    ppm_image **contour_map = data->contour_map;
    unsigned char **grid = data->grid;
    ppm_image *image = data->new_image;

    int p = image->x / step_x;
    int q = image->y / step_y;
    int start_i = data->thread_id * p / data->P;
    int end_i = min((data->thread_id + 1) * p / data->P, p);

    for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

void rescale_image(thread_data *data) {
    ppm_image *image = data->image;
    ppm_image *new_image = data->new_image;
    uint8_t sample[3];

    // we only rescale downwards
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        data->new_image = image;
        return;
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;


	int start_i = data->thread_id * (double) new_image->x / data->P;
    int end_i = min((data->thread_id + 1) * (double) new_image->x / data->P, new_image->x);
    
    // use bicubic interpolation for scaling
    for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }
    return;
}

void *thread_function(void *arg){
    thread_data *data = (thread_data *)arg;
    // 1. Rescale - All threads must finish before any can continue
    rescale_image(data);

    // Wait for all threads to finish resclaing
    pthread_barrier_wait(data->barrier);

    // 2. Sample - Each thread samples a portion of the grid
    sample_grid(data);

    // Wait for all threads to finish sampling
    pthread_barrier_wait(data->barrier);

    // 3. March - Each thread works on a different part of the grid
    march(data);

    // Wait for all threads to finish marching
    pthread_barrier_wait(data->barrier);

    pthread_exit(NULL);

}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    int P = atoi(argv[3]);
    // Arguments for threads
    thread_data *arguments = (thread_data *)malloc( P * sizeof(thread_data));
    if (!arguments) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    // Initialize threads
    pthread_t *threads = (pthread_t *)malloc(P * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    // Alloc memory for image
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;
    // Alloc memory for image data
    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image->data) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    
    // Initialize barrier 
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, P);

    int p = new_image->x / step_x;
    int q = new_image->y / step_y;
    
    // Initialize grid
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // Create threads and give them their arguments
    int ret;
    void *status;
    for(int i = 0; i < P; i++){
    arguments[i].thread_id = i;
    arguments[i].step_x = step_x;
    arguments[i].step_y = step_y;
    arguments[i].P = P;
    arguments[i].grid = grid;
    arguments[i].image = image;
    arguments[i].new_image = new_image;
    arguments[i].contour_map = contour_map;
    arguments[i].barrier = &barrier;

    ret = pthread_create(&threads[i], NULL, thread_function, &arguments[i]);

    if (ret) {
        fprintf(stderr, "Error creating thread %d!\n", i);
        exit(-1);
    }
    }

    // Wait for all threads to finish
    for(int i  = 0; i < P; i++){
        ret = pthread_join(threads[i], &status);
        if (ret) {
            fprintf(stderr, "Error joining thread %d!\n", i);
            exit(-1);
        }
    }

    // Write output
    write_ppm(arguments[0].new_image, argv[2]);

    // Free resources
    free_resources(new_image, contour_map, grid, step_x);
    free(image->data);
    free(image);
    free(arguments);
    free(threads);
    pthread_barrier_destroy(&barrier);

    return 0;
}