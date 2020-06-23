//
//  parallelForest.cpp
//  Forest
//
//  Created by Matteo Muraca
//  Copyright © 2020 Matteo Muraca. All rights reserved.
//
//  Parallel version of Forest - a green cellular automata
//  Programmed in C++ using MPI libraries
//
//  Forest rules:
//  1.  A burning tree turns into an empty cell.
//  2.  A non-burning tree with one burning neighbour turns into a burning tree.
//  3.  A tree ignites with probability 1/f due to lightning.
//  4.  An empty space grows a new tree with probability 1/p.
//
//  A cell can be:
//  1. Tree
//  0. Ground
// -1. Burning tree
//
//  You should use a divider of dim as number of processes.


#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <mpi.h>
using namespace std;

//Game settings
const int GROUND = 0;
const int TREE = 1;
const int BURNING = -1;
//matrix dimension
const int dim = 200;
//lightining probability
const int l = 2000;
//growth probability
const int g = 200;
//seed for random
unsigned int seed;

//Allegro settings
const int width = 1500;
const int height = 1500;
ALLEGRO_DISPLAY * display;

//parallel processing variables
const int root = 0;
int numOfProcesses;

//fill vector with value
void fillVector(int* vector, int dimension, int value) {
    for(int i=0; i<dimension; i++)
        vector[i] = value;
}

//2d array (matrix) contiguous allocation
template <class T>
T** matrixAllocation(int N, int M) {
    T* tmp = (T*) malloc(sizeof(T)*N*M);
    T** mat = (T**) malloc(sizeof(T*)*N);
    for(int i = 0 ; i < N ; i++) {
        mat[i] = &(tmp[M*i]);
    }
    return mat;
}

//2d array (matrix) delete
template <class T>
void deleteMatrix(T** mat) {
    if(mat == nullptr)
        return;
    free(mat[0]);
    free(mat);
}

//initialize the main matrix with non-burning cells (tree or ground)
void initModel(int** cells) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            cells[i][j] = rand_r(&seed) % 2;
        }
    }
}

//draw the matrix
void drawCells(int ** cells) {
    al_clear_to_color(al_map_rgb(0,0,0));
    
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            switch (cells[i][j]) {
                case GROUND:
                    al_draw_filled_rectangle(j * width/dim, i * width/dim,
                                             j * width/dim + width/dim,
                                             i * height/dim + height/dim, al_map_rgb(73,32,0));
                    break;
                case TREE:
                    al_draw_filled_rectangle(j * width/dim, i * width/dim,
                                             j * width/dim + width/dim,
                                             i * height/dim + height/dim, al_map_rgb(0, 150, 0));
                    
                    break;
                case BURNING:
                    al_draw_filled_rectangle(j * width/dim, i * width/dim,
                                             j * width/dim + width/dim,
                                             i * height/dim + height/dim, al_map_rgb(200, 0, 0));
                    break;
                default:
                    al_draw_filled_rectangle(j * width/dim, i * width/dim,
                                             j * width/dim + width/dim,
                                             i * height/dim + height/dim, al_map_rgb(255, 255, 255));
                    break;
            }
            
        }
        
    }
    al_flip_display();
    al_rest(0.5);
}


//if you press esc, the game will finish
bool continueProcessing(int rank) {
    int buf = 0;
    if(rank==root) {
        ALLEGRO_KEYBOARD_STATE key_state;
        al_get_keyboard_state(&key_state);
        if(!al_key_down(&key_state, ALLEGRO_KEY_ESCAPE)) //root will check if the key esc is pressed
            buf = INT_MAX; //if so, this variable's value will become INT_MAX
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&buf, 1, MPI_INT, root, MPI_COMM_WORLD); //the variable will be broadcasted to all processors
    return buf == INT_MAX; //if it's INT_MAX, the program will be stopped
}

int computeTree(int ** subMatrix, int x, int y, int * upperVector, int * lowerVector) {
    
    if(subMatrix[x][y]==GROUND){ //if the cell is ground, a tree can grow with p = 1/g
        int randomN = rand_r(&seed) % g;
        if(randomN==0)
            return TREE;
    }
    else if(subMatrix[x][y]==TREE){
        //if a tree next to the current one is burning, the tree will be burning
        if(x>0){
            if(subMatrix[x-1][y]==BURNING)
                return BURNING;
        }
        else if(upperVector[y]==BURNING)
            return BURNING;
        
        if(x<dim/numOfProcesses - 1){
            if(subMatrix[x+1][y]==BURNING)
                return BURNING;
        }
        else if(lowerVector[y]==BURNING){
            return BURNING;
        }
        if(y>0){
            if(subMatrix[x][y-1]==BURNING)
                return BURNING;
        }
        if(y<dim-1){
            if(subMatrix[x][y+1]==BURNING)
                return BURNING;
        }
        int randomN = rand_r(&seed) % l; //if the tree is not already burning, it ignites with p = 1/l
        if(randomN==0)
            return BURNING;
        
        return TREE;
    }
    //if the cell was not a tree nor ground where a tree has grown
    //it will still be ground, or become ground if there was a burning tree
    return GROUND;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    
    if(dim%numOfProcesses!=0 && rank==root) {
        cout<<"Error: please use a divider of "<< dim <<" as number of processes."<<endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    
    seed = 31415 * (rank+1); //every process will have a different seed
    
    int ** cells = nullptr;
    
    if(rank==root) {
        cells = matrixAllocation<int>(dim, dim);
        initModel(cells);
        
        al_init();
        al_init_primitives_addon();
        display = al_create_display(width, height);
        al_install_keyboard();
        drawCells(cells);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int ** subMatrix1 = matrixAllocation<int>(dim/numOfProcesses, dim);
    int ** subMatrix2 = matrixAllocation<int>(dim/numOfProcesses, dim);
    
    //split the matrix into submatrices, one for each processor
    if(rank==root)
        MPI_Scatter(&cells[0][0], dim*dim/numOfProcesses, MPI_INT, &subMatrix1[0][0], dim*dim/numOfProcesses, MPI_INT, root, MPI_COMM_WORLD);
    else
        MPI_Scatter(NULL, 0, NULL, &subMatrix1[0][0], dim*dim/numOfProcesses, MPI_INT, root, MPI_COMM_WORLD);
    
    
    int * upperVector = (int*) malloc(sizeof(int)*dim);
    int * lowerVector = (int*) malloc(sizeof(int)*dim);
    
    bool whatMatrix = true; //true: compute subMatrix2 from subMatrix1, false: compute subMatrix1 from subMatrix2
    
    while(continueProcessing(rank)) {
        MPI_Request r;
        MPI_Status s;
        
        if(whatMatrix) {
            
            if(rank!=root) //send the lowerVector to the previous process
                MPI_Isend(&subMatrix1[0][0], dim, MPI_INT, rank-1, 11, MPI_COMM_WORLD, &r);
            
            if(rank!=numOfProcesses-1) //recieve the lowerVector, or fill it with 0
                MPI_Recv(&lowerVector[0], dim, MPI_INT, rank+1, 11, MPI_COMM_WORLD, &s);
            else
                fillVector(lowerVector, dim, 0);
            
            if(rank!=numOfProcesses-1) //send the upperVector to the next process
                MPI_Isend(&subMatrix1[dim/numOfProcesses - 1][0], dim, MPI_INT, rank+1, 44, MPI_COMM_WORLD, &r);
            
            if(rank!=root) //recieve the upperVector, or fill it with 0
                MPI_Recv(&upperVector[0], dim, MPI_INT, rank-1, 44, MPI_COMM_WORLD, &s);
            else
                fillVector(lowerVector, dim, 0);
            
            for(int i=0; i<dim/numOfProcesses; i++) {
                for(int j=0; j<dim; j++) {
                    //compute the new submatrix for each processor
                    subMatrix2[i][j] = computeTree(subMatrix1, i, j, upperVector, lowerVector);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            //the new matrix is stored into the root, overwriting the old one
            if(rank==root)
                MPI_Gather(&subMatrix2[0][0], dim*dim/numOfProcesses, MPI_INT, &cells[0][0], dim*dim/numOfProcesses, MPI_INT, root, MPI_COMM_WORLD);
            else
                MPI_Gather(&subMatrix2[0][0], dim*dim/numOfProcesses, MPI_INT, NULL, 0, MPI_INT, root, MPI_COMM_WORLD);
            
        }
        else {
            
            if(rank!=root) //send the lowerVector to the previous process
                MPI_Isend(&subMatrix2[0][0], dim, MPI_INT, rank-1, 22, MPI_COMM_WORLD, &r);
            
            if(rank!=numOfProcesses-1) //recieve the lowerVector, or fill it with 0
                MPI_Recv(&lowerVector[0], dim, MPI_INT, rank+1, 22, MPI_COMM_WORLD, &s);
            else
                fillVector(lowerVector, dim, 0);
            
            if(rank!=numOfProcesses-1) //send the upperVector to the next process
                MPI_Isend(&subMatrix2[dim/numOfProcesses - 1][0], dim, MPI_INT, rank+1, 33, MPI_COMM_WORLD, &r);
            
            if(rank!=root) //recieve the upperVector, or fill it with 0
                MPI_Recv(&upperVector[0], dim, MPI_INT, rank-1, 33, MPI_COMM_WORLD, &s);
            else
                fillVector(lowerVector, dim, 0);
            
            for(int i=0; i<dim/numOfProcesses; i++) {
                for(int j=0; j<dim; j++) {
                    //compute the new submatrix for each processor
                    subMatrix1[i][j] = computeTree(subMatrix2, i, j, upperVector, lowerVector);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            //the new matrix is stored into the root, overwriting the old one
            if(rank==root)
                MPI_Gather(&subMatrix1[0][0], dim*dim/numOfProcesses, MPI_INT, &cells[0][0], dim*dim/numOfProcesses, MPI_INT, root, MPI_COMM_WORLD);
            else
                MPI_Gather(&subMatrix1[0][0], dim*dim/numOfProcesses, MPI_INT, NULL, 0, MPI_INT, root, MPI_COMM_WORLD);
            
        }
        
        whatMatrix = !whatMatrix;
        
        if(rank==root)
            drawCells(cells);
    
    }
    
    deleteMatrix(cells);
    deleteMatrix(subMatrix1);
    deleteMatrix(subMatrix2);
    free(upperVector);
    free(lowerVector);
    
    MPI_Finalize();
    return 0;
}

