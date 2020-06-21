//
//  serialForest.cpp
//  Forest
//
//  Created by Matteo Muraca
//  Copyright © 2020 Matteo Muraca. All rights reserved.
//
//  Forest - a green cellular automata
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


#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
using namespace std;

//Game settings
const int GROUND = 0;
const int TREE = 1;
const int BURNING = -1;
//matrix dimension
const int dim = 300;
//lightining probability
const int l = 2000;
//growth probability
const int g = 200;
//seed for random
unsigned int seed = 31415; //I liked π

//Allegro settings
const int width = 1500;
const int height = 1500;
ALLEGRO_DISPLAY * display;

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

int computeTree(int ** cells, int x, int y) {
    
    if(cells[x][y]==GROUND){ //if the cell is ground, a tree can grow with p = 1/g
        int randomN = rand_r(&seed) % g;
        if(randomN==0)
            return TREE;
    }
    else if(cells[x][y]==TREE){
        //if a tree next to the current one is burning, the tree will be burning
        if(x>0){
            if(cells[x-1][y]==BURNING)
                return BURNING;
        }
        if(x<dim - 1){
            if(cells[x+1][y]==BURNING)
                return BURNING;
        }
        if(y>0){
            if(cells[x][y-1]==BURNING)
                return BURNING;
        }
        if(y<dim-1){
            if(cells[x][y+1]==BURNING)
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
    
    int ** cells1 = matrixAllocation<int>(dim, dim);
    int ** cells2 = matrixAllocation<int>(dim, dim);
    initModel(cells1);
    
    al_init();
    al_init_primitives_addon();
    display = al_create_display(width, height);
    al_install_keyboard();
    drawCells(cells1);
    
    ALLEGRO_KEYBOARD_STATE key_state;
    al_get_keyboard_state(&key_state);
    
    bool whatMatrix = true; //true: compute cells2 from cells1, false: compute cells1 from cells2
    
    while(!al_key_down(&key_state, ALLEGRO_KEY_ESCAPE)) {
        
        if(whatMatrix) {
            for(int i=0; i<dim; i++) {
                for(int j=0; j<dim; j++) {
                    //compute the new matrix
                    cells2[i][j] = computeTree(cells1, i, j);
                }
            }
            drawCells(cells2);
        }
        else {
            for(int i=0; i<dim; i++) {
                for(int j=0; j<dim; j++) {
                    //compute the new matrix
                    cells1[i][j] = computeTree(cells2, i, j);
                }
            }
            drawCells(cells1);
        }
        
        whatMatrix = !whatMatrix; //swap
        
        al_get_keyboard_state(&key_state);
    }
    
    deleteMatrix(cells1);
    deleteMatrix(cells2);
    
    return 0;
}


