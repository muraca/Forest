//
//  forest.c
//  Forest
//
//  Created by Matteo Muraca
//  Copyright © 2020 Matteo Muraca. All rights reserved.
//


/*
 Forest rules:
 1.  A burning tree turns into an empty cell.
 2.  A non-burning tree with one burning neighbour turns into a burning tree.
 3.  A tree with no burning neighbour ignites with probability 1/f due to lightning.
 4.  An empty space grows a new tree with probability 1/p.
 
 A cell can be:
 1. Tree
 0. Ground
 -1. Burning tree
 
 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

const int lines = 10;
const int columns = 10;
const int GROUND = 0;
const int TREE = 1;
const int BURNING = -1;
const int f = 2000;
const int p = 200;

int computeTree(int** cells, int x, int y, int xSize, int ySize, unsigned int* seed) {
    if(cells[x][y]==GROUND){
        int random = rand_r(seed) % p;
        if(random==0)
            return TREE;
    }
    else if(cells[x][y]==TREE){
        if(x>0){
            if(cells[x-1][y]==BURNING)
                return BURNING;
        }
        if(x<xSize-1){
            if(cells[x+1][y]==BURNING)
                return BURNING;
        }
        if(y>0){
            if(cells[x][y-1]==BURNING)
                return BURNING;
        }
        if(y<ySize-1){
            if(cells[x][y+1]==BURNING)
                return BURNING;
        }
        int random = rand_r(seed) % f;
        if(random==0)
            return BURNING;
        return TREE;
    }
    return GROUND;
}

void printCell(int** cells, int x, int y){
    if(cells[x][y]==GROUND)
        printf("G");
    else if(cells[x][y]==TREE)
        printf("T");
    else
        printf("B");
}

int mainForest(int argc, const char * argv[]) {
    unsigned int *seed = malloc(sizeof(int)); *seed = 50;
    
    int** cells = (int **)malloc(lines * sizeof(int *));
    for (int i = 0; i < lines; ++i) {
        cells[i] = (int *)malloc(columns * sizeof(int));
        for (int j = 0; j < columns; ++j) {
            cells[i][j] = rand_r(seed) % 2;
        }
    }
    
    while(TREE==1){
        for (int i = 0; i < lines; ++i) {
            for (int j = 0; j < columns; ++j) {
                printCell(cells, i, j);
            }
            printf("\n");
        }
        
        printf("Computing new matrix\n");
        int** newCells = (int **)malloc(lines * sizeof(int *));
        for (int i = 0; i < lines; ++i) {
            newCells[i] = (int *)malloc(columns * sizeof(int));
            for (int j = 0; j < columns; ++j) {
                newCells[i][j] = computeTree(cells, i, j, lines, columns, seed);
            }
        }
        printf("Done! Trashing old matrix...\n");
        for(int i=0; i<lines; ++i)
            free(cells[i]);
        free(cells);
        cells = newCells;
        sleep(1);
    }
    
}
