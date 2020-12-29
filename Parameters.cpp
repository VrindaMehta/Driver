#include "basic.h"
#include <iostream>
#include <iomanip>
#include <string.h>

void params() {

    int i , j, l;

    for(i = 0; i < n; i++) {   //To add the driver
      h_x[i] = 1;
      for(j = 0; j < n; j++) {
       if(flag[j + i * n] != 0) {
          J_x[j + i * n] = 1;
        }
        else{
          J_x[j + i*n] = 0;
        }
      }
    } 
    for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
        for (l = 0; l < n; l++){
          K_x[l + j*n + i*n*n] = 0;
        }
      }
    }
    for (i = 0; i+2 < n; i+=2){
      K_x[i+2 + (i+1)*n + (i)*n*n] = 1;
      
    }

    
}
