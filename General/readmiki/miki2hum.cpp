

//include standard libraries
#include <sstream>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h> 
#include <fstream>
#include <vector>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]){
    
  //get lenght 
  int lenght=argc;
  int string;
  

  //loop and print human
  for(int i=2;i<lenght;i++){
    string=atoi(argv[i]);
    putchar(string);
  }

  exit(0);
}
