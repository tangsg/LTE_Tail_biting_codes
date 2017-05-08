/*
Convolutional tail biting encoder decoder test
    Copyright (C) 2017  sreekanth
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * Author: sreekanth dama
 * Contact: sreekanth@iith.ac.in
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define length1 25
void tail_enc( int MIB[],int len,float *C0,float *C1,float *C2);
void viterbi_wava(float p_zerobits[],float p_onebits[],float p_twobits[],unsigned int length,int *outputBitStream);

int main()
{
    printf("Hello, World!\n");
    
    int bits[length1]={1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1};
    float C0[length1],C1[length1],C2[length1];
    int out[length1];
    
    tail_enc(bits,length1,C0,C1,C2);
    viterbi_wava(C0,C1,C2,length1,out);
    int ii=0;
    for (ii=0; ii<length1; ii++){
		//fprintf(f1,"%c \n",*(outputBitStream+ii));
		printf("%d",*(out+ii));
	}
    return 0;
}
