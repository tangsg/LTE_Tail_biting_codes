/*   Tailbiting Encoder
 * 	Written By Dama Sreekanth
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers/tail_enc.h"

void tail_enc( int MIB[],int len,float *C0,float *C1,float *C2)
{

	int mem[7];
	mem[5]=MIB[len-6];mem[4]=MIB[len-5];mem[3]=MIB[len-4];mem[2]=MIB[len-3];mem[1]=MIB[len-2];mem[0]=MIB[len-1];mem[6]=0; //initialize the memory
	int jj=0;
	for (jj=0;jj<len;jj++){
		mem[6]=mem[5];mem[5]=mem[4];mem[4]=mem[3];mem[3]=mem[2];mem[2]=mem[1];mem[1]=mem[0];
		mem[0]=MIB[jj];
		*(C0+jj) = 1/sqrt(2)*(2*(mem[0]^mem[2]^mem[3]^mem[5]^mem[6])-1);
		*(C1+jj) = 1/sqrt(2)*(2*(mem[0]^mem[1]^mem[2]^mem[3]^mem[6])-1);
		*(C2+jj) = 1/sqrt(2)*(2*(mem[0]^mem[1]^mem[2]^mem[4]^mem[6])-1);
	}
}
