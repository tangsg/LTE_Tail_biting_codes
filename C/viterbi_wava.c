/*
Wrap Around Viterbi Algorithm
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



void viterbi_wava(float p_zerobits[],float p_onebits[],float p_twobits[],unsigned int length,int *outputBitStream)
{
	int location=0;int currentstate=0;
	int jj=0;int ii=0;int kk=0;
	int LUT[64][64] ={[0 ... 63][0 ... 31]=0,[0 ...63][32 ... 63]=1}; //lut for traceback unit
			//codewords
	int code[128][3]={{0,0,0},{1,1,1},{1,0,0},{0,1,1},{0,0,1},{1,1,0},{1,0,1},{0,1,0}};
	float code1[128][3]={[0 ... 127][0 ... 2]=0};
	for (jj=0;jj<8;jj++){code[8+jj][0]=(code[jj][0]^1);code[8+jj][1]=(code[jj][1]^1);code[8+jj][2]=(code[jj][2]^0);}
	for (jj=0;jj<16;jj++){code[16+jj][0]=(code[jj][0]^1);code[16+jj][1]=(code[jj][1]^1);code[16+jj][2]=(code[jj][2]^1);}
	for (jj=0;jj<32;jj++){code[32+jj][0]=(code[jj][0]^0);code[32+jj][1]=(code[jj][1]^1);code[32+jj][2]=(code[jj][2]^1);}
	for (jj=0;jj<64;jj++){code[64+jj][0]=(code[jj][0]^1);code[64+jj][1]=(code[jj][1]^1);code[64+jj][2]=(code[jj][2]^1);}
	for (jj=0;jj<128;jj++){code1[jj][0]=1/sqrt(2)*(2*code[jj][0]-1);code1[jj][1]=1/sqrt(2)*(2*code[jj][1]-1);code1[jj][2]=1/sqrt(2)*(2*code[jj][2]-1);}
	//FILE *f;  //debugger
	//f = fopen("viterbi_debug.txt","w");
	//for (kk=0; kk<128; kk++){fprintf(f,"%f %f %f\n",code1[kk][0],code1[kk][1],code1[kk][2]);}
			// initialize the state and survpath
	float state_trellis[64][length+1];
	int surv_path[64][length];
	for (jj=0;jj<64;jj++){surv_path[jj][0]=0;state_trellis[jj][0]=0;}
	float branchmetric[4]={0};
			// start decoding
	for (kk=0;kk<4;kk++){
	for (ii=0;ii<length;ii++){
		float eucl_dist[128];
		for (jj=0;jj<128;jj++){eucl_dist[jj]=-(code1[jj][0]*p_zerobits[ii])-(code1[jj][1]*p_onebits[ii])-(code1[jj][2]*p_twobits[ii]);}
		//fprintf(f,"\n eucle distance for instance %d \n",ii+1);//debugger
		//for (kk=0; kk<128; kk++){fprintf(f,"%f",eucl_dist[kk]);}
		//fprintf(f,"\n branch metrics \n");
		for (jj=0;jj<32;jj++){
			branchmetric[0]=eucl_dist[2*jj]-eucl_dist[64+2*jj]+state_trellis[2*jj][ii];
			branchmetric[1]=eucl_dist[2*jj+1]-eucl_dist[64+2*jj+1]+state_trellis[2*jj+1][ii];
			branchmetric[2]=eucl_dist[64+2*jj]-eucl_dist[2*jj]+state_trellis[2*jj][ii];
			branchmetric[3]=eucl_dist[64+2*jj+1]-eucl_dist[2*jj+1]+state_trellis[2*jj+1][ii];
			//fprintf(f,"\n %f %f %f %f \n",branchmetric[0],branchmetric[1],branchmetric[2],branchmetric[3]);//debugger
			int idx1=0;idx1=(branchmetric[0]<branchmetric[1])?0:1;int idx2=0;idx2=(branchmetric[2]<branchmetric[3])?0:1;
			state_trellis[jj][ii+1]=branchmetric[idx1];state_trellis[32+jj][ii+1]=branchmetric[idx2+2];
			surv_path[jj][ii]=idx1+2*(jj);surv_path[32+jj][ii]=idx2+2*(jj);
			//fprintf(f,"\n %f %f ",state_trellis[jj][ii+1],state_trellis[32+jj][ii+1]);//debugger
		}
	}
	float minima=0;
	minima=state_trellis[0][length];
	for(jj=1;jj<64;jj++){if ( state_trellis[jj][length] < minima ) {minima = state_trellis[jj][length];location = jj;}}
			// traceback unit
			//printf("%d %f",location,minima);
	int prevstate=0;currentstate=location;
	for (ii=length-1;ii>=0;ii--){
		prevstate=surv_path[currentstate][ii];
		outputBitStream[ii]=LUT[prevstate][currentstate];
		currentstate=prevstate;
	}
	
	if (currentstate==location){break;}
	for (ii=0;ii<64;ii++){state_trellis[ii][0]=state_trellis[ii][length];}
	}
	//FILE *f1;
	//f1 = fopen("viterbi_output.txt","w");
	//fprintf(f1,"ending state is %d \n",location);
	//fprintf(f1,"starting state is %d \n",currentstate);
	for (ii=0; ii<length; ii++){
		//fprintf(f1,"%c \n",*(outputBitStream+ii));
		printf("%d",*(outputBitStream+ii));
	}
	//fclose(f1);
	//fclose(f);
}
