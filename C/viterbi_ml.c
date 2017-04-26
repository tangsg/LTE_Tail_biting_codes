/*
ML Viterbi Algorithm for tail biting codes
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



void viterbi_ml(float p_zerobits[],float p_onebits[],float p_twobits[],unsigned int length,int *outputBitStream)
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
	float state_trellis[64][length+1];float state_ml[64]={0};
	int surv_path[64][length][64];
	for (jj=0;jj<64;jj++){surv_path[jj][0][jj]=0;state_trellis[jj][0]=0;}
	float branchmetric[4]={0};
			// start decoding
	for (kk=0;kk<64;kk++){//printf("\n");
	for (ii=0;ii<length;ii++){
		float eucl_dist[128];
		for (jj=0;jj<128;jj++){eucl_dist[jj]=-(code1[jj][0]*p_zerobits[ii])-(code1[jj][1]*p_onebits[ii])-(code1[jj][2]*p_twobits[ii]);}
		//fprintf(f,"\n eucle distance for instance %d \n",ii+1);//debugger
		//for (kk=0; kk<128; kk++){fprintf(f,"%f",eucl_dist[kk]);}
		//fprintf(f,"\n branch metrics \n");
		
		if (ii<6){
			
			int x = ((int)ceil((float)(kk+1)/(1<<ii))+1)%2;//printf("%d",x);
			for (jj=0;jj<32;jj++){
				//printf(" %d ",2*jj+x);
			branchmetric[0]=eucl_dist[2*jj+x]-eucl_dist[64+2*jj+x]+state_trellis[2*jj+x][ii];
			branchmetric[2]=eucl_dist[64+2*jj+x]-eucl_dist[2*jj+x]+state_trellis[2*jj+x][ii];
			//fprintf(f,"\n %f %f %f %f \n",branchmetric[0],branchmetric[1],branchmetric[2],branchmetric[3]);//debugger
			state_trellis[jj][ii+1]=branchmetric[0];state_trellis[32+jj][ii+1]=branchmetric[2];
			surv_path[jj][ii][kk]=2*jj+x;surv_path[32+jj][ii][kk]=2*jj+x;
			//fprintf(f,"\n %f %f ",state_trellis[jj][ii+1],state_trellis[32+jj][ii+1]);//debugger
		}
		}
		else
		{
			for (jj=0;jj<32;jj++){
			branchmetric[0]=eucl_dist[2*jj]-eucl_dist[64+2*jj]+state_trellis[2*jj][ii];
			branchmetric[1]=eucl_dist[2*jj+1]-eucl_dist[64+2*jj+1]+state_trellis[2*jj+1][ii];
			branchmetric[2]=eucl_dist[64+2*jj]-eucl_dist[2*jj]+state_trellis[2*jj][ii];
			branchmetric[3]=eucl_dist[64+2*jj+1]-eucl_dist[2*jj+1]+state_trellis[2*jj+1][ii];
			//fprintf(f,"\n %f %f %f %f \n",branchmetric[0],branchmetric[1],branchmetric[2],branchmetric[3]);//debugger
			int idx1=0;idx1=(branchmetric[0]<branchmetric[1])?0:1;int idx2=0;idx2=(branchmetric[2]<branchmetric[3])?0:1;
			state_trellis[jj][ii+1]=branchmetric[idx1];state_trellis[32+jj][ii+1]=branchmetric[idx2+2];
			surv_path[jj][ii][kk]=idx1+2*(jj);surv_path[32+jj][ii][kk]=idx2+2*(jj);}
			//fprintf(f,"\n %f %f ",state_trellis[jj][ii+1],state_trellis[32+jj][ii+1]);//debugger
		}
		
	}
	state_ml[kk]=state_trellis[kk][length];//printf("%d %f \n",kk,state_ml[kk]);
	}
	float minima=0;
	minima=state_ml[0];location=0;
	for(jj=1;jj<64;jj++){if ( state_ml[jj] < minima ) {minima = state_ml[jj];location = jj;}}
			// traceback unit
			//printf("%d %f \n",location,minima);
	int prevstate=0;currentstate=location;
	for (ii=length-1;ii>=0;ii--){
		prevstate=surv_path[currentstate][ii][location];
		outputBitStream[ii]=LUT[prevstate][currentstate];
		currentstate=prevstate;
	}
	
	//FILE *f1;
	//f1 = fopen("viterbi_output.txt","w");
	//fprintf(f1,"ending state is %d \n",location);
	//fprintf(f1,"starting state is %d \n",currentstate);
	//for (ii=0; ii<length; ii++){
		////fprintf(f1,"%c \n",*(outputBitStream+ii));
		//printf("%d",*(outputBitStream+ii));
	//}
	//fclose(f1);
	//fclose(f);
}
