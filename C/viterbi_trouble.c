/*  Viterbi Decoder
 * alll viterbi
 * Written by Dama Sreekanth
 */
#include <omp.h>
#define chunk 1
#define N     26
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "source/headers/viterbi_wava.h"
#include "source/headers/viterbi_ml.h"
#include "source/headers/viterbi_new_met.h"
#include "source/headers/tail_enc.h"
//#include "source/headers/viterbi_zero.h"

/* Tail biting convolution */
int main()
{
	float pi = 3.14159265358979;int average = 100000;
	int length=80;
	int MIB[80]={0};int ii;float gen_rand_uniform();
	//float EbnodB[10]={0.891250938133746,0.794328234724282,0.707945784384138,0.630957344480193,0.562341325190349,0.501187233627272,0.446683592150963,0.398107170553497,0.354813389233575,0.316227766016838};
	int jj;printf("\n");float EbnodB[26]={0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
	int nerror[26]={0},nerror1[26]={0},nerror2[26]={0},nerror3[26]={0},diff[80]={0},diff1[80]={0},diff2[80]={0},diff3[80]={0};int kk=0;
	float FER[26]={0},FER1[26]={0},FER2[26]={0},FER3[26]={0};
	void viterbi_wava();void viterbi_ml();void tail_enc();void viterbi_new_met();//void viterbi_zero();
	float te =0;float te1=0;float te2=0;//float te3=0;
							float *C0=0,*C1=0,*C2=0;//,*cc1=0,*cc2=0,*cc3=0;
							C0 = (float *)malloc((length)*sizeof(float));
							C1 = (float *)malloc((length)*sizeof(float));
							C2 = (float *)malloc((length)*sizeof(float));
							//cc1 = (float *)malloc((length+7)*sizeof(float));
							//cc2 = (float *)malloc((length+7)*sizeof(float));
							//cc3 = (float *)malloc((length+7)*sizeof(float));
							int *output,*output1,*output2;//,*output3;
							output = (int *)malloc(length*sizeof(int));
							output1 = (int *)malloc(length*sizeof(int));
							output2 = (int *)malloc(length*sizeof(int));
							//output3 = (int *)malloc((length+7)*sizeof(int));
FILE *f;f = fopen("BER_FER.txt","w");
	// Run the program for all SNRS
	#pragma omp parallel shared(MIB,C0,C1,C2,length,FER,FER1,FER2,nerror,nerror1,nerror2) private(jj,kk,ii)
	{
			#pragma omp for schedule(dynamic,chunk) nowait
			for (jj=0;jj<N;jj++){
				//printf("\n SNR %d  %f\n",jj,EbnodB[jj]);FER[jj]=0;FER1[jj]=0;FER2[jj]=0;

					for (kk=0;kk<average;kk++){ //average over n
							/*printf("%d\n",kk);*/te=0;te1=0;te2=0;//te3=0;

							for (ii=0;ii<length;ii++){float u1 = gen_rand_uniform();int b = (u1>0.5)?1:0;MIB[ii]=b;/*printf("%d",MIB[ii]);*/}





							tail_enc(MIB,length,C0,C1,C2);//conv_enc(MIB,length,cc1,cc2,cc3);



							/*FILE *f;
							f = fopen("encode_out.txt","w");
							for (ii=0; ii<length; ii++){
							fprintf(f,"%f %f %f \n",*(C0+ii),*(C1+ii),*(C2+ii));
							}
							fclose(f);*/


							//* channel

							for (ii=0;ii<length;ii++){
								float temp1 = gen_rand_uniform();float temp2 = gen_rand_uniform();float r = sqrt(-2*log(temp1));
								float n_re = r*cos(2*pi*temp2);float n_im = r*sin(2*pi*temp2);
								*(C0+2*ii) = *(C0+2*ii)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_re;*(C0+2*ii+1) = *(C0+2*ii+1)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_im;
								//*(cc1+2*ii) = *(cc1+2*ii)+(1/sqrt(2))*EbnodB[jj]*n_re;*(cc1+2*ii+1) = *(cc1+2*ii+1)+(1/sqrt(2))*EbnodB[jj]*n_im;}
							}


							for (ii=0;ii<length;ii++){
								float temp1 = gen_rand_uniform();float temp2 = gen_rand_uniform();float r = sqrt(-2*log(temp1));
								float n_re = r*cos(2*pi*temp2);float n_im = r*sin(2*pi*temp2);
								*(C1+2*ii) = *(C1+2*ii)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_re;*(C1+2*ii+1) = *(C1+2*ii+1)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_im;
								//*(cc2+2*ii) = *(cc2+2*ii)+(1/sqrt(2))*EbnodB[jj]*n_re;*(cc2+2*ii+1) = *(cc2+2*ii+1)+(1/sqrt(2))*EbnodB[jj]*n_im;}

						}

							for (ii=0;ii<length;ii++){
								float temp1 = gen_rand_uniform();float temp2 = gen_rand_uniform();float r = sqrt(-2*log(temp1));
								float n_re = r*cos(2*pi*temp2);float n_im = r*sin(2*pi*temp2);
								*(C2+2*ii) = *(C2+2*ii)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_re;*(C2+2*ii+1) = *(C2+2*ii+1)+(1/sqrt(2))*pow(10,-EbnodB[jj]/20)*n_im;
								//*(cc3+2*ii) = *(cc3+2*ii)+(1/sqrt(2))*EbnodB[jj]*n_re;*(cc3+2*ii+1) = *(cc3+2*ii+1)+(1/sqrt(2))*EbnodB[jj]*n_im;}

						}


							viterbi_wava(C0,C1,C2,length,output);

							viterbi_ml(C0,C1,C2,length,output1);

							viterbi_new_met(C0,C1,C2,length,output2);
							//viterbi_zero(cc1,cc2,cc3,length+6,output3);



							// bit error

							for (ii=0;ii<length;ii++){diff[ii]=(*(output+ii)+MIB[ii])%2;nerror[jj]=nerror[jj]+diff[ii];te=te+diff[ii];
								diff1[ii]=(*(output1+ii)+MIB[ii])%2;nerror1[jj]=nerror1[jj]+diff1[ii];te1=te1+diff1[ii];
								diff2[ii]=(*(output2+ii)+MIB[ii])%2;nerror2[jj]=nerror2[jj]+diff2[ii];te2=te2+diff2[ii];
								//diff3[ii]=(*(output3+ii)+MIB[ii])%2;nerror3[jj]=nerror3[jj]+diff3[ii];te3=te3+diff3[ii];
								/*printf("%d%d%d\n",diff[ii],diff1[ii],diff2[ii]);*/}//printf("\n no of errors = %d",te3);
					 te = ((te/80)>0)?1:0;te1 = ((te1/80)>0)?1:0;te2 = ((te2/80)>0)?1:0;//te3 = ((te3/80)>0)?1:0;
					FER[jj]=FER[jj]+te;FER1[jj]=FER1[jj]+te1;FER2[jj]=FER2[jj]+te2;//FER3[jj]=FER3[jj]+te3;
				}
			//for (ii=0;ii<average;ii++){error[jj]=error[jj]+nerror[jj][ii];}printf("\nNoise %f\nno of errors = %d\n",(1/sqrt(2))*EbnodB[jj],error[jj]/average);
			//printf("good bye length %.2f average ",(float)length*average);
			//free(C0,C1,C2)
			//printf("%f  %f  %f",FER[jj],FER1[jj],FER2[jj]);

			/*printf("\nNoise %f\nBER of WAVA = %f\nBER of ML = %f \nBER of new = %f\n",
			(1/sqrt(2))*EbnodB[jj],((float)nerror[jj])/((float)(length*average)),((float)nerror1[jj])/((float)(length*average)),
			((float)nerror2[jj])/((float)(length*average)));//,((float)nerror3[jj])/((float)(length*average)));
			printf("\nFER of WAVA = %f\nFER of ML = %f \nFER of new = %f \n",
			FER[jj]/average,FER1[jj]/average,FER2[jj]/average);//,FER3[jj]/average);*/
			fprintf(f,"\nBER FER WAVA ML New Method\n%f%f %f%f %f%f %f %f",((float)nerror[jj])/((float)(length*average)),FER[jj]/average,
			((float)nerror1[jj])/((float)(length*average)),FER1[jj]/average,
			((float)nerror2[jj])/((float)(length*average)),FER2[jj]/average);
			//((float)nerror3[jj])/((float)(length*average)),FER3[jj]/average);
			//fprintf(f,"\nFER\n%f %f %f",FER[jj]/average,FER1[jj]/average,FER2[jj]/average);
			//printf("good bye");

			}
fclose(f);
			// end of the code
		}
			//printf("\ngood bye\n");
	return 0;
}



float gen_rand_uniform (void) {

	return (double)rand()/((double)RAND_MAX + (double)1);
}
