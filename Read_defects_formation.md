#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#define N_DIR    4      // sampling number over the displacement direction)
#define N_TIMING 4      // sampling number over the timing
#define D_TIMING 50      // difference in neighboring timings

//#define PKAID1 171
#define PKA_ENE_MIN 6000   // integeri
#define PKA_ENE_INC 6000   // integer
#define PKA_ENE_MAX 6000 //121    // integer

#define LAM_EXE "/home/bin/lammps"
#define PROJ "data.PKA_D"

#define DATA_FILE "/home2/mosab/2017.01.04_V_MO/2017_12_12_W_STRAIN_damage_evolution/hydro/15x15x15/1.6%/ALL_DEFECTS_DATA/dump.voronoi_"
#define NA 6750

#define RUN_NUM 10
#define RUN_INT 100
int dir_num=0;//N_DIR-1;
int time_num=0;//N_TIMING-1;
int main(void)
{

	int i,p,t,k,t1,t2,n,s,discard_nrun,j;
	char cline[500],data_file[500],file_name[500];
	int na, atom_N,m;
        FILE *fr;
        int   a=PKA_ENE_MAX/PKA_ENE_INC;
        double v[16][20]; 
	double  D=N_TIMING*N_DIR;
	double pka_ene,v1;
	int in=0;

	for(i=0;i<16;i++)
        {
             
		for(k=0;k<20;k++) { v[i][k]=0;}  
               
        }

	for(n=PKA_ENE_MIN;n<=PKA_ENE_MAX;n=n+PKA_ENE_INC)
	{	
     for(s=0;s<N_DIR;s++)
	{
	        for(t=0;t<N_TIMING;t++)
        		{
               			discard_nrun=t*D_TIMING;


		
        				{

				pka_ene=n;
				sprintf(file_name,"PKA_D%.4d_E%.3d_T%.4d",s,n,discard_nrun);

				sprintf(data_file,"%sPKA_D%.4d_E%.3d_T%.4d",DATA_FILE,s,n,discard_nrun);
                               // printf("%s\n",file_name);
    					if( (fr=fopen(data_file, "r"))==NULL )
     					  {
           				  printf("error in open data file: %s%m\n",data_file); exit(1);
              				  printf("Error: %d (%s) %s\n",errno, strerror(errno),data_file);exit(1);
     				          }
                                        for(t2=0;t2<101;t2=t2+1){
					
					
					v1=0;
					
			           	fgets(cline,500,fr);
               				fscanf(fr,"%d",&t1);
              
         			        for(i=0;i<9;i++)    {fgets(cline,500,fr);}
               				for(i=0;i<NA+1;i++) { fscanf(fr,"%d", &atom_N);{if(atom_N==0)  v1=v1+1;}}
					if (t1==0 || t1==1000 || t1==2000 || t1==3000 || t1==4000 || t1==5000 || t1==6000 || t1==7000 || t1==8000 || t1==9000 || t1==10000 || t1==11000
						|| t1==12000 || t1==13000 || t1==14000 || t1==15000 || t1==16000 || t1==17000 || t1==18000 || t1==19000 || t1==20000 ){	
		
						v[(in+t)][(t1/1000)]+=v1;}
         				       
					}
              				fclose(fr);

              			}
			} in=in+4;
	}
}
	
	for(i=0;i<20;i++) { printf("\t%d",(i*1000)); }
	printf("\n");
//	for(p=0;p<N_DIR;p++){
//	for(i=0;i<a;i++)
 	for(i=0;i<16;i++)

	{
		printf("%d", (PKA_ENE_INC*i)+PKA_ENE_MIN);
		for(k=0;k<21;k++) { v[i][k]=v[i][k]; printf("\t%lf", v[i][k]); }
		printf("\n");
	}
			
return 0 ;

}



