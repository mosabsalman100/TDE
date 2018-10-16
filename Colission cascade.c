//for thsi part apply the recenter_fix  command for lammps  in the "input.voronoi_ref_fix_recenter" .
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define N_DIR    4	// sampling number over the displacement direction)
#define N_TIMING 4 	// sampling number over the timing
#define D_TIMING 50	// difference in neighboring timings

#define PKA_DIRECTION_FOLDER "/home2/mosab/2017.01.04_V_MO/2017_12_12_W_STRAIN_damage_evolution/hydro/15x15x15/resources"

#define LAMMPS_EXE "mpirun -np 16 /home/bin/lammps_mpi"			// lammps execution command
#define INP_LAMMPS_REF       "/home2/mosab/2017.01.04_V_MO/2017_12_12_W_STRAIN_damage_evolution/hydro/15x15x15/resources/input.voronoi_ref_recenter3"// reference input file (to be read and then be modified)(here we introduced voroni method to get the number of Interstitial and Vacancy and know the TDE 
#define NLINE_INP_LAMMPS_REF 55						// the number of lines in the reference input file

//#define PROJ1 "W-10x10x15_NVE-30K_165eV_d3711_dt_Xmax01_testtest"
#define RES_FILE1 "/home2/mosab/2017.01.04_V_MO/2017_12_12_W_STRAIN_damage_evolution/hydro/15x15x15/resources/restart.strain_W_16"	// restart file name
#define POT_FILE1 "/home/mosab/resources/W_BN.eam.fs"		// potential file name
#define ATOM_SYMBOL1 "W"	// atomic symbol
#define ATOM_MASS1 183.84	// atomic mass
#define TMAX1 0.001		// dt/reset: maximum time step in ps
#define XMAX1 0.01		// dt/reset: maximum displacement allowed for a step in A
#define DT_FREQ1 10  		// dt/reset: the frequency of timestep redefinition
#define NRUN1 20000		// the number of run after the recoil initiation
#define RES_FREQ1 200		// the frequency of restart file generation
#define THERMO_FREQ1 5		// the frequency of the thermo information pring
#define DISCARD_NRUN1 0		// the number of run to change the timing of recoil event
#define DISCARD_TSTEP1 TMAX1	// the timestep for the run to change the timing


#define CON_ELEMENTARY_CHARGE 1.60217662E-19
#define CON_AVOGADRO 6.022140857E23
#define RANDOM_SEED 10

/////////////////
#define PKAID1 6750// 512
#define PKA_ENE_MIN 6000	// integer
#define PKA_ENE_INC 8000	// integer
#define PKA_ENE_MAX 7000 //121	// integer
/////////////////

double func_calcuate_pka_velocity(double pka_mass, double pka_ene, double pka_vec[], double pka_vel[]);
void func_write_execution_command(char proj1[],int discard_nrun);
void func_create_input_files(char fold1[], char proj1[], double pka_ene, double pka_vel[], int pkaid, int discard_nrun);



int main(void)
{
	int i,j,n,s,t;
	int ns;
	int discard_nrun;
	char cline1[500];
	char fold1[300];
	char command1[300];
	double pka_vx,pka_vy,pka_vz,pka_speed;
	double pka_mass,pka_ene;
	double pka_vec[N_DIR][3],pka_vel[3];
	char fname_pkadir[300];
	char proj1[300];
	char fname1[300];
	FILE *fr,*fw;
	
	srand(RANDOM_SEED);

	sprintf(fname_pkadir,"%s/PKA-%.5d.dat",PKA_DIRECTION_FOLDER,N_DIR);
	if( (fr=fopen(fname_pkadir,"r"))==NULL )	{  printf("error in pka-direction file: %s\n",fname_pkadir);  exit(1);  }
        for(s=0;s<N_DIR;s++)
        {
                fscanf(fr,"%d",&ns);
                if(s!=ns)       {  printf("error: mismatch in ns number\n");  exit(1);  }
                for(i=0;i<3;i++)        fscanf(fr,"%lf",&pka_vec[s][i]);
	}
	fclose(fr);
	
	printf("mkdir ALL_DEFECTS_DATA\n");
   	for(t=0;t<N_TIMING;t++)
       	{
         	discard_nrun=t*D_TIMING;
              	printf("#Timing:\t%d\n",discard_nrun);
		sprintf(fold1,"TIMING_%.4d",discard_nrun);
		sprintf(command1,"mkdir %s",fold1);
		system(command1);

		for(s=0;s<(N_DIR);s++)
		{
			printf("#DIR\t%.3d\n",s);
			printf("#PKA-displacement-direction\t%lf\t%lf\t%lf\n",pka_vec[s][0],pka_vec[s][1],pka_vec[s][2]);

			for(n=PKA_ENE_MIN;n<=PKA_ENE_MAX;n=n+PKA_ENE_INC)
			{
				pka_ene=n;
				printf("#Energy: %lf\n",pka_ene);
				pka_speed=func_calcuate_pka_velocity(ATOM_MASS1, pka_ene, pka_vec[s], pka_vel);

				sprintf(proj1,"PKA_D%.4d_E%.3d_T%.4d",s,n,discard_nrun);
				func_create_input_files(fold1, proj1, pka_ene, pka_vel, PKAID1, discard_nrun);

				func_write_execution_command(proj1,discard_nrun);
			}
		}
	}


	return discard_nrun;
}


void func_create_input_files(char fold1[], char proj1[], double pka_ene, double pka_vel[], int pkaid, int discard_nrun)
{
	int i,j,k;
	char cline1[500],fname1[300];
	FILE *fw,*fr;

    	sprintf(fname1,"%s/input.%s",fold1,proj1);
      	if( (fw=fopen(fname1,"w"))==NULL )      {  printf("error in open %s file\n",fname1);  exit(1);  }
      	fprintf(fw,"# pka-energy[eV]:\t%lf\n",pka_ene);
     	fprintf(fw,"# %lf\t%lf\t%lf\n",pka_vel[0],pka_vel[1],pka_vel[2]);

   	if( (fr=fopen(INP_LAMMPS_REF,"r"))==NULL )      {  printf("error in open INP_LAMMPS_REF file\n");  exit(1);  }

      	for(i=0;i<NLINE_INP_LAMMPS_REF;i++)
      	{
          	fgets(cline1,500,fr);
              	if(i==1)        fprintf(fw,"variable proj1  string       %s\n",proj1);
              	else if(i==2)   fprintf(fw,"variable pkaid1 equal        %d\n",pkaid);
             	else if(i==3)   fprintf(fw,"variable pkavx1 equal        %lf\n",pka_vel[0]);
              	else if(i==4)   fprintf(fw,"variable pkavy1 equal        %lf\n",pka_vel[1]);
            	else if(i==5)   fprintf(fw,"variable pkavz1 equal        %lf\n",pka_vel[2]);
              	else if(i==6)   fprintf(fw,"variable res_file1 string    %s\n",RES_FILE1);
              	else if(i==7)   fprintf(fw,"variable pot_file1 string    %s\n",POT_FILE1);
              	else if(i==8)   fprintf(fw,"variable atom_symbol1 string %s\n",ATOM_SYMBOL1);
              	else if(i==9)   fprintf(fw,"variable atom_mass1   equal  %lf\n",ATOM_MASS1);
             	else if(i==10)  fprintf(fw,"variable tmax1 equal         %lf\n",TMAX1);
             	else if(i==11)  fprintf(fw,"variable xmax1 equal         %lf\n",XMAX1);
           	else if(i==12)  fprintf(fw,"variable nrun1 equal         %d\n",NRUN1);
             	else if(i==13)  fprintf(fw,"variable res_freq1 equal     %d\n",RES_FREQ1);
             	else if(i==14)  fprintf(fw,"variable thermo_freq1 equal  %d\n",THERMO_FREQ1);
             	else if(i==15)  fprintf(fw,"variable dt_freq1     equal  %d\n",DT_FREQ1);
            	else if(i==16)  fprintf(fw,"variable discard_nrun1  equal %d\n",discard_nrun);
             	else if(i==17)  fprintf(fw,"variable discard_tstep1 equal %lf\n",DISCARD_TSTEP1);
              	else if(i==18)   fprintf(fw,"variable pka_ene equal       %lf\n",pka_ene);

             	else    fprintf(fw,"%s",cline1);
      	} 
     	fclose(fr);
     	fclose(fw);
}


void func_write_execution_command(char proj1[],int discard_nrun)
{
	int n;
	
	printf("echo \"%s\"\n",proj1);
     	printf("date\n");
//	printf("mkdir ALL_DEFECTS_DATA\n");

	printf("cd TIMING_%.4d /\n",discard_nrun);
      	printf("%s <input.%s >output.voronoi_%s\n",LAMMPS_EXE,proj1,proj1);
      	printf("mkdir result_%s\n",proj1);
	
      	printf("mv restart.*000 result_%s\n",proj1,proj1);
//      printf("rm restart.*\n");
//	printf("cd  result_%s/\n",proj1);
//	for(n=1000;n<(NRUN1+1000);n=n+1000)
//	{
//      printf("%s -r  restart.%s.%d  data.%s.%d\n",LAMMPS_EXE,proj1,n,proj1,n);
//	printf("mv  data.%s.%d ../../ALL_DEFECTS_DATA/\n",proj1,n);
//		}

//	printf("cd ../\n");	
//      printf("mv output.voronoi_%s result_%s\n",proj1,proj1);
	printf("mv dump.voronoi_%s //home2/mosab/2017.01.04_V_MO/2017_12_12_W_STRAIN_damage_evolution/hydro/15x15x15/1.6%/ALL_DEFECTS_DATA/\n",proj1);
     // printf("mv dump.PKA_%s   result_%s\n",proj1,proj1);
        printf("rm log.*\n");
        printf("rm output.voronoi_*\n");
        printf("rm dump.PKA_*\n");



       //	printf("mv log.%s    result_%s\n",proj1,proj1);
	//printf("mv log.lammps log.%s_ini\n",proj1);
	//printf("mv log.%s_ini result_%s\n",proj1,proj1);
	printf("mv input.%s result_%s\n",proj1,proj1);
        printf("cd ../\n");
}

double func_calcuate_pka_velocity(double pka_mass, double pka_ene, double pka_vec[], double pka_vel[])
{
	int i;
	double pka_speed;

        pka_speed=sqrt((2.0*(pka_ene*CON_ELEMENTARY_CHARGE))/(0.001*pka_mass/CON_AVOGADRO));    // in [m/s]
        pka_speed=pka_speed*(1.0E10*1.0E-12);   // in [A/ps]
        for(i=0;i<3;i++)        pka_vel[i]=pka_speed*pka_vec[i]/sqrt(pow(pka_vec[0],2.0)+pow(pka_vec[1],2.0)+pow(pka_vec[2],2.0));

	return pka_speed;
}

