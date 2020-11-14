#include "stdio.h"
#include "string.h"
#include "math.h"
#include <stdlib.h>
#include "R.h" 

void bounded_monotone_fndr(double *vector,int size);
void bounded_monotone_fdr(double *vector,int size);
int least_concave_majorant(double *x,double *F, int first_point, int final_point, double *F_g, double *f_g, int size, double *f);
int greatest_convex_minorant(double *x,double *F, int first_point, int final_point, double *F_g, double *f_g, int size, double *f);
double grenander_error(double *F, double *F_g, double *x, int size);
double choose_fit(double *x, int grid_size, double *F, double *f,	double *F_g, double *f_g, int flag, int *flag_choice, int *index_p_prime);
int small_grid_calculation(double *f, double *x, int grid_size, double *f_small, double *x_small, int N);
int explained_vs_unexplained_variance(double *x, double *x_small, int p_prime_index, int n_small, double *f, double *f_g);
double goodness_calc(double *f, double *f_g, double *x, int n);
void correcting_slopes(double *f_g, double *x, int size);
void calculate_range(double x_min, double x_max, double *x, int *range, int n);

int global_flag=0;

extern void get_grids_full(char **file,
					int *columns,
					double *full_tests,				
					double *rows,
					int *stat_column,
					int *grid_size,
					double *grid_unif,
					double *grid_log,
					double *F_log,
					double *f_log,
					double *F_unif,
					double *f_unif,
					double *min_stat,
					double *max_stat,
					int *zero_count,
					int *print_messages,
					int *rownames_flag)
{		
	/* This function loads the p-value distribution in a uniform and a logarithmic grid*/
	/* This function is used when the data is given as a text file and NOT as a data frame*/
		
	zero_count[0]=0;
		
	//1. Definition of the bins in the Uniform Grid
	int i;
	double delta=(double)(max_stat[0]-min_stat[0])/(grid_size[0]-1);
	for(i=0;i<grid_size[0];i++)
		grid_unif[i]=min_stat[0]+i*delta;
	
	//2. Definition of the bins in the Logarithmic grid
	double log_max_ratio=log10(max_stat[0]/min_stat[0]);
	double delta_log=log_max_ratio/(grid_size[0]-1);
	for(i=0;i<grid_size[0];i++)
		grid_log[i]=min_stat[0]*pow(10,i*delta_log); 
		
	//3. Initialize distributions for the Probability Density Function (PDF, lower-case f) and
	//... the Cumulative Distribution Function (CDF, capital F) at both grids
	for(i=0;i<grid_size[0];i++)
	{
		f_unif[i]=0;
		f_log[i]=0;
		F_unif[i]=0;
		F_log[i]=0;
	}
	
	FILE *fp;
    int j;
    int line_count=1;
    double progress=0;
	double stat_value, log_stat_value;
	char word[1000], header[10000];
	int bin;
	double stat_count_fu=0,stat_count_flog=0;
	double decimal_part;
			
	fp=fopen(*file,"r");
	for(i=0;i<(columns[0]-rownames_flag[0]);i++)
		fscanf(fp,"%s ",header);
	
	
	//4. Load the different distributions
	while(feof(fp)==0)
	{		
		//4.1 Reads not important columns before the p-value
		for(j=0;j<(stat_column[0]-1);j++)
			fscanf(fp,"%s ",word);		
				 
		//4.2 Reads the p-value       
        fscanf(fp,"%lf ",&stat_value);
						
		//4.3 Update the counter of exact 0
		// ... we need to do this outside the regular grid, because we can't use 0 in log space
		if(stat_value==0)
			zero_count[0]++;
		
		//4.4 Transformation in the log space	
		log_stat_value=log10(stat_value/min_stat[0]);
			
		//4.5 Cumulative Distribution Function (CDF) update in Uniform and Logarithmic Grids
		// ... As a first step we simply add 1 in the bin where the data belongs
		// ... After that we will need to add the value of each bin in the higher bins and normalize
		if((min_stat[0]<stat_value)&&(stat_value<max_stat[0]))
		{
			bin=(int)((stat_value-min_stat[0])/delta);
			decimal_part=((stat_value-min_stat[0])/delta)-bin;
			if(decimal_part!=0)
				bin++;
			F_unif[bin]++;
				
			if((bin>=grid_size[0])||(bin<1))
			{
				printf("Error bin uniform\n");
				printf("Stat %.16lf --> min_stat %.16lf delta %.16lf\n",stat_value,min_stat[0],delta);
				printf("Bin %d (Max %d)\n",bin,grid_size[0]-1);
				getchar();
			}
							
			bin=(int)(log_stat_value/delta_log);
			decimal_part=log_stat_value/delta_log-bin;
			if(decimal_part!=0)
				bin++;
								
			F_log[bin]++;
			
			if((bin>=grid_size[0])||(bin<1)) 
			{
				printf("Error bin log\n");	
				printf("Stat %.16lf --> log_stat %.16lf delta_log %.16lf\n",stat_value,log_stat_value,delta_log);
				printf("Bin %d (Max %d)\n",bin,grid_size[0]-1);
				getchar();
			}			
		}
		else
		{
			if(stat_value<=min_stat[0])
			{
				F_unif[0]++;
				F_log[0]++;
			}
			if(stat_value==max_stat[0])
			{
				F_unif[grid_size[0]-1]++;
				F_log[grid_size[0]-1]++;
			}
		}
			
		//4.6 Probability Density Function (PDF) update in Uniform Grid
		//... Adds 1 if the value is at the vicinity (+- 0.5*delta) of the bin point
		//... After that we will need to normalize
		if(((min_stat[0]-delta*0.5)<stat_value)&&(stat_value<(max_stat[0]-0.5*delta)))
		{
			bin=(int)((stat_value-min_stat[0]+0.5*delta)/delta);
			f_unif[bin]++;
			
			if(bin>=grid_size[0]) printf("Error bin uniform PDF\n");
		}
		else
		{
			//The last bin is set separatedly, otherwise the maximum value will end up outside the array
			if(((max_stat[0]-delta*0.5)<stat_value)&&(stat_value<(max_stat[0]+0.5*delta)))
				f_unif[grid_size[0]-1]++;
			//Counter for data outside the min-max range (+- 0.5*delta).
			//Could be necessary for exact zeroes 
			else
				stat_count_fu++;
		}			
			
		//4.7 Probability Density Function (PDF) update in Logarithmic Grid
		//... Is analogous to 4.2 but in the logarithmic space 
		if(((-delta_log*0.5)<log_stat_value)&&(log_stat_value<(log_max_ratio-0.5*delta)))
		{
			bin=(int)((log10(stat_value/min_stat[0])+0.5*delta_log)/delta_log);
			f_log[bin]++;
			
			if(bin>=grid_size[0]) printf("Error bin log PDF\n");
		}
		else
		{
			if(((log_max_ratio-delta_log*0.5)<log_stat_value)&&(log_stat_value<(log_max_ratio+0.5*delta_log)))
				f_log[grid_size[0]-1]++;
			else
				stat_count_flog++;
		}				
						
		//4.8 Reads not important columns after the p-value
		for(j=(stat_column[0]);j<columns[0];j++)
            fscanf(fp,"%s ",word);
        			
        //4.9 Update Progress. 
        // ... Every 10 million data prints the percentage of the file that has been read
		if((line_count%10000000==0)&&(print_messages[0]==1))
        {
            progress=progress+(double)line_count*100/rows[0];
            line_count=0;
            printf("\nProgress: %.4lf %%",progress);
            fflush(stdout);
        }
		
		line_count++;		
		
	}
	if(print_messages[0]==1)
		printf("\nProgress: 100.00 %%\n");
	fclose(fp);
	
		
	//5. Normalizes and corrects distributions
	double F_unif_acum=0,F_log_acum=0;
	double area_unif=0,area_log=0;
	
	//5.1 CDF are corrected (we had a PDF-ish version)
	//... and normalized (dividing by the number of data) 		
	for(i=0;i<grid_size[0];i++)
	{
		F_unif_acum+=F_unif[i]/full_tests[0];
		F_unif[i]=F_unif_acum;
		F_log_acum+=F_log[i]/full_tests[0];
		F_log[i]=F_log_acum;
	}	
	
	//5.2 Normalization of PDF in the uniform grid
	//... Each bin is divided by the number of valid tests, and the bin width
	//... This forces the area of the PDF to be equal to 1
	if((grid_unif[0]-0.5*delta)<0)
	{
		f_unif[0]=f_unif[0]/((full_tests[0]-stat_count_fu)*(grid_unif[0]+0.5*delta));
		area_unif+=f_unif[0]*(grid_unif[0]+0.5*delta);
	}
	else
	{
		f_unif[0]=f_unif[0]/((full_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[0]*delta;	
	}
	for(i=1;i<(grid_size[0]-1);i++)
	{
		f_unif[i]=f_unif[i]/((full_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[i]*delta;		
	}
	if((grid_unif[grid_size[0]-1]+0.5*delta)>1)
	{
		f_unif[grid_size[0]-1]=f_unif[grid_size[0]-1]/((full_tests[0]-stat_count_fu)*(1-(grid_unif[grid_size[0]-1]-0.5*delta)));
		area_unif+=f_unif[grid_size[0]-1]*(1-(grid_unif[grid_size[0]-1]-0.5*delta));	
	}
	else
	{
		f_unif[grid_size[0]-1]=f_unif[grid_size[0]-1]/((full_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[grid_size[0]-1]*delta;	
	}
	//This is redundant
	area_unif=(1.0/area_unif);
	for(i=0;i<grid_size[0];i++)
		f_unif[i]=f_unif[i]*area_unif;	
	
	//5.3 Normalization of PDF in the logarithmic grid
	// ... In this case the width of the bins varies
	f_log[0]=(f_log[0]/(full_tests[0]-stat_count_flog))/(grid_log[1]-grid_log[0]);
	area_log+=f_log[0]*(grid_log[1]-grid_log[0]);
	for(i=1;i<grid_size[0]-1;i++)
	{
		f_log[i]=(f_log[i]/(full_tests[0]-stat_count_flog))*2/(grid_log[i+1]-grid_log[i-1]);
		area_log+=f_log[i]*0.5*(grid_log[i+1]-grid_log[i-1]);
	}
	f_log[grid_size[0]-1]=(f_log[grid_size[0]-1]/(full_tests[0]-stat_count_flog))/(grid_log[grid_size[0]-1]-grid_log[grid_size[0]-2]);
	area_log+=f_log[grid_size[0]-1]*(grid_log[grid_size[0]-1]-grid_log[grid_size[0]-2]);
}

extern void get_distribution_df(double *stats_array,
						int *data_size,
						int *processed_data_size,
						double *grid,
						int *grid_size, 
						double *f,
						double *F)
{
	/* This function loads the p-value distribution at each unique reported p-value (data grid)*/
	/* This function is used when data is given via a dataframe, NOT a file path*/
	
	int i,bin=0,count=0,neg_count=0;
		
	//1. Initialize PDF	
	for(i=0;i<grid_size[0];i++)
		f[i]=0;
	
	//2. Load CDF in Data Grid
	for(i=0;((i<processed_data_size[0])&&(bin<grid_size[0]));i++)
	{			
		if(stats_array[i]<=grid[bin])
			count++;		
		else
		{
			F[bin]=(double)count/data_size[0];			
			bin++;			
			i--;
		}
	}
	if(i>=processed_data_size[0])
	{
		for(;bin<grid_size[0];bin++)
			F[bin]=(double)count/data_size[0];
	}
	
	//3. Load PDF in Data Grid
	for(i=0;stats_array[i]<(grid[0]-0.5*(grid[1]-grid[0]));i++)	neg_count++;
	bin=0;	
	count=0;
	for(;(i<processed_data_size[0])&&(bin<(grid_size[0]-1));i++)
	{		
		if(stats_array[i]<(grid[bin]+0.5*(grid[bin+1]-grid[bin])))
			count++;		
		else
		{
			f[bin]=count;
			bin++;
			count=0;
			i--;
		}
	}
	count=0;
	for(;(stats_array[i]<(grid[grid_size[0]-1]+0.5*(grid[grid_size[0]-1]-grid[grid_size[0]-2])))&&(i<processed_data_size[0]);i++)
		count++;	
	f[grid_size[0]-1]=count;
	//neg_count is zero always?
	for(;i<processed_data_size[0];i++) neg_count++;
	
	//4. Normalize PDF	
	if((grid[0]-0.5*(grid[1]-grid[0]))<0)
		f[0]=f[0]/((data_size[0]-neg_count)*(grid[0]+0.5*(grid[1]-grid[0])));
	else
		f[0]=f[0]/((data_size[0]-neg_count)*(grid[1]-grid[0]));
	for(i=1;i<(grid_size[0]-1);i++)
		f[i]=f[i]*2.0/((data_size[0]-neg_count)*(grid[i+1]-grid[i-1]));
	if((grid[grid_size[0]-1]+0.5*(grid[grid_size[0]-1]-grid[grid_size[0]-2]))>1)
		f[grid_size[0]-1]=f[grid_size[0]-1]/((data_size[0]-neg_count)*(1-0.5*(grid[grid_size[0]-1]+grid[grid_size[0]-2])));	
	else
		f[grid_size[0]-1]=f[grid_size[0]-1]/((data_size[0]-neg_count)*(grid[grid_size[0]-1]-grid[grid_size[0]-2]));	
}

extern void get_grids_perm(char **file,
					int *columns,
					double *perm_tests,				
					double *rows,
					int *stat_columns,
					int *stat_columns_array,
					int *grid_size,
					double *grid_unif,
					double *grid_log,
					double *F_log,
					double *f_log,
					double *F_unif,
					double *f_unif,
					int *zero_count,
					int *print_messages,
					int *rownames_flag)
					
{	
	/* This function loads the p-value distribution in a uniform and a logarithmic grid*/
	/* This function is used when the permuted data is given as a text file and NOT as a data frame*/
	/* ... and the full data is also given as a text file*/
	
	/* It follows the same structure of get_grids_full with two differences:*/
	/* 	1. The grids have been already defined (so grids are an input here)*/
	/*  2. We can have several columns in the file with valid data (permutations)*/
	
	//Local variables (unnecessary)
	int columns_array[stat_columns[0]],i;
		
	for(i=0;i<stat_columns[0];i++)
		columns_array[i]=stat_columns_array[i]-1;
		
	double min_stat=grid_unif[0], max_stat=grid_unif[grid_size[0]-1];
	double delta=(double)(max_stat-min_stat)/(grid_size[0]-1);
	double log_max_ratio=log10(max_stat/min_stat);
	double delta_log=log_max_ratio/(grid_size[0]-1);
		
	//1. Initialize distributions
	for(i=0;i<grid_size[0];i++)
	{
		f_unif[i]=0;
		f_log[i]=0;
		F_unif[i]=0;
		F_log[i]=0;
	}
	zero_count[0]=0;
			
	FILE *fp;
    int j;
    int line_count=1;
    double progress=0;
    double stat_value, log_stat_value;
	char word[1000], header[10000];
	int bin,next_column;
	double stat_count_fu=0,stat_count_flog=0;
	double decimal_part;
				
	fp=fopen(*file,"r");
	for(i=0;i<(columns[0]-rownames_flag[0]);i++)
		fscanf(fp,"%s ",header);
	
	//2. Load distributions
	while(feof(fp)==0)
	{
		//Counter of the valid columns that have been read.
		next_column=0;
		for(j=0;j<columns[0];j++)
		{
			//Check if the column to read corresponds to a valid p-value column
			if(j==columns_array[next_column])
			{
				fscanf(fp,"%lf ",&stat_value);
				
				next_column++;
									
				if(stat_value==0) zero_count[0]++;					
					
				log_stat_value=log10(stat_value/min_stat);
			
				//2.1 CDF in Uniform and Logarithmic Grids
				if((min_stat<stat_value)&&(stat_value<max_stat))
				{
					bin=(int)((stat_value-min_stat)/delta);
					decimal_part=((stat_value-min_stat)/delta)-bin;
					if(decimal_part!=0)
						bin++;
					F_unif[bin]++;
					
					if(bin>=grid_size[0]) printf("Error bin uniform\n");
				
					bin=(int)(log_stat_value/delta_log);
					decimal_part=log_stat_value/delta_log-bin;
					if(decimal_part!=0)
						bin++;
					
					F_log[bin]++;
					
					if(bin>=grid_size[0]) printf("Error bin log\n");
				}
				else
				{
					if(stat_value<=min_stat)
					{							
						F_unif[0]++;
						F_log[0]++;
					}
					if(stat_value==max_stat)
					{
						F_unif[grid_size[0]-1]++;
						F_log[grid_size[0]-1]++;
					}
				}
			
				//2.2 PDF in Uniform Grid
				if(((min_stat-delta*0.5)<stat_value)&&(stat_value<(max_stat-0.5*delta)))
				{
					bin=(int)((stat_value-min_stat+0.5*delta)/delta);
					f_unif[bin]++;
					
					if(bin>=grid_size[0]) printf("Error bin uniform PDF\n");
				}
				else
				{
					if(((max_stat-delta*0.5)<stat_value)&&(stat_value<(max_stat+0.5*delta)))
						f_unif[grid_size[0]-1]++;
					else
						stat_count_fu--;
				}			
			
				//2.3 PDF in Logarithmic Grid
				if(((-delta_log*0.5)<log_stat_value)&&(log_stat_value<(log_max_ratio-0.5*delta)))
				{
					bin=(int)((log10(stat_value/min_stat)+0.5*delta_log)/delta_log);
					f_log[bin]++;
					
					if(bin>=grid_size[0]) printf("Error bin log PDF\n");
				}
				else
				{
					if(((log_max_ratio-delta_log*0.5)<log_stat_value)&&(log_stat_value<(log_max_ratio+0.5*delta_log)))
						f_log[grid_size[0]-1]++;
					else
						stat_count_flog--;
				}				
						
				//2.4 Update Progress
				if((line_count%10000000==0)&&(print_messages[0]==1))
				{
					progress=progress+(double)line_count*100/rows[0];
					line_count=0;
					printf("\nProgress: %.4lf %%",progress);
					fflush(stdout);
				}
        
				line_count++;				
			}
			else
				fscanf(fp,"%s ",word);
		}        
	}
	if(print_messages[0]==1)
		printf("\nProgress: 100.00 %%\n");
	fclose(fp);
	
	//3. Normalizing distributions
	double F_unif_acum=0,F_log_acum=0;
		
	//3.1 CDFs
	for(i=0;i<grid_size[0];i++)
	{
		F_unif_acum+=F_unif[i]/perm_tests[0];
		F_unif[i]=F_unif_acum;
		F_log_acum+=F_log[i]/perm_tests[0];
		F_log[i]=F_log_acum;		
	}
	
	//3.2 Uniform PDF
	double area_unif=0;
	if((grid_unif[0]-0.5*delta)<0)
	{
		f_unif[0]=f_unif[0]/((perm_tests[0]-stat_count_fu)*(grid_unif[0]+0.5*delta));
		area_unif+=f_unif[0]*(grid_unif[0]+0.5*delta);
	}
	else
	{
		f_unif[0]=f_unif[0]/((perm_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[0]*delta;	
	}
	for(i=1;i<(grid_size[0]-1);i++)
	{
		f_unif[i]=f_unif[i]/((perm_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[i]*delta;	
	}
	if((grid_unif[grid_size[0]-1]+0.5*delta)>1)
	{
		f_unif[grid_size[0]-1]=f_unif[grid_size[0]-1]/((perm_tests[0]-stat_count_fu)*(1-(grid_unif[grid_size[0]-1]-0.5*delta)));
		area_unif+=f_unif[grid_size[0]-1]*(1-(grid_unif[grid_size[0]-1]-0.5*delta));	
	}
	else
	{
		f_unif[grid_size[0]-1]=f_unif[grid_size[0]-1]/((perm_tests[0]-stat_count_fu)*delta);
		area_unif+=f_unif[grid_size[0]-1]*delta;	
	}
	//This is redundant
	area_unif=(1.0/area_unif);
	for(i=0;i<grid_size[0];i++)
		f_unif[i]=f_unif[i]*area_unif;	
			
	//3.3 Logarithmic PDF
	f_log[0]=(f_log[0]/(perm_tests[0]-stat_count_flog))/(grid_log[1]-grid_log[0]);
	for(i=1;i<grid_size[0]-1;i++)
		f_log[i]=(f_log[i]/(perm_tests[0]-stat_count_flog))*2/(grid_log[i+1]-grid_log[i-1]);
	f_log[grid_size[0]-1]=(f_log[grid_size[0]-1]/(perm_tests[0]-stat_count_flog))/(grid_log[grid_size[0]-1]-grid_log[grid_size[0]-2]);
}

extern void get_grids_perm_from_data(char **file,
					int *columns,	
					double *perm_tests,			
					double *rows,
					int *stat_columns,
					int *stat_columns_array,
					int *grid_size,
					double *grid,
					double *F_o,
					double *f_o,
					int *zero_count,
					int *print_messages,
					int *rownames_flag)
{
	/* This function loads the distribution of permuted p-values*/
	/* .. when the permuted data is given by text file */
	/* ..but the full data has been given as data frame*/
	
	/* Main difference with get_grids_perm is that we calculate the distributions in only one grid (data grid)*/
	/* .. instead of two grids (uniform and logarithmic)*/
	
	/* The algorithm to build the distributions is more complex than previous cases for two reasons:*/
	/* 		1. Grid is variable*/
	/*		2. Data are not sorted*/
	/* So, in order to save some time, first we build an auxilary grid*/
	/* ... and from there we can find the correct bin, without exploring the whole grid*/
	
	//Local variables (unnecessary)
	int columns_array[stat_columns[0]],i;
	
	for(i=0;i<stat_columns[0];i++)
		columns_array[i]=stat_columns_array[i]-1;
	
	double grid_aux[grid_size[0]];
	double delta=1.0/grid_size[0];
	int j, first_value[grid_size[0]];
		
	//1. Construction of an auxiliar uniform grid between 0 and 1
	for(i=0;i<grid_size[0];i++)
	{
		grid_aux[i]=i*delta;
		F_o[i]=0;
		f_o[i]=0;
	}
		
	//2. Array of correspondence with the data grid
	first_value[0]=0;
	j=0;
	for(i=1;i<grid_size[0];i++)
	{		
		while((grid_aux[i]>grid[j])&&(j<grid_size[0]))
			j++;
			
		first_value[i]=j;
	}
	
	
	FILE *fp;
	int line_count=1;
	double progress=0;
    double stat_value;
    char word[1000], header[10000];
	int bin_aux,bin,next_column;
	double f_count=0;
	int flag;
	
	fp=fopen(*file,"r");
	for(i=0;i<(columns[0]-rownames_flag[0]);i++)
		fscanf(fp,"%s ",header);
	
	//3. Load the distribution
	while(feof(fp)==0)
	{
		next_column=0;
		for(j=0;j<columns[0];j++)
		{
			if(j==columns_array[next_column])
			{
				fscanf(fp,"%lf ",&stat_value);
				next_column++;
				
				if(stat_value==0) zero_count[0]++;
				if(stat_value<(grid[0]-0.5*(grid[1]-grid[0]))) f_count--;
				if(stat_value>(grid[grid_size[0]-1]+0.5*(grid[grid_size[0]-1]-grid[grid_size[0]-2]))) f_count--;
					
				//3.1 Obtain the auxiliar bin
				bin_aux=(int)stat_value/delta;
				bin=first_value[bin_aux];
					
				//3.2 Look for data bin, and charge distributions
				flag=0;
				while((flag==0)&&(bin<grid_size[0]))
				{
					if(stat_value<=grid[bin])
					{
						F_o[bin]++;
						while((flag==0)&&(bin>0))
						{
							if(((grid[bin]-0.5*(grid[bin]-grid[bin-1]))<stat_value)&&((grid[bin]+0.5*(grid[bin+1]-grid[bin]))>stat_value))
							{
								f_o[bin]++;
								flag=1;
							}
							else
								bin--;
						}
						if(bin==0)
							flag=1;							
					}
					else
						bin++;
				}
					
				//3.3 Update Progress
				if((line_count%10000000==0)&&(print_messages[0]==1))
				{
					progress=progress+(double)line_count*100/rows[0];
					line_count=0;
					printf("\nProgress: %.4lf %%",progress);
					fflush(stdout);
				}
        
				line_count++;
				
			}
			else
				fscanf(fp,"%s ",word);
		}		
	}
	if(print_messages[0]==1)
		printf("\nProgress: 100.00 %%\n");
	fclose(fp);
	
	//4. Normalizing PDF
	if((grid[0]-0.5*(grid[1]-grid[0]))<0)
		f_o[0]=f_o[0]*2.0/((perm_tests[0]-f_count)*(grid[0]+0.5*(grid[1]-grid[0])));
	else
		f_o[0]=f_o[0]*2.0/((perm_tests[0]-f_count)*(grid[1]-grid[0]));
	for(i=0;i<(grid_size[0]-1);i++)
		f_o[i]=f_o[i]*2.0/((perm_tests[0]-f_count)*(grid[i+1]-grid[i-1]));
	if((grid[grid_size[0]-1]+0.5*(grid[grid_size[0]-1]-grid[grid_size[0]-2]))>1)
		f_o[grid_size[0]-1]=f_o[grid_size[0]-1]*2.0/((perm_tests[0]-f_count)*(1-0.5*(grid[grid_size[0]-1]+grid[grid_size[0]-2])));	
	else
		f_o[grid_size[0]-1]=f_o[grid_size[0]-1]*2.0/((perm_tests[0]-f_count)*(grid[grid_size[0]-1]-grid[grid_size[0]-2]));
		
	
	//5. Correcting and normalizing CDF
	double acum=0;
	for(i=0;i<grid_size[0];i++)
	{
		acum+=F_o[i]/perm_tests[0];
		F_o[i]=acum;
	}	
}

extern void get_grids_perm_df(double *perm_stats_array,
						int *perm_tests,
						int *grid_size,
						double *grid_log,
						double *grid_unif,
						double *F_o_log,
						double *f_o_log,
						double *F_o_unif,
						double *f_o_unif)
{
	/* This function loads the permuted p-value distributions when...*/
	/* ... permuted data is given by data frame*/
	/* ... full data is given by text file (Uniform and Logarithmic Grid)*/
	
	/* We will take advantage of the fact that permuted p-values will be in crecent order*/
	
	int i;
	int bin=0;
	
	//1. Initializing distributions
	for(i=0;i<grid_size[0];i++)
	{
		F_o_log[i]=0;
		f_o_log[i]=0;
		F_o_unif[i]=0;
		f_o_unif[i]=0;		
	}	
	
	//2. Loading CDF in Logarithmic Grid
	for(i=0;((i<perm_tests[0])&&(bin<grid_size[0]));i++)
	{
		if(perm_stats_array[i]<=grid_log[bin])
			F_o_log[bin]++;
		else
		{
			bin++;
			i--;
		}
	}
	
	//3. Loading CDF in Uniform Grid
	int aux_count=0;
	bin=0;
	for(i=0;((i<perm_tests[0])&&(bin<grid_size[0]));i++)
	{
		if(perm_stats_array[i]<=grid_unif[bin])
		{
			F_o_unif[bin]++;
			aux_count++;
		}
		else
		{
			bin++;
			i--;
		}
	}
	
	//4. Loading PDF in Uniform Grid
	int count_f_unif=0;
	double delta=grid_unif[1]-grid_unif[0];
	bin=0;
	for(i=0;perm_stats_array[i]<(grid_unif[0]-delta*0.5);i++) 
		count_f_unif++;
	for(;((i<perm_tests[0])&&(bin<(grid_size[0]-1)));i++)
	{
		if((perm_stats_array[i]>(grid_unif[bin]-delta*0.5))&&((perm_stats_array[i]<(grid_unif[bin]+delta*0.5))))
			f_o_unif[bin]++;		
		else
		{
			bin++;
			i--;
		}
	}
	bin=grid_size[0]-1;
	for(;i<perm_tests[0];i++)
	{
		if(perm_stats_array[i]<grid_unif[bin])
		{
			f_o_unif[bin]++;
		}
		else
		{
			i--;
			break;
		}
	}
	for(;i<perm_tests[0];i++)
		count_f_unif++;
	
	//5. Loading PDF in Logarithmic Grid
	int count_f_log=0;
	double delta_left=0,delta_right=grid_log[1]-grid_log[0];
	bin=0;
	for(i=0;perm_stats_array[i]<grid_log[0];i++) 
		count_f_log++;
	for(;((i<perm_tests[0])&&(bin<(grid_size[0]-1)));i++)
	{		
		if((perm_stats_array[i]>(grid_log[bin]-delta_left*0.5))&&((perm_stats_array[i]<(grid_log[bin]+delta_right*0.5))))
			f_o_log[bin]++;		
		else
		{
			bin++;
			i--;
			delta_left=grid_log[bin]-grid_log[bin-1];
			delta_right=grid_log[bin+1]-grid_log[bin];
		}
	}
	bin=grid_size[0]-1;
	for(;i<perm_tests[0];i++)
	{
		if(perm_stats_array[i]<grid_log[bin])
		{
			f_o_log[bin]++;
		}
		else
			break;
	}
	for(;i<perm_tests[0];i++)
		count_f_log++;
	
	//6. Correcting and normalizing CDFs
	double acum=0;
	for(i=0;i<grid_size[0];i++)
	{
		acum+=F_o_log[i]/perm_tests[0];
		F_o_log[i]=acum;
	}
	acum=0;
	for(i=0;i<grid_size[0];i++)
	{
		acum+=F_o_unif[i]/perm_tests[0];
		F_o_unif[i]=acum;
	}
	
	//7. Normalizing PDFs
	f_o_log[0]=f_o_log[0]*2.0/((perm_tests[0]-count_f_log)*(grid_log[1]-grid_log[0]));
	for(i=0;i<(grid_size[0]-1);i++)
		f_o_log[i]=f_o_log[i]*2.0/((perm_tests[0]-count_f_log)*(grid_log[i+1]-grid_log[i-1]));
	f_o_log[grid_size[0]-1]=f_o_log[grid_size[0]-1]*2.0/((perm_tests[0]-count_f_log)*(grid_log[grid_size[0]-1]-grid_log[grid_size[0]-2]));
	
	double area_unif=0;
	if((grid_unif[0]-0.5*delta)<0)
	{
		f_o_unif[0]=f_o_unif[0]/((perm_tests[0]-count_f_unif)*(grid_unif[0]+0.5*delta));
		area_unif+=f_o_unif[0]*(grid_unif[0]+0.5*delta);
	}
	else
	{
		f_o_unif[0]=f_o_unif[0]/((perm_tests[0]-count_f_unif)*delta);
		area_unif+=f_o_unif[0]*delta;	
	}
	for(i=1;i<(grid_size[0]-1);i++)
	{
		f_o_unif[i]=f_o_unif[i]/((perm_tests[0]-count_f_unif)*delta);
		area_unif+=f_o_unif[i]*delta;	
	}
	if((grid_unif[grid_size[0]-1]+0.5*delta)>1)
	{
		f_o_unif[grid_size[0]-1]=f_o_unif[grid_size[0]-1]/((perm_tests[0]-count_f_unif)*(1-(grid_unif[grid_size[0]-1]-0.5*delta)));
		area_unif+=f_o_unif[grid_size[0]-1]*(1-(grid_unif[grid_size[0]-1]-0.5*delta));	
	}
	else
	{
		f_o_unif[grid_size[0]-1]=f_o_unif[grid_size[0]-1]/((perm_tests[0]-count_f_unif)*delta);
		area_unif+=f_o_unif[grid_size[0]-1]*delta;	
	}
	//This is redundant
	area_unif=(1.0/area_unif);
	for(i=0;i<grid_size[0];i++)
		f_o_unif[i]=f_o_unif[i]*area_unif;	
}


extern void interpolate_fdrs(
					  char **file, 
					  int *columns,
                      int *flag_rownames,
					  double *rows,
					  int *stat_column, 
					  int *sampled_stats,
					  double *sampled_stats_array, 
					  double *Fdr, 
					  double *Fndr, 
					  double *Fdr_ST, 
					  double *Fndr_ST, 
					  double *Fdr_BH_unif, 
					  double *Fdr_BH_perm, 
					  double *Fdr_MV,
                      double *BF,
					  char **name,
					  int *alt_methods,
					  double *zero_Fdr,
                      double *report_threshold,
                      int *print_messages,
                      double *significance_threshold,
                      double *significant_tests)
{
	/* This function interpolates the fdrs at each p-values from the fdrs at the level of grids*/
	
	/* In the case of p-values given by text file, we calculate fdrs in a grid (uniform and logarithmic)*/
	/* in order to report a fdr for each individual p-value, we need this function*/
	
	int flag;
	int column_stat=stat_column[0]-1;
	
	int i;
	//Length of alt_methods vector: Fdr_ST, Fndr_ST,Fdr_BH_unif,Fdr_BH_perm,Fdr_MV,BF:
	int methods=6;
	char **method_names;
    //We start putting always fdr and Fndr: thats why correct_methods=2
	int correct_methods=2;
	double **Fdr_local,*Fdr_value,*correct_zero_fdr;
	FILE *fp,*fp2;
    int j;
    int line_count=1;
    double progress=0;
    double stat_value;
	char word[columns[0]][1000];
	char output_file[1000];
	int bin_log;
	double min_stat=sampled_stats_array[0], max_stat=sampled_stats_array[sampled_stats[0]-1];
	double log_max_ratio=log10(max_stat/min_stat);
	double delta_log=log_max_ratio/(sampled_stats[0]-1);
    int correct_fdrs=0;
    
    //1. We make fdr and fndr monotone
    bounded_monotone_fdr(Fdr,sampled_stats[0]);
    bounded_monotone_fndr(Fndr,sampled_stats[0]);
    
    //2. We load fdrs and fndrs from other methods if requested by user
	for(i=0;i<methods;i++)
		if(alt_methods[i]==1) correct_methods++;
	
	//2.1 We save memory
    int fdr_indexes[6];    
	Fdr_local=(double **)malloc(correct_methods*sizeof(double *));
	Fdr_value=(double *)malloc(correct_methods*sizeof(double *));
	method_names=(char **)malloc(correct_methods*sizeof(char *));
    correct_zero_fdr=(double *)malloc(correct_methods*sizeof(double *));

	for(i=0;i<correct_methods;i++)
	{
		Fdr_local[i]=(double *)malloc(sampled_stats[0]*sizeof(double));
		method_names[i]=(char *)malloc(100*sizeof(char));
	}

	int correct_methods_2=0;
	//The first two -Fdr / Fndr- are always processed:
    
    //2.2. Load Fdr 
    sprintf(method_names[correct_methods_2]," Fdr_SAB");
    correct_zero_fdr[correct_methods_2]=zero_Fdr[0];
    for(i=0;i<sampled_stats[0];i++)
        Fdr_local[correct_methods_2][i]=Fdr[i];
    correct_methods_2++;
    fdr_indexes[correct_fdrs]=0;
    significant_tests[correct_methods_2]=0;
    correct_fdrs++;
    
	//2.3 Load Fndr
    sprintf(method_names[correct_methods_2]," Fndr_SAB");
    correct_zero_fdr[correct_methods_2]=zero_Fdr[1];
    for(i=0;i<sampled_stats[0];i++)
        Fdr_local[correct_methods_2][i]=Fndr[i];
    significant_tests[correct_methods_2]=0;
    correct_methods_2++;

	//2.4 Load Fdr Storey-Tibshirani if requested
	if(alt_methods[0]==1)
	{
		printf("Load FdrST\n");
		sprintf(method_names[correct_methods_2]," Fdr_ST");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[2];
        bounded_monotone_fdr(Fdr_ST,sampled_stats[0]);
		for(i=0;i<sampled_stats[0];i++)
			Fdr_local[correct_methods_2][i]=Fdr_ST[i];
        fdr_indexes[correct_fdrs]=correct_methods_2;
        significant_tests[correct_methods_2]=0;
        correct_fdrs++;
        correct_methods_2++;
	}
    
    //2.4 Load Fndr Storey-Tibshirani if requested
	if(alt_methods[1]==1)
	{
		printf("Load FndrST\n");
		sprintf(method_names[correct_methods_2]," Fndr_ST");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[3];
        bounded_monotone_fndr(Fndr_ST,sampled_stats[0]);
		for(i=0;i<sampled_stats[0];i++)
			Fdr_local[correct_methods_2][i]=Fndr_ST[i];
		significant_tests[correct_methods_2]=0;
		correct_methods_2++;
	}

    //2.4 Load Fdr Benjamini-Hochberg (uniform null distribution) if requested
	if(alt_methods[2]==1)
	{
		printf("Load FdrBHunif\n");
		sprintf(method_names[correct_methods_2]," Fdr_BH_unif");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[4];
        bounded_monotone_fdr(Fdr_BH_unif,sampled_stats[0]);
		for(i=0;i<sampled_stats[0];i++)
			Fdr_local[correct_methods_2][i]=Fdr_BH_unif[i];
        fdr_indexes[correct_fdrs]=correct_methods_2;
        significant_tests[correct_methods_2]=0;
        correct_fdrs++;
        correct_methods_2++;
	}

    //2.5 Load Fdr Benjamini-Hochberg (permutations null distribution) if requested
	if(alt_methods[3]==1)
	{
		printf("Load FdrBHperm\n");
		sprintf(method_names[correct_methods_2]," Fdr_BH_perm");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[5];
        bounded_monotone_fdr(Fdr_BH_perm,sampled_stats[0]);
		for(i=0;i<sampled_stats[0];i++)
			Fdr_local[correct_methods_2][i]=Fdr_BH_perm[i];
        fdr_indexes[correct_fdrs]=correct_methods_2;
        significant_tests[correct_methods_2]=0;
        correct_fdrs++;
        correct_methods_2++;
	}
	
	//2.6 Load Fdr Millstein and Volfson if requested
	if(alt_methods[4]==1)
	{
		printf("Load FdrMV\n");
		sprintf(method_names[correct_methods_2]," Fdr_MV");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[6];
        bounded_monotone_fdr(Fdr_MV,sampled_stats[0]);
		for(i=0;i<sampled_stats[0];i++)
			Fdr_local[correct_methods_2][i]=Fdr_MV[i];
        fdr_indexes[correct_fdrs]=correct_methods_2;
        significant_tests[correct_methods_2]=0;
        correct_fdrs++;
        correct_methods_2++;
	}
    
    //2.7 Load Fdr Bonferroni if requested
    if(alt_methods[5]==1)
    {
		printf("Load BF\n");
        sprintf(method_names[correct_methods_2]," BF");
        correct_zero_fdr[correct_methods_2]=zero_Fdr[7];
        //Bonferroni is monotone and bounded by construction
        for(i=0;i<sampled_stats[0];i++)
            Fdr_local[correct_methods_2][i]=BF[i];
        fdr_indexes[correct_fdrs]=correct_methods_2;
        significant_tests[correct_methods_2]=0;
        correct_fdrs++;
        correct_methods_2++;
    }   
    
    //3. Read the original file, and write another one with additional columns corresponding to Fdrs and Fndrs
    fp=fopen(*file,"r");
	//sprintf(command,"mkdir -p outputs_fdrs_%s",*name);
	//system(command);
	sprintf(output_file,"%s/fdrs.txt",*name);
	fp2=fopen(output_file,"w");
	    
	//3.1 Read and write the headers
	for(j=0;j<(columns[0]-flag_rownames[0]);j++)
	{
		fscanf(fp,"%s ",word[j]);
		fprintf(fp2,"%s ",word[j]);
	}
    for(j=0;j<correct_methods;j++)
        fprintf(fp2,"%s ",method_names[j]);
    fprintf(fp2,"\n");
	    
    while(feof(fp)==0)
	{
		//3.2 Read columns BEFORE p-value
		for(j=0;j<column_stat;j++)
        	fscanf(fp,"%s ",word[j]);
        		
        //3.3 Read p-value
		fscanf(fp,"%lf ",&stat_value);
					
		//3.4 Read columns AFTER p-value
		for(j=(column_stat+1);j<columns[0];j++)
			fscanf(fp,"%s ",word[j]);
		
		//3.5 Log transform
		bin_log=(int)((log10(stat_value/min_stat))/delta_log);

		//3.6 Interpolate Fdrs and Fndrs for every different method
		if(bin_log<0)
		{			
			//Fdr at zero is calculated apart from the rest
			// Is the same at every method? Dont think so... 
			// Should be different for Fndr
            // Joaquin: Corrected, now I'm passing a vector to the function instead the (wrong) escalar from before
			for(j=0;j<correct_methods;j++)
				Fdr_value[j]=zero_Fdr[j];
		}
		else
		{
			for(j=0;j<correct_methods;j++)
			{
				if(bin_log<(sampled_stats[0]-1))
					Fdr_value[j]=Fdr_local[j][bin_log]+(Fdr_local[j][bin_log+1]-Fdr_local[j][bin_log])*(stat_value-sampled_stats_array[bin_log])/(sampled_stats_array[bin_log+1]-sampled_stats_array[bin_log]);
				else
					Fdr_value[j]=Fdr_local[j][bin_log];							
			}
		}				
            
        //3.7 Check if the Fdr is below the report threshold for any method
        flag=0;
        for(j=0;j<correct_fdrs;j++)
        {
            if(Fdr_value[fdr_indexes[j]]<=report_threshold[0])
            {
                flag=1;
                break;
            }            
        }
        
        //3.8 Update counter of significant tests
        for(j=0;j<correct_methods;j++)
        {
			if(Fdr_value[j]<=significance_threshold[0])
				significant_tests[j]++;
		}
			
		//3.9 Write the data (if fdr is below report threshold for any method)
		if(flag==1)
		{
			for(j=0;j<column_stat;j++)
				fprintf(fp2,"%s ",word[j]);
			fprintf(fp2,"%.10e ",stat_value);
			for(j=(column_stat+1);j<columns[0];j++)
				fprintf(fp2,"%s ",word[j]);
			for(j=0;j<correct_methods;j++)
				fprintf(fp2,"%.10lf ",Fdr_value[j]);
			fprintf(fp2,"\n");
		}	
				
		//3.10 Update progress
		if((line_count%10000000==0)&&(print_messages[0]==1))
        {
            progress=progress+(double)line_count*100/rows[0];
            line_count=0;
            printf("\nProgress: %.4lf %%",progress);
            fflush(stdout);
        }
		
		line_count++;
	}
	if(print_messages[0]==1)
	{
		printf("\nProgress: 100.00 %%\n");
		fflush(stdout);
	}
	fclose(fp);
	fclose(fp2);
	free(Fdr_local);	
}

extern void monotone_fdrs(double *Fdr, double *Fndr, double *Fdr_ST,double *Fndr_ST, double *Fdr_BH_unif, double *Fdr_BH_perm, double *Fdr_MV, int *alt_methods,int *sampled_stats)
{
	/* This function monotone Fdrs and Fndrs for the different reqested methods*/
	/* We will only need this (instead of interpolate_fdrs) when data is given via data frame*/
	
	bounded_monotone_fdr(Fdr,sampled_stats[0]);
    bounded_monotone_fndr(Fndr,sampled_stats[0]);
    if(alt_methods[0]==1)  bounded_monotone_fdr(Fdr_ST,sampled_stats[0]);
    if(alt_methods[1]==1)  bounded_monotone_fndr(Fndr_ST,sampled_stats[0]);
    if(alt_methods[2]==1)  bounded_monotone_fdr(Fdr_BH_unif,sampled_stats[0]);
    if(alt_methods[3]==1)  bounded_monotone_fdr(Fdr_BH_perm,sampled_stats[0]);
    if(alt_methods[4]==1)  bounded_monotone_fdr(Fdr_MV,sampled_stats[0]);
}

void bounded_monotone_fdr(double *vector,int size)
{
    int i,j;
    
    if(vector[0]<0)
        vector[0]=0;
    if(vector[0]>1)
        vector[0]=1;
    
    for(i=1;i<size;i++)
    {
        if(vector[i]<0)
            vector[i]=0;
        if(vector[i]>1)
            vector[i]=1;
        for(j=1;j<=(i);j++)
        {
            if(vector[i-j]>vector[i])
                vector[i-j]=vector[i];
            else
                break;
        }
    }
}

void bounded_monotone_fndr(double *vector,int size)
{
    int i,j;
    
    if(vector[0]<0)
        vector[0]=0;
    if(vector[0]>1)
        vector[0]=1;
    
    for(i=(size-1);i>=0;i--)
    {
        if(vector[i]<0)
            vector[i]=0;
        if(vector[i]>1)
            vector[i]=1;
        for(j=1;j<(size-i);j++)
        {
            if(vector[i+j]>vector[i])
                vector[i+j]=vector[i];
            else
                break;
        }
    }
}

int least_concave_majorant(double *x,double *F, int first_point, int final_point, double *F_g, double *f_g, int size, double *f)
{
	double slope,slope_aux;
	double initial_x;
	double initial_y;
	int	initial_point=first_point,last_point,i;
	
	if(global_flag==1)
		printf("initial_point: %d \n",initial_point);
	if(initial_point>0)
	{
		initial_x=x[initial_point-1];
		initial_y=F[initial_point-1];
	}
	else
	{
		initial_x=x[initial_point];
		initial_y=F[initial_point];
	}
		
	while(initial_point<=final_point)
	{
		if(initial_point==0)
			slope=f[0];
		else
			slope=(F[initial_point]-initial_y)/(x[initial_point]-initial_x);
				
		last_point=initial_point;
		for(i=initial_point+1;i<=final_point;i++)
		{
			slope_aux=(F[i]-initial_y)/(x[i]-initial_x);
						
			if(i>=size)
			{
				printf("Problem least_concave_majorant slope aux\n");
				printf("initial_point %d  final_point %d  size %d\n", initial_point,final_point,size);
				printf("i %d \n\n",i);
			}
			
			if(slope_aux>slope)
			{
				slope=slope_aux;
				last_point=i;					
			}
		}
		for(i=initial_point;i<=last_point;i++) 
		{
			F_g[i]=initial_y+slope*(x[i]-initial_x);
			f_g[i]=slope;
			
			if(i>=size)
			{
				printf("Problem least_concave_majorant slope \n");
				printf("initial_point %d  final_point %d  size %d\n", initial_point,final_point,size);
				printf("i %d \n\n",i);
			}
		}
		initial_point=last_point+1;
		initial_y=F[last_point];
		initial_x=x[last_point];
	}
	
	/*if(flag_ini_zero==1)
	{
		F_g[0]=0;
		f_g[0]=f_g[1];
	}*/
	
	if(global_flag==1)
		printf("First point: F=%lf f=%lf\n",F_g[0],f_g[0]);
	
	return initial_point;
}

int greatest_convex_minorant(double *x,double *F, int first_point, int final_point, double *F_g, double *f_g, int size, double *f)
{
	double slope,slope_aux;
	int initial_point=first_point,last_point,i;
	double initial_y;
	double	initial_x;	
	
	if(initial_point>0)
	{
		initial_x=x[initial_point-1];
		initial_y=F[initial_point-1];
	}
	else
	{
		initial_x=x[initial_point];
		initial_y=F[initial_point];
		
		/*if(x[0]==0)
		{
			flag_ini_zero=1;
			initial_point=1;
		}*/
	}
	
	while(initial_point<=final_point)
	{			
		if(initial_point==0)
			slope=f[0];
		else
			slope=(F[initial_point]-initial_y)/(x[initial_point]-initial_x);
			
		last_point=initial_point;
		
		for(i=initial_point+1;i<=final_point;i++)
		{
			slope_aux=(F[i]-initial_y)/(x[i]-initial_x);
			
			if(i>=size)
			{
				printf("Problem greatest_convex_minorant slope aux\n");
				printf("initial_point %d  final_point %d  size %d\n", initial_point,final_point,size);
				printf("i %d \n\n",i);
			}
					
			if(slope_aux<slope)
			{
				slope=slope_aux;
				last_point=i;				
			}				
		}
		for(i=initial_point;i<=last_point;i++)
		{
			F_g[i]=initial_y+slope*(x[i]-initial_x);
			f_g[i]=slope;
			
			if(i>=size)
			{
				printf("Problem greatest_convex_minorant slope aux\n");
				printf("initial_point %d  final_point %d  size %d\n", initial_point,final_point,size);
				printf("i %d \n\n",i);
			}		
		}
		initial_point=last_point+1;
		initial_y=F[last_point];
		initial_x=x[last_point];
	}
	
	if(global_flag==1)
		printf("Last point: F=%lf f=%lf\n",F_g[last_point],f_g[last_point]);
	
	/*if(flag_ini_zero==1)
	{
		F_g[0]=0;
		f_g[0]=f_g[1];
	}*/
	
	return initial_point;
}

double grenander_error(double *F, double *F_g, double *x, int size)
{
	double error=0;
	int i;
		
	error+=(F[0]-F_g[0])*(-1+2*((F[0]-F_g[0])>0))*(x[1]-x[0]);
	for(i=1;i<(size-1);i++)
		error+=2*(F[i]-F_g[i])*(-1+2*((F[i]-F_g[i])>0))*(x[i+1]-x[i-1]);	
	error+=(F[size-1]-F_g[size-1])*(-1+2*((F[size-1]-F_g[size-1])>0))*(x[size-1]-x[size-2]);	
	
	return error;
}

extern void estimate_p_star(double *x,
					int *grid_size,
					double *F, 
					double *F_o, 
					double *f, 
					double *f_o, 
					double *F_g, 
					double *F_o_g, 
					double *f_g,
					double *f_o_g,
					int *index_p_star,
					double *fitting_goodness,
					int *flag_grenander_choice)
{
	
	double goodness;
	int i,flag_good_fit=1,index_p_prime[2];

    index_p_prime[0]=-1;
    index_p_prime[1]=-1;

	//1. Fit full data
	goodness=choose_fit(x,grid_size[0],F,f,F_g,f_g,0,flag_grenander_choice,index_p_prime);	
	if(goodness<0)
		flag_good_fit=0;
	fitting_goodness[0]=goodness;
	
	//2. Fit perm data
	goodness=choose_fit(x,grid_size[0],F_o,f_o,F_o_g,f_o_g,1,flag_grenander_choice,index_p_prime);
	if(goodness<0)
		flag_good_fit=0;
	fitting_goodness[1]=goodness;
	
	//3. Estimating p_star (scaling point) if R2 is good
	double min_ratio;
	min_ratio=f_g[0]/f_o_g[0];
	index_p_star[0]=0;
	for(i=1;i<grid_size[0];i++)
	{
		if(f_g[i]/f_o_g[i]<min_ratio)
		{
			min_ratio=f_g[i]/f_o_g[i];
			index_p_star[0]=i;
		}
	}
	index_p_star[0]++;
}

double choose_fit(double *x,
					int grid_size,
					double *F, 
					double *f, 
					double *F_g,  
					double *f_g,
					int flag,
					int *flag_choice,
					int *index_p_prime)
{
	int i,flag_monotone[2];
	double F_g_aux[grid_size],f_g_aux[grid_size];
	//FILE *fp;
	int p_aux_index[100],delta_int,aux_index;
	int n_small;
	int range[2];
	double x_min,x_max;
	
	//A. We agregate data of the probability density function in lesser bins
	//The rationale here, is to reduce "white noise" that comes as the consequence of
	// ... having very small counts per bin.
	//However, this new grid should contain enough bins to detect different tendencies in the distributions
	if(grid_size<5000)
		n_small=100;	
	else
	{
		if(grid_size<50000)
			n_small=(int)(grid_size/50.0);
		else
			n_small=10000;
	}	
	double f_small[n_small],F_small[n_small],x_small[n_small];	
	
	if(grid_size>100)
	{
		small_grid_calculation(f,x,grid_size,f_small,x_small,n_small);
		small_grid_calculation(F,x,grid_size,F_small,x_small,n_small);	
		x_min=x_small[3];
		x_max=x_small[n_small-4];
        calculate_range(x_min,x_max,x,range,grid_size);
	}
	
	//B. Define an array of 100 bins (this is not the "small" grid)
	if(grid_size>100)
	{
		delta_int=(int)(grid_size/99);
		for(i=0;i<99;i++)
			p_aux_index[i]=i*delta_int;
		p_aux_index[99]=grid_size-1;
	}
	
    
    //If flag_choice is zero, we need to choose the pattern: 1.\ 2.\/ 3./ or 4./\. Otherwise only the apriori know pattern is fitted
    double goodness;
    if(flag_choice[flag]==0)
    {
        //C. Fit \/
        //p_prime marks the change in concavity,
        //and is determined by minimizing the error
        int p_prime_index,p_prime_index_min,aux_index_min,initial_point,final_point;
        double error=999,error_aux,area_grenander;
        //C.1 Calculates optimum p_prime (1st level)
        //We will use the 100 values of the array previously defined, as a first aproximation of the result
        //We will skip the extremes
        if(grid_size>100)
        {
            for(aux_index=0;aux_index<99;aux_index++)
            {
                p_prime_index=p_aux_index[aux_index];
                if((x[p_prime_index]>x_min)&&(x[p_prime_index]<x_max))
                {
                    initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                        aux_index_min=aux_index;
                    }
                }
            }
        }
        //C.2. Calculates optimum p_prime (2nd level)
        //Now we look in the vicinity of the first election of p prime
        //If the number of data is low enough, we will use each value
        //Otherwise we just define another set of values (200 this time)
        int first, last;
        int p_aux_index_2[200],delta_int_2;
        if(grid_size<100)
        {
            first=3;
            last=grid_size-3;
        }
        else
        {
            first=p_aux_index[aux_index_min-1]+1;
            last=p_aux_index[aux_index_min+1];
            
            if(x[first]<x_min) first=range[0]+1;
            if(x[last]>x_max) last=range[1]-1;
        }
        if(grid_size<=10000)
        {
            for(p_prime_index=first;p_prime_index<last;p_prime_index++)
            {
                initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                error_aux=grenander_error(F,F_g,x,grid_size);
                if(error_aux<error)
                {
                    error=error_aux;
                    p_prime_index_min=p_prime_index;
                }
            }
        }
        else
        {
            delta_int_2=(int)((last-first)/199);
            for(i=0;i<199;i++)
                p_aux_index_2[i]=first+i*delta_int_2;
            p_aux_index_2[199]=last-1;
            
            for(aux_index=0;aux_index<200;aux_index++)
            {
                p_prime_index=p_aux_index_2[aux_index];
                initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                error_aux=grenander_error(F,F_g,x,grid_size);
                if(error_aux<error)
                {
                    error=error_aux;
                    p_prime_index_min=p_prime_index;
                    aux_index_min=aux_index;
                }
            }
        }
        
        //C.3 Calculating best distribution (once p_prime is determined)
        //printf("p_prime %lf\n",x[p_prime_index_min]);
        initial_point=least_concave_majorant(x,F,0,p_prime_index_min,F_g,f_g,grid_size,f);
        final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
        //printf("El ultimo: %lf\n",F_g[grid_size-1]);
        //getchar();
        //C.4 Correcting the slopes (Grenander's PDFs)
        //f_g[0]=f_g[1]; //No sé porque estaba esto, vamos a tratar el primer punto como el resto
        //for(i=1;i<(grid_size-2);i++)
        //f_g[i]=0.5*(f_g[i]+f_g[i+1]);
        correcting_slopes(f_g,x,grid_size);
        //C.5 Renormalizing Grenander's PDFs
        area_grenander=f_g[0]*(x[1]-x[0]);
        for(i=1;i<(grid_size-1);i++)
            area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
        area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
        area_grenander=1.0/area_grenander;
        for(i=0;i<grid_size;i++)
            f_g[i]=f_g[i]*area_grenander;
        
        //D. Calculates explained variance and unexplainded variance at both sides of p_prime
        double f_g_small[n_small],F_g_small[n_small];
        double error_options[2], error_grenander_choice[2];
        if(grid_size>100)
        {
            small_grid_calculation(f_g,x,grid_size,f_g_small,x_small,n_small);
            //flag_monotone[0]=explained_vs_unexplained_variance(x,x_small,p_prime_index_min,n_small,f_small,f_g_small);
            small_grid_calculation(F_g,x,grid_size,F_g_small,x_small,n_small);
            flag_monotone[0]=explained_vs_unexplained_variance(x,x_small,p_prime_index_min,n_small,F_small,F_g_small);
        }
        else
        {
            //flag_monotone[0]=explained_vs_unexplained_variance(x,x,p_prime_index_min,grid_size,f,f_g);
            flag_monotone[0]=explained_vs_unexplained_variance(x,x,p_prime_index_min,grid_size,F,F_g);
        }
        
        
        if(flag_monotone[0]==1)
        {
            //D.1 If the condition is fulfilled we have to change to a \ fitting
            initial_point=least_concave_majorant(x,F,0,grid_size-1,F_g,f_g,grid_size,f);
            
            //D.2 Correcting the slopes (Grenander's PDFs)
            //f_g[0]=f_g[1];
            //for(i=1;i<(grid_size-2);i++)
            //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            
            //D.3 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;
            
            
        }
        
        //if(flag_monotone[0]==1)
          //  printf(" -->  Monotically decreasing distribution ");
        //else
          //  printf(" -->  Valley distribution (min. at %lf) ", x[p_prime_index_min]);
        
        //E.Assigning error to first option \ or \/
        
        error_options[0]=0;
        error_options[0]=grenander_error(f,f_g,x,grid_size);
        error_grenander_choice[0]=grenander_error(F,F_g,x,grid_size);
        //printf("with error: %.12lf \n",error_grenander_choice[0]);
        
        /*DEBUG*/
        ///////////////////////
        //if(flag==0)
          //  fp=fopen("first_fit_full.txt","wt");
        //else
          //  fp=fopen("first_fit_perm.txt","wt");
        //fprintf(fp,"x f f_g F F_g\n");
        //for(i=0;i<grid_size;i++)
          //  fprintf(fp,"%.32lf %.32lf %.32lf %lf %lf\n",x[i],f[i],f_g_aux[i],F[i],F_g_aux[i]);
        //fclose(fp);
        ///////////////////////
        /*END OF DEBUG*/
        
        //F. Save f_g and F_g
        for(i=0;i<grid_size;i++)
        {
            F_g_aux[i]=F_g[i];
            f_g_aux[i]=f_g[i];
            F_g[i]=0;
            f_g[i]=0;
            index_p_prime[flag]=p_prime_index_min;
        }
        
        /* Symmetric option*/
        
        //G. Fit /\
        //p_prime marks the change in concavity,
        //and is determined by minimizing the error
        error=999;
        //G.1 Calculates optimum p_prime (1st level)
        //We will use the 100 values of the array previously defined, as a first aproximation of the result
        //We skip the extremes
        if(grid_size>100)
        {
            for(aux_index=0;aux_index<99;aux_index++)
            {
                p_prime_index=p_aux_index[aux_index];
                
                if((x[p_prime_index]>x_min)&&(x[p_prime_index]<x_max))
                {
                    
                    initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                        aux_index_min=aux_index;
                    }
                }
            }
        }
        //G.2. Calculates optimum p_prime (2nd level)
        //Now we look in the vicinity of the first election of p prime
        //If the number of data is low enough, we will use each value
        //Otherwise we just define another set of values (200 this time)
        if(grid_size<=100)
        {
            first=3;
            last=grid_size-3;
        }
        else
        {
            first=p_aux_index[aux_index_min-1]+1;
            last=p_aux_index[aux_index_min+1];
            
            if(x[first]<x_min) first=range[0];
            if(x[last]>x_max) last=range[1];
        }
        if(grid_size<=10000)
        {
            for(p_prime_index=first;p_prime_index<last;p_prime_index++)
            {
                initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                error_aux=grenander_error(F,F_g,x,grid_size);
                if(error_aux<error)
                {
                    error=error_aux;
                    p_prime_index_min=p_prime_index;
                }
            }
        }
        else
        {
            delta_int_2=(int)((last-first)/199);
            for(i=0;i<199;i++)
                p_aux_index_2[i]=first+i*delta_int_2;
            p_aux_index_2[199]=last-1;
            
            for(aux_index=0;aux_index<200;aux_index++)
            {
                p_prime_index=p_aux_index_2[aux_index];
                initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                error_aux=grenander_error(F,F_g,x,grid_size);
                if(error_aux<error)
                {
                    error=error_aux;
                    p_prime_index_min=p_prime_index;
                    aux_index_min=aux_index;
                }
            }
        }
        //G.3 Calculating best distribution (once p_prime is determined)
        //printf("p_prime %lf\n",x[p_prime_index_min]);
        initial_point=greatest_convex_minorant(x,F,0,p_prime_index_min,F_g,f_g,grid_size,f);
        final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
        //G.4 Correcting the slopes (Grenander's PDFs)
        //f_g[0]=f_g[1];
        //for(i=1;i<(grid_size-2);i++)
        //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
        correcting_slopes(f_g,x,grid_size);
        //G.5 Renormalizing Grenander's PDFs
        area_grenander=f_g[0]*(x[1]-x[0]);
        for(i=1;i<(grid_size-1);i++)
            area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
        area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
        area_grenander=1.0/area_grenander;
        for(i=0;i<grid_size;i++)
            f_g[i]=f_g[i]*area_grenander;
        
        //H. Calculates explained variance and unexplainded variance at both sides of p_prime
        if(grid_size>100)
        {
            //global_flag=1;
            small_grid_calculation(f_g,x,grid_size,f_g_small,x_small,n_small);
            //flag_monotone[1]=explained_vs_unexplained_variance(x,x_small,p_prime_index_min,n_small,f_small,f_g_small);
            small_grid_calculation(F_g,x,grid_size,F_g_small,x_small,n_small);
            flag_monotone[1]=explained_vs_unexplained_variance(x,x_small,p_prime_index_min,n_small,F_small,F_g_small);
            //global_flag=0;
            //flag_monotone[1]=0;
        }
        else
        {
            //flag_monotone[1]=explained_vs_unexplained_variance(x,x,p_prime_index_min,grid_size,f,f_g);
            flag_monotone[1]=explained_vs_unexplained_variance(x,x,p_prime_index_min,grid_size,F,F_g);
        }
        
        if(flag_monotone[1]==1) 
        {
            //H.1 If the condition is fulfilled we have to change to a / fitting
            final_point=greatest_convex_minorant(x,F,0,grid_size-1,F_g,f_g,grid_size,f);
            
            //H.2 Correcting the slopes (Grenander's PDFs) 
            //f_g[0]=f_g[1];
            //for(i=1;i<(grid_size-2);i++)
            //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            
            //H.3 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;	
        }
        
        //if(flag_monotone[1]==1)
          //  printf(" -->  Monotically increasing distribution ");
        //else
          //  printf(" -->  Mountain distribution (max. at %lf) ",x[p_prime_index_min]);
        
        //I.Assigning error to second option  
        error_options[1]=0;	
        error_options[1]=grenander_error(f,f_g,x,grid_size);
        error_grenander_choice[1]=grenander_error(F,F_g,x,grid_size);
        //printf("with error: %.12lf\n",error_grenander_choice[1]);
        
        /*DEBUG*/
        ///////////////////////
        //if(flag==0)
          //  fp=fopen("second_fit_full.txt","wt");
        //else
          //  fp=fopen("second_fit_perm.txt","wt");
        //fprintf(fp,"x f f_g F F_g\n");
        //for(i=0;i<grid_size;i++)
          //  fprintf(fp,"%e %e %e %lf %lf\n",x[i],f[i],f_g[i],F[i],F_g[i]);
        //fclose(fp);
        ///////////////////////
        /*END OF DEBUG*/
        
        //J. Choose better option
        if(error_grenander_choice[0]<error_grenander_choice[1])
        {
            //printf(" We choose first option\n");
            for(i=0;i<grid_size;i++)
            {
                F_g[i]=F_g_aux[i];
                f_g[i]=f_g_aux[i];
            }
            flag_choice[flag]=2-flag_monotone[0];
        }
        else
        {
            //printf(" We choose second option\n");
            flag_choice[flag]=4-flag_monotone[1];
        }
        
        //K. Goodness of fit
        if(((flag_monotone[0]==1)&&(error_grenander_choice[0]<error_grenander_choice[1]))||((flag_monotone[1]==1)&&(error_grenander_choice[0]>error_grenander_choice[1])))
        {
            if(grid_size>100)
                small_grid_calculation(f_g,x,grid_size,f_g_small,x_small,n_small);
        }
        
        if(grid_size>100)
            goodness=goodness_calc(f_small,f_g_small,x_small,n_small);
        else
            goodness=goodness_calc(f,f_g,x,grid_size);
        
    }else{
        int p_prime_index,p_prime_index_min,aux_index_min,initial_point,final_point;
        double error=999,error_aux,area_grenander;
        int first, last;
        int p_aux_index_2[200],delta_int_2;

        if(flag_choice[flag]==1){
            //D.1 If the condition is fulfilled we have to change to a \ fitting
            initial_point=least_concave_majorant(x,F,0,grid_size-1,F_g,f_g,grid_size,f);
            
            //D.2 Correcting the slopes (Grenander's PDFs)
            //f_g[0]=f_g[1];
            //for(i=1;i<(grid_size-2);i++)
            //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            
            //D.3 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;
            
            
        }else if(flag_choice[flag]==2)
        {
            //C. Fit \/
            //p_prime marks the change in concavity,
            //and is determined by minimizing the error
            
            //C.1 Calculates optimum p_prime (1st level)
            //We will use the 100 values of the array previously defined, as a first aproximation of the result
            //We will skip the extremes
            if(grid_size>100)
            {
                for(aux_index=0;aux_index<99;aux_index++)
                {
                    p_prime_index=p_aux_index[aux_index];
                    if((x[p_prime_index]>x_min)&&(x[p_prime_index]<x_max))
                    {
                        initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                        final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                        error_aux=grenander_error(F,F_g,x,grid_size);
                        if(error_aux<error)
                        {
                            error=error_aux;
                            p_prime_index_min=p_prime_index;
                            aux_index_min=aux_index;
                        }
                    }
                }
            }
            //C.2. Calculates optimum p_prime (2nd level)
            //Now we look in the vicinity of the first election of p prime
            //If the number of data is low enough, we will use each value
            //Otherwise we just define another set of values (200 this time)

            if(grid_size<100)
            {
                first=3;
                last=grid_size-3;
            }
            else
            {
                first=p_aux_index[aux_index_min-1]+1;
                last=p_aux_index[aux_index_min+1];
                
                if(x[first]<x_min) first=range[0]+1;
                if(x[last]>x_max) last=range[1]-1;
            }
            if(grid_size<=10000)
            {
                for(p_prime_index=first;p_prime_index<last;p_prime_index++)
                {
                    initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                    }
                }
            }
            else
            {
                delta_int_2=(int)((last-first)/199);
                for(i=0;i<199;i++)
                    p_aux_index_2[i]=first+i*delta_int_2;
                p_aux_index_2[199]=last-1;
                
                for(aux_index=0;aux_index<200;aux_index++)
                {
                    p_prime_index=p_aux_index_2[aux_index];
                    initial_point=least_concave_majorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                        aux_index_min=aux_index;
                    }
                }
            }
            
            //C.3 Calculating best distribution (once p_prime is determined)
            //printf("p_prime %lf\n",x[p_prime_index_min]);
            initial_point=least_concave_majorant(x,F,0,p_prime_index_min,F_g,f_g,grid_size,f);
            final_point=greatest_convex_minorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
            //printf("El ultimo: %lf\n",F_g[grid_size-1]);
            //getchar();
            //C.4 Correcting the slopes (Grenander's PDFs)
            //f_g[0]=f_g[1]; //No sé porque estaba esto, vamos a tratar el primer punto como el resto
            //for(i=1;i<(grid_size-2);i++)
            //f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            //C.5 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;
        }else if(flag_choice[flag]==3){
            //H.1 If the condition is fulfilled we have to change to a / fitting
            final_point=greatest_convex_minorant(x,F,0,grid_size-1,F_g,f_g,grid_size,f);
            
            //H.2 Correcting the slopes (Grenander's PDFs)
            //f_g[0]=f_g[1];
            //for(i=1;i<(grid_size-2);i++)
            //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            
            //H.3 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;
        }else if(flag_choice[flag]==4){//G. Fit /\
            //p_prime marks the change in concavity,
            //and is determined by minimizing the error
            error=999;
            //G.1 Calculates optimum p_prime (1st level)
            //We will use the 100 values of the array previously defined, as a first aproximation of the result
            //We skip the extremes
            if(grid_size>100)
            {
                for(aux_index=0;aux_index<99;aux_index++)
                {
                    p_prime_index=p_aux_index[aux_index];
                    
                    if((x[p_prime_index]>x_min)&&(x[p_prime_index]<x_max))
                    {
                        
                        initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                        final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                        error_aux=grenander_error(F,F_g,x,grid_size);
                        if(error_aux<error)
                        {
                            error=error_aux;
                            p_prime_index_min=p_prime_index;
                            aux_index_min=aux_index;
                        }
                    }
                }
            }
            //G.2. Calculates optimum p_prime (2nd level)
            //Now we look in the vicinity of the first election of p prime
            //If the number of data is low enough, we will use each value
            //Otherwise we just define another set of values (200 this time)
            if(grid_size<=100)
            {
                first=3;
                last=grid_size-3;
            }
            else
            {
                first=p_aux_index[aux_index_min-1]+1;
                last=p_aux_index[aux_index_min+1];
                
                if(x[first]<x_min) first=range[0];
                if(x[last]>x_max) last=range[1];
            }
            if(grid_size<=10000)
            {
                for(p_prime_index=first;p_prime_index<last;p_prime_index++)
                {
                    initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                    }
                }
            }
            else
            {
                delta_int_2=(int)((last-first)/199);
                for(i=0;i<199;i++)
                    p_aux_index_2[i]=first+i*delta_int_2;
                p_aux_index_2[199]=last-1;
                
                for(aux_index=0;aux_index<200;aux_index++)
                {
                    p_prime_index=p_aux_index_2[aux_index];
                    initial_point=greatest_convex_minorant(x,F,0,p_prime_index,F_g,f_g,grid_size,f);
                    final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
                    error_aux=grenander_error(F,F_g,x,grid_size);
                    if(error_aux<error)
                    {
                        error=error_aux;
                        p_prime_index_min=p_prime_index;
                        aux_index_min=aux_index;
                    }
                }
            }
            //G.3 Calculating best distribution (once p_prime is determined)
            //printf("p_prime %lf\n",x[p_prime_index_min]);
            initial_point=greatest_convex_minorant(x,F,0,p_prime_index_min,F_g,f_g,grid_size,f);
            final_point=least_concave_majorant(x,F,initial_point,grid_size-1,F_g,f_g,grid_size,f);
            //G.4 Correcting the slopes (Grenander's PDFs)
            //f_g[0]=f_g[1];
            //for(i=1;i<(grid_size-2);i++)
            //	f_g[i]=0.5*(f_g[i]+f_g[i+1]);
            correcting_slopes(f_g,x,grid_size);
            //G.5 Renormalizing Grenander's PDFs
            area_grenander=f_g[0]*(x[1]-x[0]);
            for(i=1;i<(grid_size-1);i++)
                area_grenander+=f_g[i]*(x[i+1]-x[i-1])*0.5;
            area_grenander+=f_g[grid_size-1]*(x[grid_size-1]-x[grid_size-2]);
            area_grenander=1.0/area_grenander;
            for(i=0;i<grid_size;i++)
                f_g[i]=f_g[i]*area_grenander;
        }
        goodness=0;
    }

	if((flag_choice[flag]==1)||(flag_choice[flag]==3))
		index_p_prime[flag]=-1;
	return goodness;
}

int small_grid_calculation(double *f, double *x, int grid_size, double *f_small, double *x_small, int N)
{
	int i,delta_int,aux_i;
	int index[N+1];
	double delta,resto,resto_acum=0;
	double x1,x2;
	double x2_small;
	double distance[N];
	
	delta=(double)grid_size/(N);
	delta_int=(int)delta;
	resto=delta-delta_int;
	for(i=0;i<N;i++)
	{
		index[i]=i*delta_int;
		resto_acum+=resto;
		
		if(resto_acum>1)
		{
			resto_acum=resto_acum-1;
			index[i]++;
		}
	}
	index[N]=grid_size-1;
		
	for(i=0;i<N;i++)
	{
		x_small[i]=0.5*(x[index[i]+1]+x[index[i]]);
		f_small[i]=0;
	}
		
	for(i=0;i<N;i++)
		distance[i]=x[index[i+1]]-x[index[i]];
	
	/*DEBUG*/
	//FILE *fp;
	//fp=fopen("small_grid.txt","wt");
	//fprintf(fp,"Aux_index Real_index x distance\n");
	//for(i=0;i<N;i++)
	//	fprintf(fp,"%d %d %lf %lf\n",i,index[i],x[index[i]],distance[i]);
	//fclose(fp);
	/* */
		
	aux_i=0;
	for(i=0;aux_i<(N-1);i++)
	{	
		if(i>0)	
			x1=x[i]-0.5*(x[i]-x[i-1]);
		else
			x1=x[0];
		
		if(i<(grid_size-1))
			x2=x[i]+0.5*(x[i+1]-x[i]);
		else
			x2=x[grid_size-1];
		
		if(aux_i<(N-1))
			x2_small=x_small[aux_i]+0.5*(x_small[aux_i+1]-x_small[aux_i]);
		else
			x2_small=x[grid_size-1];
		
		if(x2<x2_small)	
			f_small[aux_i]+=f[i]*(x2-x1);
		else
		{			
			f_small[aux_i]+=(x2_small-x1)*f[i];
			aux_i++;
			f_small[aux_i]+=(x2-x2_small)*f[i];
		}
	}
	
	for(i=0;i<N;i++)
		f_small[i]/=distance[i];
		
	/*DEBUG*/
	//fp=fopen("original.txt","wt");
	//fprintf(fp,"x p(x)\n");
	//for(i=0;i<grid_size;i++)
	//	fprintf(fp,"%lf %lf\n",x[i],f[i]);
	//fclose(fp);
	//fp=fopen("reduced.txt","wt");
	//fprintf(fp,"x p(x)\n");
	//for(i=0;i<N;i++)
	//	fprintf(fp,"%lf %lf\n",x_small[i],f_small[i]);
	//fclose(fp);
    return(1);
}

int explained_vs_unexplained_variance(double *x, double *x_small, int p_prime_index, int n_small, double *f, double *f_g)
{
	int i,new_p_prime_index,flag_decision;
	double sigma_exp[2], sigma_unexp[2];
	double R_square[2];
	double f_line;
	
	//A. Position of p_prime in the reduced grid
	for(i=0;((x[p_prime_index]>x_small[i])&&(i<n_small));)
		i++;		
	if(x[p_prime_index]>x_small[n_small-1]) new_p_prime_index=n_small-1;
	else
	{
		if(x[p_prime_index]<x_small[0]) new_p_prime_index=0;
		else
		{			
			if((x[p_prime_index]-x_small[i-1])<(x_small[i]-x[p_prime_index]))
				new_p_prime_index=i-1;
			else
				new_p_prime_index=i;
		}
	}	
	
	//B Calculates the average (with Grenander data?)
	//This is not the exact average of the grenander function, as we are not considering the spacing between points (change it?)
	/*average[0]=0;
	for(i=0;i<new_p_prime_index;i++)
		average[0]+=f_g[i];
	average[0]/=new_p_prime_index;
	average[1]=0;
	for(i=new_p_prime_index;i<n_small;i++)
		average[1]+=f_g[i];
	average[1]/=n_small-new_p_prime_index;*/
	
	//C. Calculates explained variance
	/*sigma_exp[0]=0;
	for(i=0;i<new_p_prime_index;i++)
		sigma_exp[0]+=(f_g[i]-average[0])*(f_g[i]-average[0]);
	sigma_exp[0]=sqrt(sigma_exp[0]/(new_p_prime_index-1));
	sigma_exp[1]=0;
	for(i=new_p_prime_index;i<n_small;i++)
		sigma_exp[1]+=(f_g[i]-average[1])*(f_g[i]-average[1]);
	sigma_exp[1]=sqrt(sigma_exp[1]/(n_small-new_p_prime_index-1));*/
	
	sigma_exp[0]=0;
	for(i=0;i<new_p_prime_index;i++)
	{
		f_line=0+(f[new_p_prime_index]-0)/(x_small[new_p_prime_index]-0)*(x_small[i]-0);
		sigma_exp[0]+=(f_g[i]-f_line)*(f_g[i]-f_line);
	}		
	sigma_exp[0]=sqrt(sigma_exp[0]/(new_p_prime_index-1));
	sigma_exp[1]=0;
	for(i=new_p_prime_index;i<n_small;i++)
	{
		f_line=f[new_p_prime_index]+(1-f[new_p_prime_index])/(x_small[n_small-1]-x_small[new_p_prime_index])*(x_small[i]-x_small[new_p_prime_index]);
		sigma_exp[1]+=(f_g[i]-f_line)*(f_g[i]-f_line);
	}
	sigma_exp[1]=sqrt(sigma_exp[1]/(n_small-new_p_prime_index-1));
	
	
	//D. Calculates unexplained variance
	sigma_unexp[0]=0;
	for(i=0;i<new_p_prime_index;i++)
		sigma_unexp[0]+=(f_g[i]-f[i])*(f_g[i]-f[i]);
	sigma_unexp[0]=sqrt(sigma_unexp[0]/(new_p_prime_index-1));
	sigma_unexp[1]=0;
	for(i=new_p_prime_index;i<n_small;i++)
		sigma_unexp[1]+=(f_g[i]-f[i])*(f_g[i]-f[i]);
	sigma_unexp[1]=sqrt(sigma_unexp[1]/(n_small-new_p_prime_index-1));
	
	if(global_flag==1)
	{
		printf("Explained variance: %lf %lf\n",sigma_exp[0],sigma_exp[1]);
		printf("Unexplained variance: %lf %lf\n",sigma_unexp[0],sigma_unexp[1]);
	}
	
	for(i=0;i<2;i++)
		R_square[i]=1-sigma_unexp[i]/(sigma_unexp[i]+sigma_exp[i]);
	
	if ((sigma_exp[0]<2*sigma_unexp[0])||(sigma_exp[1]<2*sigma_unexp[1]))
		flag_decision=1;
	else
		flag_decision=0;
		
	return flag_decision;
}

double goodness_calc(double *f, double *f_g, double *x, int n)
{
	double goodness=0;
	int i;
	
	goodness+=(f_g[0]-f[0])*(-1+2*((f_g[0]-f[0])>0))*(x[1]-x[0]);
	for(i=1;i<(n-1);i++)
		goodness+=(f_g[i]-f[i])*(-1+2*((f_g[i]-f[i])>0))*0.5*(x[i+1]-x[i-1]);		
	goodness+=(f_g[n-1]-f[n-1])*(-1+2*((f_g[n-1]-f[n-1])>0))*(x[n-1]-x[n-2]);
	goodness=1-0.5*goodness;
	
	return goodness;
}

extern void estimating_grid_log_aux(double *grid, 
									double *grid_aux, 
									int *grid_size,
									double *y,
									double *y_aux)
{
	int i,c;
	double delta=1e-100;
	
	for(i=0;i<grid_size[0];i++)
		grid_aux[i]=-1;
	
	c=1;
	grid_aux[0]=grid[0];
	y_aux[0]=y[0];
	for(i=1;i<grid_size[0];i++)
	{
		if((grid[i]-grid_aux[c-1])>delta)
		{
			grid_aux[c]=grid[i];
			y_aux[c]=y[i];
			c++;
		}
	}
	printf("c=%d\n",c);
}

extern void interpolate_MV(double *grid,
							double *grid_aux,
							double *pi_hat_long,
							double *pi_hat_aux,
							int *grid_size,
							int *grid_size_aux)
{
	int i,c=1;
	
	pi_hat_long[0]=pi_hat_aux[0];
		
	for(i=1;(i<grid_size[0])&&(c<grid_size_aux[0]);i++)
	{
		if(grid[i]<=grid_aux[c])
			pi_hat_long[i]=pi_hat_aux[c-1]+(pi_hat_aux[c]-pi_hat_aux[c-1])/(grid_aux[c]-grid_aux[c-1])*(grid[i]-grid_aux[c-1]);
		else
		{
			c++;
			i--;
		}
	}
}


void correcting_slopes(double *f_g,double *x,int size)
{
	int i;
	
	for(i=1;i<(size-1);i++)
		f_g[i]=(f_g[i]*(x[i]-x[i-1])+f_g[i+1]*(x[i+1]-x[i]))/(x[i+1]-x[i-1]);
	
}

void calculate_range(double x_min, double x_max, double *x, int *range, int n)
{
	int i;
	for(i=0;x[i]<x_min;i++)
	range[0]=i;
	for(i=(n-1);x[i]>x_max;i--)
	range[1]=i;	
}
