/* main code to compute a power spectral density.
using without any arguments or -h gives usage information.
*/

#include "station_analysis_tools.h"


int main (int argc, char *argv[]) {
	/* local variable declaration */
	int i, ii, l, j, acc_flag=1;
	int start_sample, samples_per_second, temp_npts;
	int n_psds;
	float n_second_sub_record, n_sec_shift;
	FILE *paramfile;
	FILE *output_file;
	char paramfile_name[100];
	char buff[100];
	char output_file_name[100];
	char sac_file_name[100];
	char pz_file_name[100];
	float *temp_array_f;
	double *signal;
	double *temp_signal, *temp_psd;
	double *mean_psd, *sd_psd;
	SAC_HD sac_header;
	fftw_complex *temp_spect;
	POLE_ZERO pz;
	errno_t err;
	/* temporary arrays now used in signal_to_ground_motion to "simplify" memory management */
	/*
	 * double *signal, double *untapered, SPECTRUM *signal_spectrum, SPECTRUM *response_spectrum,
	 * SPECTRUM *untapered_spectrum, double *fftw_dbl, fftw_complex *fftw_cmplx,
	 * fftw_plan plan_forward, fftw_plan plan_backward
	 */
	double *stga_signal, *stga_fftw_dbl;
	fftw_complex *stga_signal_spectrum, *stga_response_spectrum;
	fftw_complex *stga_fftw_cmplx;
	fftw_plan stga_plan_forward, stga_plan_backward;

	//check for -h usage. Should help prevent misuse.
	if (argc==2) {
		if (!(strcmp(argv[1],"-h"))) {
			fprintf(stderr,"USAGE: %s [parameter_file]\n", argv[0]);
			fprintf(stderr,"parameter file format:\n");
			fprintf(stderr,"sac_file\n");
			fprintf(stderr,"polezero_file\n");
			fprintf(stderr,"number_of_seconds_in_each_record number_of_seconds_offset_between_records acc_flag\n");
			fprintf(stderr,"output_file\n\n\n");
			fprintf(stderr,"The full .SAC file will be divided into N records chosen by number_of_seconds_in_each_record and number_of_seconds_offset_between_records.\n");
			fprintf(stderr,"As an example, putting 28800 and 3600 for those parameters for a 24 hour 1 Hz SAC file computes 17 psds estimates each with length of 8 hours overlapping 1 hour\n");
			fprintf(stderr,"acc_flag can be 0, 1, or 2. This determines the number of derivatives computed - 0 = displacement, 1 = velocity, 2 = acceleration\n");
			exit(1);
		}
	}

	/* open and read the parameter file */
	/* grab command line arguments */
	if(argc!=2) {
		fprintf(stderr,"USAGE: %s [parameter_file]\n", argv[0]);
		fprintf(stderr,"parameter file format:\n");
		fprintf(stderr,"sac_file\n");
		fprintf(stderr,"polezero_file\n");
		fprintf(stderr,"number_of_seconds_in_each_record number_of_seconds_offset_between_records acc_flag\n");
		fprintf(stderr,"output_file\n\n\n");
		fprintf(stderr,"The full .SAC file will be divided into N records chosen by number_of_seconds_in_each_record and number_of_seconds_offset_between_records.\n");
		fprintf(stderr,"As an example, putting 28800 and 3600 for those parameters for a 24 hour 1 Hz SAC file computes 17 psds estimates each with length of 8 hours overlapping 1 hour\n");
		fprintf(stderr,"acc_flag can be 0, 1, or 2. This determines the number of derivatives computed - 0 = displacement, 1 = velocity, 2 = acceleration\n");
		exit(1);
	}

	// copy the argument to the parameter filename
//	strcpy(paramfile_name, argv[1]);
	strcpy_s(paramfile_name,100, argv[1]);

	//open the parameter file
	
	if ((err = fopen_s(&paramfile,paramfile_name,"r"))) {
		fprintf(stderr,"Error: file %s not found\n", paramfile_name);
		exit(1);
	}
        
	// read the parameter file and store its contents
	l=0;
	while (fgets(buff,100,paramfile)) {
		if (l == 0) {
			sscanf_s(buff,"%s", sac_file_name, (unsigned)_countof(sac_file_name));
			l++;
		} else if (l == 1) {
			sscanf_s(buff,"%s",pz_file_name, (unsigned)_countof(pz_file_name));
			l++;
		} else if (l == 2) {
			sscanf_s(buff,"%f %f %d",&n_second_sub_record, &n_sec_shift, &acc_flag);
			l++;
     	 	} else if (l == 3) {
			sscanf_s(buff,"%s", output_file_name, (unsigned)_countof(pz_file_name));
			l++;
		} 
	}

	fclose(paramfile);

	//run some checks on the parameters
	//first check is for a reasonable unit to compute
	if (acc_flag != 0 && acc_flag != 1 && acc_flag != 2) {
		fprintf(stderr,"Error, acc_flag must equal 0, 1, or 2. Determines the number of differentiations from the displacement response.\n");
		exit(1);
	}
	
	//second check that the sub record and time shift are non-zero
	if (n_second_sub_record <= 0.0 ) {
		fprintf(stderr,"Error, n_second_sub_record must be greater than 0.\n");
		exit(1);
	}
	if (n_sec_shift <= 0.0 ) {
		fprintf(stderr,"Error, n_sec_shift must be greater than 0.\n");
		exit(1);
	}

	//read the sac files
	if ( ( temp_array_f = malloc(sizeof(*temp_array_f) * NPTS_SAC_FILE_MAX) ) == NULL) {
		fprintf(stderr,"Error allocating memory. Reduce the value of NPTS_SAC_FILE_MAX in station_analysis_tools.h and recompile.\n");
		exit(1);
	}
	if ( ( signal = malloc(sizeof(*signal) * NPTS_SAC_FILE_MAX) ) == NULL) {
		fprintf(stderr,"Error allocating memory. Reduce the value of NPTS_SAC_FILE_MAX in station_analysis_tools.h and recompile.\n");
		exit(1);
	}

	if (!(read_sac ( sac_file_name, temp_array_f, &(sac_header), NPTS_SAC_FILE_MAX))) {
		fprintf(stderr,"file %s not found\n", sac_file_name );
		exit(1);
	}
	for (i=0; i<sac_header.npts; i++) {
		signal[i] = (double) temp_array_f[i];
	}
	
	free(temp_array_f);
	
	//read the pole zero files
	if ( read_sac_pole_zero_file(pz_file_name, &pz) == -1) {
		fprintf(stderr,"Error, file %s not found\n",pz_file_name);
		exit(1);
	}

	//invert delta to samples per second
	samples_per_second = floor(1.0/sac_header.delta + 0.5);
	
	//hard wired for 8 hour data segments. (8*3600 = seconds in 8 hours)
	temp_npts = (int) floor(n_second_sub_record * samples_per_second+0.5);
	
	//allocate memory for working arrays used in the next few steps.
	if ((temp_signal = malloc(sizeof(*temp_signal) * temp_npts)) == NULL) {
		fprintf(stderr,"Error allocating temp_signal\n");
		exit(1);
	}
	if ((temp_psd = malloc(sizeof(*temp_psd) * (temp_npts / 2+1))) == NULL) {
		fprintf(stderr,"Error allocating temp_psd\n");
		exit(1);
	}

	if ((mean_psd = malloc(sizeof(*mean_psd) * (temp_npts / 2+1))) == NULL) {
		fprintf(stderr,"Error allocating mean_psd\n");
		exit(1);
	}
	if ((sd_psd = malloc(sizeof(*sd_psd) * (temp_npts / 2+1))) == NULL) {
		fprintf(stderr,"Error allocating sd_psd\n");
		exit(1);
	}

	if ((temp_spect = fftw_malloc(sizeof(*temp_spect) * (temp_npts/2+1))) == NULL) {
		fprintf(stderr,"Error allocating temp_spect\n");
		exit(1);
	}
	
	if ((stga_signal = fftw_malloc(sizeof(*stga_signal) * temp_npts)) == NULL) {
		fprintf(stderr,"Error allocating stga_signal\n");
		exit(1);
	}
	if ((stga_fftw_dbl = fftw_malloc(sizeof(*stga_fftw_dbl) * temp_npts)) == NULL) {
		fprintf(stderr,"Error allocating stga_fftw_dbl\n");
		exit(1);
	}
	if ((stga_fftw_cmplx = fftw_malloc(sizeof(*stga_fftw_cmplx) * (temp_npts/2+1))) == NULL) {
		fprintf(stderr,"Error allocating stga_fftw_cmplx\n");
		exit(1);
	}
	
	if ((stga_signal_spectrum = fftw_malloc(sizeof(*stga_signal_spectrum) * (temp_npts/2+1))) == NULL) {
		fprintf(stderr,"Error allocating stga_signal_spectrum\n");
		exit(1);
	}
	if ((stga_response_spectrum = fftw_malloc(sizeof(*stga_response_spectrum) * (temp_npts/2+1))) == NULL) {
		fprintf(stderr,"Error allocating stga_response_spectrum\n");
		exit(1);
	}
		
	
	stga_plan_forward = fftw_plan_dft_r2c_1d(temp_npts, stga_fftw_dbl, stga_fftw_cmplx, FFTW_ESTIMATE);
	stga_plan_backward = fftw_plan_dft_c2r_1d(temp_npts, stga_fftw_cmplx, stga_fftw_dbl, FFTW_ESTIMATE);
	
	n_psds = ( sac_header.npts - temp_npts ) / (int) (n_sec_shift * (float) samples_per_second) + 1;

	//clear out the mean and sd...
	for (ii=0; ii<temp_npts/2+1; ii++) {
		mean_psd[ii] = 0.0;
		sd_psd[ii] = 0.0;
	}
	
	for (j = 0; j < n_psds; j++ ) {
		//hard wired to step 1 hour for each data segment
		start_sample = j * n_sec_shift * samples_per_second;
		//make sure we don't go out of bounds
		if (start_sample + temp_npts < sac_header.npts) {
			//fill with data from the current segment
			for (ii=0; ii<temp_npts; ii++) { 
				temp_signal[ii] = signal[start_sample+ii];
			}
			//convert to acceleration in m/s/s
			signal_to_ground_acceleration(temp_signal, temp_signal, temp_npts, sac_header.delta, &(pz), acc_flag,
										  stga_signal, stga_signal_spectrum, stga_response_spectrum,
										  stga_fftw_dbl, stga_fftw_cmplx,
										  stga_plan_forward, stga_plan_backward);
						
			//fft
			for (ii=0; ii<temp_npts; ii++) {
				stga_fftw_dbl[ii] = temp_signal[ii];
			}
			fftw_execute(stga_plan_forward);
			for (ii=0; ii<temp_npts/2+1; ii++) {
				temp_spect[ii][0] = stga_fftw_cmplx[ii][0];
				temp_spect[ii][1] = stga_fftw_cmplx[ii][1];
			}
			//compute the power 
			calculate_psd(temp_spect, temp_psd, sac_header.delta, temp_npts);
			//for each period, compute the mean and standard deviation recursively
			if (j == 0) {
				for (ii=0; ii<temp_npts/2+1; ii++) {
					mean_psd[ii] = temp_psd[ii];
					sd_psd[ii] = 0.0;
				}
//				printf("%lf +/- %lf\n",mean_psd[20], sd_psd[20]);
			} else {	
				for (ii=0; ii<temp_npts/2+1; ii++) {
					mean_psd[ii] = recursive_mean(temp_psd[ii], mean_psd[ii], j+1);
					sd_psd[ii] = recursive_standard_deviation(temp_psd[ii], mean_psd[ii], sd_psd[ii], j+1);
				}
//				printf("%lf +/- %lf\n",mean_psd[20], sd_psd[20]);
			}
		}
		
	} // done with each segment

	//simple message to tell the user how many psds were calculated
	fprintf(stdout,"%d psds computed and averaged\n",n_psds);

	free(signal);
	free(temp_signal);
	fftw_free(temp_spect);
	fftw_free(stga_signal_spectrum);
	fftw_free(stga_response_spectrum);
	fftw_free(stga_signal);
	fftw_free(stga_fftw_dbl);
	fftw_free(stga_fftw_cmplx);
	fftw_destroy_plan(stga_plan_forward);
	fftw_destroy_plan(stga_plan_backward);
	
	
	//output
	if ((err=  fopen_s(&output_file,output_file_name,"w"))) {
		fprintf(stderr, "Error opening %s to write\n",output_file_name);
		exit(1);
	}
	for (i=1; i<temp_npts/2+1; i++) {
		//printf("%3.5lf %3.5lf %3.5lf\n" , 1.0 / ( 1.0 / sac_header.delta * (i) / temp_npts), mean_psd[i], sd_psd[i]);
		fprintf(output_file,"%3.5lf %3.5lf %3.5lf\n" , 1.0 / ( 1.0 / sac_header.delta * (i) / temp_npts), mean_psd[i], sd_psd[i]);
	}
	fclose(output_file);

	
	//free up the memory like a responsible citizen - the OS should do this anyways...

	free(mean_psd);
	free(sd_psd);
		
	return 0;
}
