#include<stdio.h>
#include<stdlib.h>
#define __USE_C99_MATH
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<mpi.h>
#include<sys/time.h>

long int product(int *array, int n) {
    long int product = 1;
    for(int i=0; i<n; i++) {
        product *= array[i];
    }
    return product;
}

int *read_dims(char *filename) {
    FILE *file = fopen(filename,"r");
    
    if(file == NULL) {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

    char firstline[500];
    fgets(firstline, 500, file);
    
    int line_length = strlen(firstline);

    int num_dims = 0;
    for(int i=0; i<line_length; i++) {
        if(firstline[i] == ' ') {
            num_dims++;
        }
    }
    
    int *dims = malloc((num_dims+1)*sizeof(int));
    dims[0] = num_dims;
    const char s[2] = " ";
    char *token;
    token = strtok(firstline, s);
    int i = 0;
    while( token != NULL ) {
        dims[i+1] = atoi(token);
        i++;
        token = strtok(NULL, s);
    }
    fclose(file);
    return dims;
}

double * read_array(char *filename, int *dims, int num_dims) {
    FILE *file = fopen(filename,"r");

    if(file == NULL) {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

    char firstline[500];
    fgets(firstline, 500, file);

    //Ignore first line and move on since first line contains 
    //header information and we already have that. 

    long int total_elements = product(dims, num_dims);

    double *one_d = malloc(sizeof(double) * total_elements);

    for(int i=0; i<total_elements; i++) {
        fscanf(file, "%lf", &one_d[i]);
    }
    fclose(file);
    return one_d;
}

int write_array(char *filename, int total_elements, double *output){
    FILE *file = fopen(filename,"w");

    if(file == NULL) {
        printf("Unable to open file: %s", filename);
        return -1;
    }

    if (file != NULL) {
        fprintf(file, "%d ", total_elements);
        fprintf(file, "\n");
    }

    for(int i=0; i<total_elements; i++) {
        fprintf(file, "%.7f ", output[i]);
    }

    fclose(file);
    return 1;
}

void bubblingSort(double *arr, int n){
	int i, j;
	printf("bubbling!\n");
	for (i = 0; i < n; ++i){
		for (j = 0; j < n - i - 1; ++j){
			if (arr[j] > arr[j + 1]) {
				double temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
			}
		}
	}
}


void Merge(double *sourceArr,double *tempArr, int startIndex, int midIndex, int endIndex){
    int i = startIndex, j=midIndex+1, k = startIndex;
    while(i!=midIndex+1 && j!=endIndex+1) {
        if(sourceArr[i] > sourceArr[j])
            tempArr[k++] = sourceArr[j++];
        else
            tempArr[k++] = sourceArr[i++];
    }
    while(i != midIndex+1)
        tempArr[k++] = sourceArr[i++];
    while(j != endIndex+1)
        tempArr[k++] = sourceArr[j++];
    for(i=startIndex; i<=endIndex; i++)
        sourceArr[i] = tempArr[i];
}
 
void traditional_MergeSort(double *sourceArr, double *tempArr, int startIndex, int endIndex) {
    int midIndex;
    if(startIndex < endIndex) {
        midIndex = startIndex + (endIndex-startIndex) / 2;
        traditional_MergeSort(sourceArr, tempArr, startIndex, midIndex);
        traditional_MergeSort(sourceArr, tempArr, midIndex+1, endIndex);
        Merge(sourceArr, tempArr, startIndex, midIndex, endIndex);
    }
}


void MergeSort(double *arr, int start, int end, int s, int all_size){	
	//printf("merging!\n");
	double *temp = malloc(sizeof(double) * all_size);;
	if (start >= end)		
		return;	
	int mid = end - s;	
 
	int length = 0; 
	int i_start = start;
	int i_end = mid;
	int j_start = mid + 1;
	int j_end = end;
	while (i_start <= i_end && j_start <= j_end){		
		if (arr[i_start] < arr[j_start]){			
			temp[length] = arr[i_start];
			length++;			
			i_start++;		
		}		
		else{			
			temp[length] = arr[j_start];
			length++;		
			j_start++;		
		}	
	}	
	while (i_start <= i_end){	
		temp[length] = arr[i_start];
		i_start++;		
		length++;	
	}	
	while (j_start <= j_end){		
		temp[length] = arr[j_start];
		j_start++;
		length++;			
	}	
	for (int i = 0; i < length; i++){		
		arr[start + i] = temp[i];
	}
}

int  main(int argc, char* argv[]) {
    //Define the number and id of threads
	int process_id, process_num;
    //Packets for broadcasting
    int matrix_properties;
	int namelen;
    double MPI_timer[4];
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    //Raw data read from thread 0, belongs to address header, no space allocated
    double *input_data = NULL;
    bool match = true;
    char input_filename[500];
    char check_filename[500];
	MPI_Init(&argc, &argv);
    MPI_timer[0] = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	MPI_Get_processor_name(processor_name, &namelen);
	if(process_id == 0) {
		if(argc != 3) {
            printf("Usage: %s <filename_input> <filename_output>\n", argv[0]);
            return -1;
        }
        int compareOutput = 1;
        strcpy(input_filename, argv[1]);
        strcpy(check_filename, argv[2]);
        int *input_dims_original = read_dims(input_filename);
        if(input_dims_original == NULL) {
            return -1;
        }
        int input_num_dims = input_dims_original[0];
        int *input_dims = input_dims_original+1;
        input_data = read_array(input_filename, input_dims, input_num_dims);
        if(input_data == NULL) {
            return -1;
        }
        //Array length
        int array_size = input_dims[0];
        matrix_properties = array_size;
        printf("have threads %d\n", process_num);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_timer[1] = MPI_Wtime();
	MPI_Bcast(&matrix_properties, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int chunk_size = matrix_properties / process_num;
    //Create memory space for input and chunck storage for each process
	double *input = malloc(sizeof(double) * matrix_properties);
    //Space for each thread to store data
	double *blocks = malloc(sizeof(double) * chunk_size);
    //Create space for data to be placed after each thread's calculation
    double *back_blocks = malloc(sizeof(double) * chunk_size);
	if(process_id == 0){
        input = input_data;
    }
    //Separate inputs into blocks and send them to each thread
	MPI_Scatter(input, chunk_size, MPI_DOUBLE, blocks, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//printf("The process %d of %d is on %s\n", process_id, process_num, processor_name);
	//bubblingSort(blocks, chunk_size);
    //Calculation
    traditional_MergeSort(blocks, back_blocks, 0, chunk_size - 1);
    //gather the results of each thread's operations
	MPI_Gather(back_blocks, chunk_size, MPI_DOUBLE, input, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    //Final redistribution of data at process 0
	if(process_id == 0){
		for (int p = 2; p <= process_num ; p++) 
		{
			MergeSort(input, 0, (p * chunk_size) - 1, chunk_size, matrix_properties);
		}
	}
    MPI_timer[2] = MPI_Wtime();
    printf("%f\n", MPI_timer[2] - MPI_timer[1]);
    //Process 0 does the file writing
	if(process_id == 0){
		int write = write_array(check_filename, matrix_properties, input);
        if(write == 1) {
            printf("Writing successful!\n");
        }
	}
	MPI_Finalize();
	return 0;
}
