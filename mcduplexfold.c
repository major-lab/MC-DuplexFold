/****************************************************************************
 * MC-DuplexFold (mcdf) v 1.1
 * Copyright 2025 Simon Chasles <simon.chasles@umontreal.ca>
 * 
 * For basic compilation:   gcc -o mcdf -fopenmp mcduplexfold.c -lm
 * For basic execution:     ./mcdf -s {first RNA sequence} -t {second RNA sequence}
 * Sirt-1:miR34a example:   ./mcdf -s ACACCCAGCUAGGACCAUUACUGCCA -t UGGCAGUGUCUUAGCUGGUUGU -p 1
 * Sirt-1:miR34a example:   ./mcdf -s ACACCCAGCUAGGACCAUUACUGCCA -t UGGCAGUGUCUUAGCUGGUUGU -f 5
 * HNF4a:mir34a example:    ./mcdf -s aacauggccuaagggccacaucccacugcca -t UGGCAGUGUCUUAGCUGGUUGU -p 1
 * For help:                ./mcdf -h
 * Note: mcff must be executable in the current directory.
 ****************************************************************************/

#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <stdbool.h>
#include <inttypes.h>
#include "MRG32k3aOMP.h"
#include "leftE.h"
#include "rightE.h"
#include "stackE.h"

#define ncm_nb 13
#define two 2
#define default_gap_id 0
#define max_seq_size 101
#define mcff_structure_nb 500
#define max_mcff_output_length 250
#define INF_3 1000
#define INF_4 10000
#define INF_5 20000
#define possible_prime_nb 186
#define max_possible_prime 51991

char *mcff_loop = "GAAA";
char *mcff_1 = "mcff -s ";
char *mcff_2 = " -m \'";
char *mcff_3 = "\' -t 2.0";
int ncm_gaps[ncm_nb][two] = {{1, 1},
							 {1, 2}, {2, 1}, {2, 2},
							 {1, 3}, {3, 1}, {2, 3}, {3, 2}, {3, 3},
							 {1, 4}, {4, 1}, {2, 4}, {4, 2}};
int ncm_table[5][5] = {{-1, -1, -1, -1, -1},
					   {-1,  0,  1,  4,  9},
					   {-1,  2,  3,  6, 11},
					   {-1,  5,  7,  8, -1},
					   {-1, 10, 12, -1, -1}};
int n = 0;
int m = 0;
int max_job_nb = 1;
long int warm_start = 1;
long int print_stats = 0;
long int simulation_nb = 100;
long int suboptimal_nb = 0;
long int K = 10;
double alpha = 0.45;
long int B = 5;
double start_RT = 3.0;
double end_RT = 0.333333333;
long int MRG32k3a_seed = 1;
int first_primes[max_seq_size] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547};
int possible_primes[possible_prime_nb] = {50021, 50023, 50033, 50047, 50051, 50053, 50069, 50077, 50087, 50093, 50101, 50111, 50119, 50123, 50129, 50131, 50147, 50153, 50159, 50177, 50207, 50221, 50227, 50231, 50261, 50263, 50273, 50287, 50291, 50311, 50321, 50329, 50333, 50341, 50359, 50363, 50377, 50383, 50387, 50411, 50417, 50423, 50441, 50459, 50461, 50497, 50503, 50513, 50527, 50539, 50543, 50549, 50551, 50581, 50587, 50591, 50593, 50599, 50627, 50647, 50651, 50671, 50683, 50707, 50723, 50741, 50753, 50767, 50773, 50777, 50789, 50821, 50833, 50839, 50849, 50857, 50867, 50873, 50891, 50893, 50909, 50923, 50929, 50951, 50957, 50969, 50971, 50989, 50993, 51001, 51031, 51043, 51047, 51059, 51061, 51071, 51109, 51131, 51133, 51137, 51151, 51157, 51169, 51193, 51197, 51199, 51203, 51217, 51229, 51239, 51241, 51257, 51263, 51283, 51287, 51307, 51329, 51341, 51343, 51347, 51349, 51361, 51383, 51407, 51413, 51419, 51421, 51427, 51431, 51437, 51439, 51449, 51461, 51473, 51479, 51481, 51487, 51503, 51511, 51517, 51521, 51539, 51551, 51563, 51577, 51581, 51593, 51599, 51607, 51613, 51631, 51637, 51647, 51659, 51673, 51679, 51683, 51691, 51713, 51719, 51721, 51749, 51767, 51769, 51787, 51797, 51803, 51817, 51827, 51829, 51839, 51853, 51859, 51869, 51871, 51893, 51899, 51907, 51913, 51929, 51941, 51949, 51971, 51973, 51977, 51991};
int chosen_prime = 0;
int hash_inverse[max_possible_prime];
int entropy_structure_nb = 500;

int cannot_cast_str_to_int(char *string, char param) {
	char *endptr;
	strtol(string, &endptr, 10);
	if (endptr == string) {
		printf("No digits were found after option -%c\n", param);
		return 1;
	} else if (*endptr != '\0') {
		printf("Invalid character %c after option -%c\n", *endptr, param);
		return 1;
	} else {
		return 0;
	}
}

int cannot_cast_str_to_double(char *string, char param) {
	char *endptr;
	strtod(string, &endptr);
	if (endptr == string) {
		printf("No digits were found after option -%c\n", param);
		return 1;
	} else if (*endptr != '\0') {
		printf("Invalid character %c after option -%c\n", *endptr, param);
		return 1;
	} else {
		return 0;
	}
}

int minimum(int n, int m) {
	if (n < m) {return n;}
	return m;
}

int maximum(int n, int m) {
	if (n > m) {return n;}
	return m;
}

_Bool is_ncm(int p, int q) {
	if (p > 4 || q > 4) {
		return 0;
	}
	return ncm_table[p][q] > -1;
}

int get_seq_index(int u, int i, int j, int v) {
	return 64 * u + 16 * i + 4 * j + v;
}

double sigmoid_of_negative(double x) {
	return 1.0 / (1.0 + exp(x));
}

int encode_sequence(int length, char seq[length + 1], int encoded_seq[length]) {
	for (int i = 0; i < length; i++) {
		switch (seq[i]) {
			case 'a':
				encoded_seq[i] = 0;
				break;
			case 'A':
				encoded_seq[i] = 0;
				break;
			case 'c':
				encoded_seq[i] = 1;
				break;
			case 'C':
				encoded_seq[i] = 1;
				break;
			case 'g':
				encoded_seq[i] = 2;
				break;
			case 'G':
				encoded_seq[i] = 2;
				break;
			case 'u':
				encoded_seq[i] = 3;
				break;
			case 'U':
				encoded_seq[i] = 3;
				break;
			case 't':
				encoded_seq[i] = 3;
				break;
			case 'T':
				encoded_seq[i] = 3;
				break;
			default:
				printf("RNA sequences should have a maximum of 100 nucleotides and only contain A, C, G, U or T.\n");
				return 1;
		}
	}
	return 0;
}

int decode_sequence(int length, char seq[length + 1], int encoded_seq[length]) {
	for (int i = 0; i < length; i++) {
		switch (encoded_seq[i]) {
			case 0:
				seq[i] = 'A';
				break;
			case 1:
				seq[i] = 'C';
				break;
			case 2:
				seq[i] = 'G';
				break;
			case 3:
				seq[i] = 'U';
				break;
			default:
				seq[i] = 'N';
		}
	}
	seq[length] = '\0';
	return 0;
}

char decode_nucleotide(int i) {
	switch (i) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'U';
		default:
			return 'N';
	}
}

int get_length_of_mcff_command(int n, int m, char seq_1[n+1], char seq_2[m+1]) {
	int command_length = 2 * (strlen(seq_1) + strlen(mcff_loop) + strlen(seq_2)) + strlen(mcff_1) + strlen(mcff_2) + strlen(mcff_3);
	return command_length;
}

void make_mcff_command(int n, int m, int l, char seq_1[n+1], char seq_2[m+1], char command[l+1]) {
	int j = 0;
	for (int i = 0; i < strlen(mcff_1); i++) {
		command[j] = mcff_1[i];
		j++;
	}
	for (int i = 0; i < strlen(seq_1); i++) {
		command[j] = seq_1[i];
		j++;
	}
	for (int i = 0; i < strlen(mcff_loop); i++) {
		command[j] = mcff_loop[i];
		j++;
	}
	for (int i = 0; i < strlen(seq_2); i++) {
		command[j] = seq_2[i];
		j++;
	}
	for (int i = 0; i < strlen(mcff_2); i++) {
		command[j] = mcff_2[i];
		j++;
	}
	for (int i = 0; i < strlen(seq_1); i++) {
		command[j] = 'p';
		j++;
	}
	for (int i = 0; i < strlen(mcff_loop); i++) {
		command[j] = '.';
		j++;
	}
	for (int i = 0; i < strlen(seq_2); i++) {
		command[j] = 'q';
		j++;
	}
	for (int i = 0; i < strlen(mcff_3); i++) {
		command[j] = mcff_3[i];
		j++;
	}
	command[j] = '\0';
}

int mcff_output_to_structure(int n, int S[n], char *mcff_output) {
	int bp_nb = 0;
	int i = n - 1;
	int j = n + 4;
	char s;
	char t;
	while (i > -1) {
		s = mcff_output[i];
		t = mcff_output[j];
		if (s == '(') {
			if (t == ')') {
				S[n - 1 - i] = j - n - 4;
				i = i - 1;
				bp_nb++;
			}
			j = j + 1;
		} else {
			S[n - 1 - i] = INF_5;
			i = i - 1;
		}
	}
	return bp_nb;
}

void mcff_line_to_structure(int n, int m, int k, char *mcff_output, unsigned short int H_str[mcff_structure_nb][n+1], double deltaGs[mcff_structure_nb]) {
	int i = n - 1;
	int j = n + 4;
	char s;
	char t;
	while (i > -1) {
		s = mcff_output[i];
		t = mcff_output[j];
		if (s == '(') {
			if (t == ')') {
				H_str[k][n - i] = j - n - 4;
				i = i - 1;
			}
			j = j + 1;
		} else {
			H_str[k][n - i] = INF_5;
			i = i - 1;
		}
	}
	int energyStart = n + m + 5;
	int energyEnd = strlen(mcff_output) - 4;
	int energyStringLength = energyEnd - energyStart;
	char energyString[energyStringLength + 1];
	int iter = 0;
	for (int position = energyStart; position < energyEnd; position++) {
		energyString[iter] = mcff_output[position];
		iter++;
	}
	double energy = atof(energyString);
	deltaGs[k] = energy;
}

void update_mcff_3(int n, int m, int try) {
	if (try == 1) {
		if (n + m <= 15) {mcff_3 = "\' -ft 600";}
		else if (n + m <= 50) {mcff_3 = "\' -ft 900";}
		else if (n + m <= 100) {mcff_3 = "\' -ft 1200";}
		else {mcff_3 = "\' -ft 1500";}
	} else if (try == 2) {
		if (n + m <= 15) {mcff_3 = "\' -ft 800";}
		else if (n + m <= 50) {mcff_3 = "\' -ft 1200";}
		else if (n + m <= 100) {mcff_3 = "\' -ft 1600";}
		else {mcff_3 = "\' -ft 2000";}
	} else if (try == 3) {
		if (n + m <= 15) {mcff_3 = "\' -ft 800";}
		else if (n + m <= 50) {mcff_3 = "\' -ft 2000";}
		else if (n + m <= 100) {mcff_3 = "\' -ft 4000";}
		else {mcff_3 = "\' -ft 8000";}
	} else if (try == 4) {
		if (n + m <= 15) {mcff_3 = "\' -ft 800";}
		else if (n + m <= 50) {mcff_3 = "\' -ft 5000";}
		else if (n + m <= 100) {mcff_3 = "\' -ft 10000";}
		else {mcff_3 = "\' -ft 20000";}
	} else {
		if (n + m <= 15) {mcff_3 = "\' -ft 800";}
		else if (n + m <= 50) {mcff_3 = "\' -ft 25000";}
		else if (n + m <= 100) {mcff_3 = "\' -ft 50000";}
		else {mcff_3 = "\' -ft 100000";}
	}
}

int get_max_repeat(int try) {
	if (try == 1) {return 100;}
	else if (try == 2) {return 250;}
	else if (try == 3) {return 500;}
	else if (try == 4) {return 1000;}
	else {return 2000;}
}

int initialize_S(int n, int m, char seq_1[n + 1], char seq_2[m + 1], int S[n], unsigned short int H_str[mcff_structure_nb][n+1], double deltaGs[mcff_structure_nb], int try) {
	FILE *mcff_stream;
	char mcff_output[mcff_structure_nb][max_mcff_output_length];
	update_mcff_3(n, m, try);
	int mcff_command_length = get_length_of_mcff_command(n, m, seq_1, seq_2);
	char command[mcff_command_length + 1];
	make_mcff_command(n, m, mcff_command_length, seq_1, seq_2, command);
	mcff_stream = popen(command, "r");
	_Bool repeated_structure, same_structure, is_empty;
	int repeated_nb;
	int max_repeat = get_max_repeat(try);
	for (int k = 0; k < mcff_structure_nb; k++) {
		repeated_nb = 0;
		do {
			fgets(mcff_output[k], max_mcff_output_length, mcff_stream);
			is_empty = mcff_output[k][0] == '\0';
			if (is_empty) {
				if (try <= 4) {return initialize_S(n, m, seq_1, seq_2, S, H_str, deltaGs, try + 1);
				} else {return -1;}
			}
			repeated_structure = 0;
			for (int l = 0; l < k; l++) {
				same_structure = 1;
				for (int i = 0; i < n + m + 4; i++) {
					if (mcff_output[l][i] != mcff_output[k][i]) {
						same_structure = 0;
						break;
					}
				}
				if (same_structure) {
					if (repeated_nb >= max_repeat) {
						if (try <= 4) {return initialize_S(n, m, seq_1, seq_2, S, H_str, deltaGs, try + 1);
						} else {return -1;}
					} else {
						repeated_nb++;
						repeated_structure = 1;
						break;
					}
				}
			}
		}
		while (repeated_structure);
	}
	pclose(mcff_stream);
	int bp_nb = mcff_output_to_structure(n, S, mcff_output[0]);
	for (int k = 0; k < mcff_structure_nb; k++) {
		mcff_line_to_structure(n, m, k, mcff_output[k], H_str, deltaGs);
	}
	return bp_nb;
}

int quick_initialize_S(int n, int m, char seq_1[n + 1], char seq_2[m + 1], int S[n]) {
	FILE *mcff_stream;
	char mcff_output[max_mcff_output_length];
	int mcff_command_length = get_length_of_mcff_command(n, m, seq_1, seq_2);
	char command[mcff_command_length + 1];
	make_mcff_command(n, m, mcff_command_length, seq_1, seq_2, command);
	mcff_stream = popen(command, "r");
	fgets(mcff_output, max_mcff_output_length, mcff_stream);
	pclose(mcff_stream);
	return mcff_output_to_structure(n, S, mcff_output);
}

int initialize_S_from_dotb(int n, int m, char seq_1[n + 1], char seq_2[m + 1], char init_dotb[n + m + 1], int S[n]) {
	int bp_nb = 0;
	int i = n - 1;
	int j = n;
	char s;
	char t;
	while (i > -1) {
		s = init_dotb[i];
		t = init_dotb[j];
		if (s == '(') {
			if (t == ')') {
				S[n - 1 - i] = j - n;
				i = i - 1;
				bp_nb++;
			}
			j = j + 1;
		} else {
			S[n - 1 - i] = INF_5;
			i = i - 1;
		}
	}
	return bp_nb;
}

void print_duplex(int n, int m, int S[n], char *seq_1, char *seq_2) {
	char pp_s[max_mcff_output_length], pp_t[max_mcff_output_length], pp_S[max_mcff_output_length];
	pp_s[0] = '\0';
	pp_t[0] = '\0';
	pp_S[0] = '\0';
	int last_i = -1, last_j = -1;
	int step_i, step_j;
	int limit_s = 0, limit_t = 0, limit_S = 0;
	char symbol = '.';
	for (int i = 0; i < n; i++) {
		if (S[i] < INF_5) {
			step_i = i - last_i;
			step_j = S[i] - last_j;
			limit_s = maximum(0, step_j - step_i);
			limit_t = maximum(0, step_i - step_j);
			limit_S = maximum(step_i, step_j) - 1;
			for (int k = 0; k < limit_s; k++) {
				symbol = '-';
				strncat(pp_s, &symbol, 1);
			}
			for (int k = last_i+1; k < i+1; k++) {
				symbol = seq_1[n-1-k];
				strncat(pp_s, &symbol, 1);
			}
			for (int k = 0; k < limit_t; k++) {
				symbol = '-';
				strncat(pp_t, &symbol, 1);
			}
			for (int k = last_j+1; k < S[i]+1; k++) {
				symbol = seq_2[k];
				strncat(pp_t, &symbol, 1);
			}
			for (int k = 0; k < limit_S; k++) {
				symbol = '.';
				strncat(pp_S, &symbol, 1);
			}
			symbol = '|';
			strncat(pp_S, &symbol, 1);
			last_i = i;
			last_j = S[i];
		}
	}
	for (int k = last_i+1; k < n; k++) {
		symbol = seq_1[n-1-k];
		strncat(pp_s, &symbol, 1);
	}
	strcat(pp_s, " <- 5\'");
	for (int k = last_j+1; k < m; k++) {
		symbol = seq_2[k];
		strncat(pp_t, &symbol, 1);
	}
	strcat(pp_t, " -> 3\'");
	limit_S = minimum(n - last_i, m - last_j) - 1;
	for (int k = 0; k < limit_S; k++) {
		symbol = '.';
		strncat(pp_S, &symbol, 1);
	}
	printf("%s\n", pp_s);
	printf("%s\n", pp_S);
	printf("%s\n", pp_t);
}

void print_matrix(int n, int m, int total_measure_nb, int H_bps[n][m], int s[n], int t[m]) {
	float percentage;
	printf("  ");
	for (int j = 0; j < m; j++) {
		printf("  %c   ", decode_nucleotide(t[j]));
	}
	printf("\n");
	for (int i = 0; i < n; i++) {
		percentage = 100 * (float) H_bps[i][0] / (float) total_measure_nb;
		if (percentage < 10.0) {printf("%c  %.2f", decode_nucleotide(s[n-i-1]), percentage);}
		else {printf("%c %.2f", decode_nucleotide(s[n-i-1]), percentage);}
		for (int j = 1; j < m; j++) {
			percentage = 100 * (float) H_bps[i][j] / (float) total_measure_nb;
			if (percentage < 10.0) {printf("  %.2f", percentage);}
			else {printf(" %.2f", percentage);}
		}
		printf("\n");
	}
}

int hash_S(int n, int S[n], int prime) {
	int total = 21001;
	for (int i = 0; i < n; i++) {
		if (S[i] < INF_5) {
			total = (total * first_primes[max_seq_size - S[i] - 2] + first_primes[S[i]] + S[i] + i + 1) % prime;
		} else {
			total = (total * first_primes[max_seq_size - i - 1] + first_primes[i] + S[i] + i + 2) % prime;
		}
	}
	return total;
}

int hash_H(int n, unsigned short int S[n+1], int prime) {
	int total = 21001;
	for (int i = 1; i <= n; i++) {
		if (S[i] < INF_5) {
			total = (total * first_primes[max_seq_size - S[i] - 2] + first_primes[S[i]] + S[i] + i) % prime;
		} else {
			total = (total * first_primes[max_seq_size - i] + first_primes[i - 1] + S[i] + i + 1) % prime;
		}
	}
	return total;
}

void init_hash_inverse() {
	for (int p = 0; p < max_possible_prime; p++) {
		hash_inverse[p] = -1;
	}
}

int determine_chosen_prime(int n, unsigned short int H_str[mcff_structure_nb][n+1]) {
	int prime, h;
	_Bool collision;
	for (int p = 0; p < possible_prime_nb; p++) {
		prime = possible_primes[p];
		collision = 0;
		init_hash_inverse();
		for (int k = 0; k < mcff_structure_nb; k++) {
			h = hash_H(n, H_str[k], prime);
			if (hash_inverse[h] > -1) {
				collision = 1;
				break;
			} else {
				hash_inverse[h] = k;
			}
		}
		if (!collision) {
			return prime;
		}
	}
	return 1;
}

int get_max_candidate_nb(int n, int m) {
	if (n > m) {return 14 * (1 + n / 4);}
	return 14 * (1 + m / 4);
}

int determine_max_job_nb(int n, int m) {
	int sum = n + m;
	if (sum <= 30) {
		return 256;
	} else if (sum < 60) {
		return 128;
	} else if (sum < 90) {
		return 100;
	} else {
		return 50;
	}
}

int initialize_candidates(int n, int m, int max_candidate_nb, int S[n], int candidates[max_candidate_nb][2]) {
	int candidate_nb = 0;
	for (int i = 0; i < n; i++) {
		if (S[i] < INF_5) {
			candidates[candidate_nb][0] = i;
			candidates[candidate_nb][1] = S[i];
			candidate_nb++;
		}
	}
	int i, j, a, b, u, v;
	i = candidates[candidate_nb - 1][0];
	j = candidates[candidate_nb - 1][1];
	for (int k = 0; k < ncm_nb; k++) {
		u = i + ncm_gaps[k][0];
		v = j + ncm_gaps[k][1];
		if (u < n && v < m) {
			candidates[candidate_nb][0] = u;
			candidates[candidate_nb][1] = v;
			candidate_nb++;
		}
	}
	i = candidates[0][0];
	j = candidates[0][1];
	for (int k = 0; k < ncm_nb; k++) {
		a = i - ncm_gaps[k][0];
		b = j - ncm_gaps[k][1];
		if (a > -1 && b > -1) {
			candidates[candidate_nb][0] = a;
			candidates[candidate_nb][1] = b;
			candidate_nb++;
		}
	}
	a = -1;
	b = -1;
	for (int u = 0; u < n; u++) {
		if (S[u] < INF_5) {
			if (a > -1 && b > -1) {
				for (int k = 0; k < ncm_nb; k++) {
					i = a + ncm_gaps[k][0];
					j = b + ncm_gaps[k][1];
					if (i < u && j < S[u]) {
						candidates[candidate_nb][0] = i;
						candidates[candidate_nb][1] = j;
						candidate_nb++;
					}
				}
			}
			a = u;
			b = S[u];
		}
	}
	u = n;
	v = m;
	for (int a = n-1; a > -1; a--) {
		if (S[a] < INF_5) {
			if (u < n && v < m) {
				for (int k = 0; k < ncm_nb; k++) {
					i = u - ncm_gaps[k][0];
					j = v - ncm_gaps[k][1];
					if (i > a && j > S[a] && !is_ncm(i - a, j - S[a])) {
						candidates[candidate_nb][0] = i;
						candidates[candidate_nb][1] = j;
						candidate_nb++;
					}
				}
			}
			u = a;
			v = S[a];
		}
	}
	return candidate_nb;
}

double get_pairing_probability(int n, int m, int i, int j, int s[n], int t[m], int S[n], double RT) {
	int c = -INF_4, d = -INF_4, a = -INF_3, b = -INF_3, u = INF_3, v = INF_3, w = INF_4, z = INF_4;
	_Bool ab_defined = 0, uv_defined = 0;
	for (int k = i - 1; k > -1; k--) {
		if (S[k] < INF_5) {
			if (ab_defined) {
				c = k;
				d = S[k];
				break;
			} else {
				a = k;
				b = S[k];
				ab_defined = 1;
			}
		}
	}
	for (int k = i + 1; k < n; k++) {
		if (S[k] < INF_5) {
			if (uv_defined) {
				w = k;
				z = S[k];
				break;
			} else {
				u = k;
				v = S[k];
				uv_defined = 1;
			}
		}
	}
	int gap1_p = a - c, gap1_q = b - d;
	int gap2_p = i - a, gap2_q = j - b;
	int gap3_p = u - i, gap3_q = v - j;
	int gap4_p = w - u, gap4_q = z - v;
	int outr_p = u - a, outr_q = v - b;
    int seq2_id, seq3_id, outr_seq_id;
	int gap1_id, gap2_id, gap3_id, gap4_id, outr_gap_id;
	double energy = 0.0;
	if (is_ncm(gap2_p, gap2_q)) {
        seq2_id = get_seq_index(s[n - 1 - i], s[n - 1 - a], t[b], t[j]);
        gap2_id = ncm_table[gap2_p][gap2_q];
		gap1_id = is_ncm(gap1_p, gap1_q) ? ncm_table[gap1_p][gap1_q] : default_gap_id;
        energy = leftE[gap1_id][gap2_id][seq2_id];
        if (is_ncm(gap3_p, gap3_q)) {
			seq3_id = get_seq_index(s[n - 1 - u], s[n - 1 - i], t[j], t[v]);
			gap3_id = ncm_table[gap3_p][gap3_q];
			gap4_id = is_ncm(gap4_p, gap4_q) ? ncm_table[gap4_p][gap4_q] : default_gap_id;
			energy += stackE[gap2_id][gap3_id][gap4_id][seq3_id];
			if (is_ncm(outr_p, outr_q)) {
				outr_seq_id = get_seq_index(s[n - 1 - u], s[n - 1 - a], t[b], t[v]);
				outr_gap_id = ncm_table[outr_p][outr_q];
				energy -= stackE[gap1_id][outr_gap_id][gap4_id][outr_seq_id];
			}
			return sigmoid_of_negative(energy / RT);
		}
        return sigmoid_of_negative(energy / RT);
    }
	if (is_ncm(gap3_p, gap3_q)) {
		seq3_id = get_seq_index(s[n - 1 - u], s[n - 1 - i], t[j], t[v]);
		gap3_id = ncm_table[gap3_p][gap3_q];
		gap4_id = is_ncm(gap4_p, gap4_q) ? ncm_table[gap4_p][gap4_q] : default_gap_id;
        energy = rightE[gap3_id][gap4_id][seq3_id];
        return sigmoid_of_negative(energy / RT);
	}
    return 0.5;
}

void simulate(int job_nb, int sim_nb, int n, int m, int chosen_prime, int iteration_nb, int bp_nb, int s[n], int t[m], int init_S[n], int H_bpss[max_job_nb][n][m], unsigned short int H_strs[max_job_nb][mcff_structure_nb][n+1], int H_strfs[max_job_nb][mcff_structure_nb], int transition_nbs[max_job_nb]) {
	const int64_t m1 = INT64_C(4294967087);
	const int64_t m2 = INT64_C(4294944443);
	const int32_t a12 = INT32_C(1403580);
	const int32_t a13 = INT32_C(810728);
	const int32_t a21 = INT32_C(527612);
	const int32_t a23 = INT32_C(1370589);
	const int64_t corr1 = (m1 * a13);
	const int64_t corr2 = (m2 * a23);
	const double norm = 0x1.000000d00000bp-32;
	int64_t s10 = __MRG32k3a_s10[sim_nb];
	int64_t s11 = __MRG32k3a_s11[sim_nb];
	int64_t s12 = __MRG32k3a_s12[sim_nb];
	int64_t s20 = __MRG32k3a_s20[sim_nb];
	int64_t s21 = __MRG32k3a_s21[sim_nb];
	int64_t s22 = __MRG32k3a_s22[sim_nb];
	int64_t p, r;
	
	int transition_nb = 0;
	double inverse_start_RT = 1 / start_RT;
	double inverse_end_RT = 1 / end_RT;
	double diff_inverse_RT = inverse_end_RT - inverse_start_RT;
	double RT_step = diff_inverse_RT / (iteration_nb - 1);
	int max_candidate_nb = get_max_candidate_nb(n, m);
	int candidates[max_candidate_nb][2];
	for (int i = 0; i < max_candidate_nb; i++) {
		for (int j = 0; j < 2; j++) {
			candidates[i][j] = 0;
		}
	}
	int S[n];
	for (int i = 0; i < n; i++) {
		S[i] = init_S[i];
	}
	int candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
	int candidate, i, j, to_swap, frequency, hash, hash_idx;
	_Bool paired, was_paired, found_match, hash_match;
	double RT, pairing_probability;
	int last_frequency = 1;
	for (int iteration = 0; iteration < iteration_nb; iteration++) {
		r = s12 - s22;
		r -= m1 * ((r - 1) >> 63);
		p = (a12 * s11 - a13 * s10 + corr1) % m1;
		s10 = s11;
		s11 = s12;
		s12 = p;
		p = (a21 * s22 - a23 * s20 + corr2) % m2;
		s20 = s21;
		s21 = s22;
		s22 = p;
		candidate = floor(candidate_nb * r * norm);
		if (bp_nb == 1) {
			r = s12 - s22;
			r -= m1 * ((r - 1) >> 63);
			p = (a12 * s11 - a13 * s10 + corr1) % m1;
			s10 = s11;
			s11 = s12;
			s12 = p;
			p = (a21 * s22 - a23 * s20 + corr2) % m2;
			s20 = s21;
			s21 = s22;
			s22 = p;
			candidate = floor(1 + (candidate_nb - 1) * r * norm);
		}
		i = candidates[candidate][0];
		j = candidates[candidate][1];
		was_paired = S[i] == j;
		RT = 1.0 / (inverse_start_RT + RT_step * iteration);
		pairing_probability = get_pairing_probability(n, m, i, j, s, t, S, RT);
		r = s12 - s22;
		r -= m1 * ((r - 1) >> 63);
		p = (a12 * s11 - a13 * s10 + corr1) % m1;
		s10 = s11;
		s11 = s12;
		s12 = p;
		p = (a21 * s22 - a23 * s20 + corr2) % m2;
		s20 = s21;
		s21 = s22;
		s22 = p;
		paired = r * norm < pairing_probability;
		if (paired != was_paired) {
			transition_nb++;
			for (int k = 0; k < n; k++) {
				if (S[k] < INF_5) {H_bpss[job_nb][k][S[k]] = H_bpss[job_nb][k][S[k]] + last_frequency;}
			}
			hash = hash_S(n, S, chosen_prime);
			hash_idx = hash_inverse[hash];
			if (hash_idx > -1) {
				hash_match = 1;
				for (int l = 0; l < n; l++) {
					if (S[l] != H_strs[job_nb][hash_idx][l+1]) {
						hash_match = 0;
						break;
					}
				}
				if (hash_match) {
					H_strfs[job_nb][hash_idx] = H_strfs[job_nb][hash_idx] + last_frequency;
				}
			}
			last_frequency = 0;
			if (paired && !was_paired) {
				S[i] = j;
				bp_nb++;
				candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
			}
			else {
				S[i] = INF_5;
				bp_nb--;
				candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
			}
		}
		last_frequency++;
	}
	for (int k = 0; k < n; k++) {
		if (S[k] < INF_5) {H_bpss[job_nb][k][S[k]] = H_bpss[job_nb][k][S[k]] + last_frequency;}
	}
	hash = hash_S(n, S, chosen_prime);
	hash_idx = hash_inverse[hash];
	if (hash_idx > -1) {
		hash_match = 1;
		for (int l = 0; l < n; l++) {
			if (S[l] != H_strs[job_nb][hash_idx][l+1]) {
				hash_match = 0;
				break;
			}
		}
		if (hash_match) {
			H_strfs[job_nb][hash_idx] = H_strfs[job_nb][hash_idx] + last_frequency;
		}
	}
	transition_nbs[job_nb] = transition_nb;
}

void quick_simulate(int job_nb, int sim_nb, int n, int m, int iteration_nb, int bp_nb, int s[n], int t[m], int init_S[n], int H_bpss[max_job_nb][n][m], int transition_nbs[max_job_nb]) {
	const int64_t m1 = INT64_C(4294967087);
	const int64_t m2 = INT64_C(4294944443);
	const int32_t a12 = INT32_C(1403580);
	const int32_t a13 = INT32_C(810728);
	const int32_t a21 = INT32_C(527612);
	const int32_t a23 = INT32_C(1370589);
	const int64_t corr1 = (m1 * a13);
	const int64_t corr2 = (m2 * a23);
	const double norm = 0x1.000000d00000bp-32;
	int64_t s10 = __MRG32k3a_s10[sim_nb];
	int64_t s11 = __MRG32k3a_s11[sim_nb];
	int64_t s12 = __MRG32k3a_s12[sim_nb];
	int64_t s20 = __MRG32k3a_s20[sim_nb];
	int64_t s21 = __MRG32k3a_s21[sim_nb];
	int64_t s22 = __MRG32k3a_s22[sim_nb];
	int64_t p, r;
	
	int transition_nb = 0;
	double inverse_start_RT = 1 / start_RT;
	double inverse_end_RT = 1 / end_RT;
	double diff_inverse_RT = inverse_end_RT - inverse_start_RT;
	double RT_step = diff_inverse_RT / (iteration_nb - 1);
	int max_candidate_nb = get_max_candidate_nb(n, m);
	int candidates[max_candidate_nb][2];
	for (int i = 0; i < max_candidate_nb; i++) {
		for (int j = 0; j < 2; j++) {
			candidates[i][j] = 0;
		}
	}
	int S[n];
	for (int i = 0; i < n; i++) {
		S[i] = init_S[i];
	}
	int candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
	int candidate, i, j, to_swap, frequency;
	_Bool paired, was_paired, found_match;
	double RT, pairing_probability;
	int last_frequency = 1;
	for (int iteration = 0; iteration < iteration_nb; iteration++) {
		r = s12 - s22;
		r -= m1 * ((r - 1) >> 63);
		p = (a12 * s11 - a13 * s10 + corr1) % m1;
		s10 = s11;
		s11 = s12;
		s12 = p;
		p = (a21 * s22 - a23 * s20 + corr2) % m2;
		s20 = s21;
		s21 = s22;
		s22 = p;
		candidate = floor(candidate_nb * r * norm);
		if (bp_nb == 1) {
			r = s12 - s22;
			r -= m1 * ((r - 1) >> 63);
			p = (a12 * s11 - a13 * s10 + corr1) % m1;
			s10 = s11;
			s11 = s12;
			s12 = p;
			p = (a21 * s22 - a23 * s20 + corr2) % m2;
			s20 = s21;
			s21 = s22;
			s22 = p;
			candidate = floor(1 + (candidate_nb - 1) * r * norm);
		}
		i = candidates[candidate][0];
		j = candidates[candidate][1];
		was_paired = S[i] == j;
		RT = 1.0 / (inverse_start_RT + RT_step * iteration);
		pairing_probability = get_pairing_probability(n, m, i, j, s, t, S, RT);
		r = s12 - s22;
		r -= m1 * ((r - 1) >> 63);
		p = (a12 * s11 - a13 * s10 + corr1) % m1;
		s10 = s11;
		s11 = s12;
		s12 = p;
		p = (a21 * s22 - a23 * s20 + corr2) % m2;
		s20 = s21;
		s21 = s22;
		s22 = p;
		paired = r * norm < pairing_probability;
		if (paired != was_paired) {
			transition_nb++;
			for (int k = 0; k < n; k++) {
				if (S[k] < INF_5) {H_bpss[job_nb][k][S[k]] = H_bpss[job_nb][k][S[k]] + last_frequency;}
			}
			last_frequency = 0;
			if (paired && !was_paired) {
				S[i] = j;
				bp_nb++;
				candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
			}
			else {
				S[i] = INF_5;
				bp_nb--;
				candidate_nb = initialize_candidates(n, m, max_candidate_nb, S, candidates);
			}
		}
		last_frequency++;
	}
	for (int k = 0; k < n; k++) {
		if (S[k] < INF_5) {H_bpss[job_nb][k][S[k]] = H_bpss[job_nb][k][S[k]] + last_frequency;}
	}
	transition_nbs[job_nb] = transition_nb;
} 

void dynamic_fold(int n, int m, int total_measure_nb, int H_bps[n][m], int final_S[n]) {
	int table[n][m][3];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			table[i][j][0] = 0;
			table[i][j][1] = -1;
			table[i][j][2] = -1;
		}
	}
	for (int i = 0; i < n; i++) {
		table[i][m-1][0] = H_bps[i][m-1];
	}
	for (int j = 0; j < m - 1; j++) {
		table[n-1][j][0] = H_bps[n-1][j];
	}
	int max_score = 0;
	int best_u = -1, best_v = -1;
	int u_limit, v_limit;
	for (int i = n-2; i > -1; i--) {
		for (int j = m-2; j > -1; j--) {
			max_score = table[i+1][j+1][0];
			best_u = i + 1;
			best_v = j + 1;
			u_limit = minimum(i + B + 2, n);
			v_limit = minimum(j + B + 2, m);
			for (int u = i+2; u < u_limit; u++) {
				if (table[u][j+1][0] > max_score) {
					max_score = table[u][j+1][0];
					best_u = u;
					best_v = j + 1;
				}
			}
			for (int v = j+2; v < v_limit; v++) {
				if (table[i+1][v][0] > max_score) {
					max_score = table[i+1][v][0];
					best_u = i + 1;
					best_v = v;
				}
			}
			table[i][j][0] = max_score + H_bps[i][j];
			table[i][j][1] = best_u;
			table[i][j][2] = best_v;
		}
	}
	max_score = table[0][0][0];
	int start_i = 0, start_j = 0;
	int score;
	for (int i = 1; i < n; i++) {
		score = table[i][0][0];
		if (score > max_score) {
			max_score = table[i][0][0];
			start_i = i;
			start_j = 0;
		}
	}
	for (int j = 1; j < m; j++) {
		score = table[0][j][0];
		if (score > max_score) {
			max_score = table[0][j][0];
			start_i = 0;
			start_j = j;
		}
	}
	score = H_bps[start_i][start_j];
	if (score > alpha * total_measure_nb) {
		final_S[start_i] = start_j;
	}
	int i = table[start_i][start_j][1];
	int j = table[start_i][start_j][2];
	int temp_i;
	while (i > -1) {
		score = H_bps[i][j];
		if (score > alpha * total_measure_nb) {
			final_S[i] = j;
		}
		temp_i = i;
		i = table[temp_i][j][1];
		j = table[temp_i][j][2];
	}
}

double compute_entropy(int n, int total_measure_nb, int H_strf[mcff_structure_nb]) {
	double entropy = 0.0;
	double probability = 0.0;
	for (int k = 0; k < entropy_structure_nb; k++) {
		if (H_strf[k] > 0) {
			probability = (float) H_strf[k] / (float) total_measure_nb;
			entropy -= probability * log(probability);
		}
	}
	return entropy;
}

void print_help() {
	printf("MC-DuplexFold infers duplex secondary structure from two RNA sequences.\n\n");
	printf("Input parameters:\n\n");
	printf("  -s\n");
	printf("   Mandatory, string (maximum 100 characters)\n");
	printf("   Specifies the first RNA sequence from 5' to 3' with the alphabet acgtuACGTU\n\n");
	printf("  -t\n");
	printf("   Mandatory, string (maximum 100 characters)\n");
	printf("   Specifies the second RNA sequence from 5' to 3' with the alphabet acgtuACGTU\n\n");
	printf("Execution parameters:\n\n");
	printf("  -i\n");
	printf("   Optional, String, default=None\n");
	printf("   If not None, uses this input structure as a warm start for the simulation.\n");
	printf("   Structure must be written in dot-bracket notation without spaces inside quotation marks. e.g. \'(((..(()).)))\'\n");
	printf("   If None, initiates the simulation according to parameter w.\n\n");
	printf("  -w\n");
	printf("   Optional, boolean, default=1\n");
	printf("   If 1, uses mcff's MFE structure as a warm start for the simulation. (Requires mcff to be executable in PATH.)\n");
	printf("   If 0, only pairs the first nucleotide of each strand to start the simulation.\n\n");
	printf("  -n\n");
	printf("   Optional, integer, default=100\n");
	printf("   Specifies the number of simulations that are executed.\n\n");
	printf("  -f\n");
	printf("   Optional, integer, default=0\n");
	printf("   Specifies the number of suboptimal structures for which we print the frequency.\n\n");
	printf("  -l\n");
	printf("   Optional, integer, default=10\n");
	printf("   Controls the number of steps in each simulation, which is given by l*((number of nucleotides in s and t)/2)^2\n\n");
	printf("  -a\n");
	printf("   Optional, float, default=0.45\n");
	printf("   Base pairs with an appearance frequency lower than this parameter are removed from the final structure.\n\n");
	printf("  -b\n");
	printf("   Optional, integer, default=5\n");
	printf("   Specifies the maximum bulge length in the final structure computation.\n\n");
	printf("  -u\n");
	printf("   Optional, float, default=3.0\n");
	printf("   Specifies the initial temperature of the simulation.\n\n");
	printf("  -v\n");
	printf("   Optional, float, default=0.333333333\n");
	printf("   Specifies the final temperature of the simulation.\n\n");
	printf("  -r\n");
	printf("   Optional, int, default=1\n");
	printf("   Specifies the seed for the random number generation.\n\n");
	printf("  -p\n");
	printf("   Optional, boolean, default=0\n");
	printf("   If 1, prints the base pair frequency matrix.\n");
	printf("  -h\n");
	printf("   Display the usage details message.\n");
}

int main(int argc, char *argv[]) {
	int return_code = 0;
	char c;
	char *endptr;
	int n = 0, m = 0, init_dotb_length = 0;
	char seq_1[max_seq_size], seq_2[max_seq_size], init_dotb[2 * max_seq_size];
	_Bool changed_s = 0, changed_t = 0, changed_i = 0;
	while ((c = getopt(argc, argv, "hs:t:i:w:n:f:r:l:a:b:u:v:p:")) != -1) {
		switch (c) {
		case 'h':
			print_help();
			return 0;
		case 's':
			n = strlen(optarg);
			for (int i = 0; i < minimum(n, max_seq_size - 1); i++) {seq_1[i] = optarg[i];}
			seq_1[n] = '\0';
			changed_s = 1;
			break;
		case 't':
			m = strlen(optarg);
			for (int i = 0; i < minimum(m, max_seq_size - 1); i++) {seq_2[i] = optarg[i];}
			seq_2[m] = '\0';
			changed_t = 1;
			break;
		case 'i':
			init_dotb_length = strlen(optarg);
			for (int i = 0; i < init_dotb_length; i++) {init_dotb[i] = optarg[i];}
			init_dotb[init_dotb_length] = '\0';
			changed_i = 1;
			break;
		case 'w':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {warm_start = strtol(optarg, &endptr, 10);}
			break;
		case 'n':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {simulation_nb = strtol(optarg, &endptr, 10);}
			break;
		case 'f':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {suboptimal_nb = strtol(optarg, &endptr, 10);}
			break;
		case 'l':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {K = strtol(optarg, &endptr, 10);}
			break;
		case 'a':
			return_code = cannot_cast_str_to_double(optarg, c);
			if (return_code) {return 1;}
			else {alpha = strtod(optarg, &endptr);}
			break;
		case 'b':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {B = strtol(optarg, &endptr, 10);}
			break;
		case 'u':
			return_code = cannot_cast_str_to_double(optarg, c);
			if (return_code) {return 1;}
			else {start_RT = strtod(optarg, &endptr);}
			break;
		case 'v':
			return_code = cannot_cast_str_to_double(optarg, c);
			if (return_code) {return 1;}
			else {end_RT = strtod(optarg, &endptr);}
			break;
		case 'r':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {MRG32k3a_seed = strtol(optarg, &endptr, 10);}
			break;
		case 'p':
			return_code = cannot_cast_str_to_int(optarg, c);
			if (return_code) {return 1;}
			else {print_stats = strtol(optarg, &endptr, 10);}
			break;
		case '?':
			if (optopt == 's')
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			else if (optopt == 't')
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint(optopt))
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			abort();
		}
	}
	if (!changed_s) {
		fprintf(stderr, "First sequence must be specified with option -s\n");
		return 1;
	} else {
		n = strlen(seq_1);
	}
	if (!changed_t) {
		fprintf(stderr, "Second sequence must be specified with option -t\n");
		return 1;
	} else {
		m = strlen(seq_2);
	}
	if (suboptimal_nb > 1 && (n < 6 || m < 6)) {
		fprintf(stderr, "Sequences too small for suboptimal analysis. Avoid -f %ld.\n", suboptimal_nb);
		return 1;
	}
	int s[n];
	int t[m];
	return_code = encode_sequence(n, seq_1, s);
	if (return_code) {return 1;}
	return_code = encode_sequence(m, seq_2, t);
	if (return_code) {return 1;}
	return_code = decode_sequence(n, seq_1, s);
	if (return_code) {return 1;}
	return_code = decode_sequence(m, seq_2, t);
	if (return_code) {return 1;}
	int init_S[n];
	max_job_nb = determine_max_job_nb(n, m);
	unsigned short int H_str[mcff_structure_nb][n+1];
	int H_strf[mcff_structure_nb];
	unsigned short int H_strs[max_job_nb][mcff_structure_nb][n+1];
	int H_strfs[max_job_nb][mcff_structure_nb];
	if (suboptimal_nb > 0) {
		for (int k = 0; k < mcff_structure_nb; k++) {
			for (int i = 0; i < n+1; i++) {
				H_str[k][i] = 0;
			}
			H_strf[k] = 0;
		}
	}
	double deltaGs[mcff_structure_nb];
	if (suboptimal_nb > 0) {
		for (int k = 0; k < mcff_structure_nb; k++) {
			deltaGs[k] = 0.0;
		}
	}
	int bp_nb = 1;
	if (!changed_i) {
		init_dotb[0] = '\0';
		if (warm_start) {
			if (suboptimal_nb > 0) {
				bp_nb = initialize_S(n, m, seq_1, seq_2, init_S, H_str, deltaGs, 1);
				if (bp_nb == -1) {
					fprintf(stderr, "mcff error. Avoid -f %ld.\n", suboptimal_nb);
					return 1;
				}
			} else {
				bp_nb = quick_initialize_S(n, m, seq_1, seq_2, init_S);
			}
		} else {
			init_S[0] = 0;
			for (int i = 1; i < n; i++) {init_S[i] = INF_5;}
		}
	} else {
		if (strlen(init_dotb) != n+m) {
			fprintf(stderr, "Invalid initiation dot-bracket structure.\n");
			return 1;
		} else {
			if (suboptimal_nb > 0) {
				bp_nb = initialize_S(n, m, seq_1, seq_2, init_S, H_str, deltaGs, 1);
				if (bp_nb == -1) {
					fprintf(stderr, "mcff error. Avoid -f %ld.\n", suboptimal_nb);
					return 1;
				}
			} else {
				bp_nb = quick_initialize_S(n, m, seq_1, seq_2, init_S);
			}
			bp_nb = initialize_S_from_dotb(n, m, seq_1, seq_2, init_dotb, init_S);
		}
	}
	double power = 2.0;
	double base = (n + m) / 2.0;
	int iteration_nb = K * pow(base, power);
	int total_iteration_nb = iteration_nb * simulation_nb;
	int total_measure_nb = total_iteration_nb + simulation_nb;
	int H_bps[n][m];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			H_bps[i][j] = 0;
		}
	}
	int H_bpss[max_job_nb][n][m];
	for (int job_nb = 0; job_nb < max_job_nb; job_nb++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				H_bpss[job_nb][i][j] = 0;
			}
		}
	}
	
	int chosen_prime;
	if (suboptimal_nb > 0) {
		chosen_prime = determine_chosen_prime(n, H_str);
		if (chosen_prime == 1) {
			fprintf(stderr, "No efficient hash map could be computed.\n");
			return 1;
		}
		for (int k = 0; k < mcff_structure_nb; k++) {
			H_str[k][0] = hash_H(n, H_str[k], chosen_prime);
		}
	}
	
	if (suboptimal_nb > 0) {
		for (int job_nb = 0; job_nb < max_job_nb; job_nb++) {
			for (int k = 0; k < mcff_structure_nb; k++) {
				for (int i = 0; i < n+1; i++) {
					H_strs[job_nb][k][i] = H_str[k][i];
				}
				H_strfs[job_nb][k] = 0;
			}
		}
	}
	
	for (int mrg_nb = 0; mrg_nb < MRG_SIZE; mrg_nb++) {
		MRG32k3a_init(MRG32k3a_seed + mrg_nb, mrg_nb);
		MRG32k3a_init6(__MRG32k3a_s10[mrg_nb], __MRG32k3a_s11[mrg_nb], __MRG32k3a_s12[mrg_nb], __MRG32k3a_s20[mrg_nb], __MRG32k3a_s21[mrg_nb], __MRG32k3a_s22[mrg_nb], mrg_nb);
	}
	int transition_nb = 0;
	int transition_nbs[max_job_nb];
	for (int job_nb = 0; job_nb < max_job_nb; job_nb++) {transition_nbs[job_nb] = 0;}
	int sim_nb = 0, jobs_to_do = simulation_nb, job_batch_nb = 0;
	int job_nb_ub = 0;
	if (suboptimal_nb > 0) {
		while (jobs_to_do > 0) {
			job_nb_ub = minimum(jobs_to_do, max_job_nb);
			#pragma omp parallel for
			for (int job_nb = 0; job_nb < job_nb_ub; job_nb++) {
				sim_nb = job_batch_nb * max_job_nb + job_nb;
				simulate(job_nb, sim_nb, n, m, chosen_prime, iteration_nb, bp_nb, s, t, init_S, H_bpss, H_strs, H_strfs, transition_nbs);
			}
			for (int job_nb = 0; job_nb < job_nb_ub; job_nb++) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						H_bps[i][j] += H_bpss[job_nb][i][j];
						H_bpss[job_nb][i][j] = 0;
					}
				}
				for (int k = 0; k < mcff_structure_nb; k++) {
					H_strf[k] += H_strfs[job_nb][k];
					H_strfs[job_nb][k] = 0;
				}
				transition_nb += transition_nbs[job_nb];
				transition_nbs[job_nb] = 0;
			}
			jobs_to_do -= job_nb_ub;
			job_batch_nb += 1;
		}
	} else {
		while (jobs_to_do > 0) {
			job_nb_ub = minimum(jobs_to_do, max_job_nb);
			#pragma omp parallel for
			for (int job_nb = 0; job_nb < job_nb_ub; job_nb++) {
				sim_nb = job_batch_nb * max_job_nb + job_nb;
				quick_simulate(job_nb, sim_nb, n, m, iteration_nb, bp_nb, s, t, init_S, H_bpss, transition_nbs);
			}
			for (int job_nb = 0; job_nb < job_nb_ub; job_nb++) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						H_bps[i][j] += H_bpss[job_nb][i][j];
						H_bpss[job_nb][i][j] = 0;
					}
				}
				transition_nb += transition_nbs[job_nb];
				transition_nbs[job_nb] = 0;
			}
			jobs_to_do -= job_nb_ub;
			job_batch_nb += 1;
		}
	}
	float transition_frequency = (float) transition_nb / (float) total_measure_nb;
	int final_S[n];
	for (int i = 0; i < n; i++) {
		final_S[i] = INF_5;
	}
	dynamic_fold(n, m, total_measure_nb, H_bps, final_S);
	int frequent_S[n];
	double subopt_percentage, other_percentage;
	int total_structures_noted = 0;
	int printed_subopt[mcff_structure_nb];
	for (int k = 0; k < mcff_structure_nb; k++) {printed_subopt[k] = 0;}
	int best_freq, best_idx, freq_sum;
	double entropy = 0.0;
	if (suboptimal_nb > 0) {
		for (int subopt = 0; subopt < suboptimal_nb; subopt++) {
			best_freq = 0;
			best_idx = -1;
			freq_sum = 0;
			for (int k = 0; k < mcff_structure_nb; k++) {
				if (printed_subopt[k] == 0) {
					freq_sum += H_strf[k];
					if (H_strf[k] > best_freq) {
						best_freq = H_strf[k];
						best_idx = k;
					}
				}
			}
			if (freq_sum == 0) {break;}
			printed_subopt[best_idx] = 1;
			for (int i = 0; i < n; i++) {
				frequent_S[i] = H_str[best_idx][i+1];
			}
			total_structures_noted += H_strf[best_idx];
			subopt_percentage = 100 * (float) H_strf[best_idx] / (float) total_measure_nb;
			printf("Structure %d: %.6f%%\n", subopt + 1, subopt_percentage);
			print_duplex(n, m, frequent_S, seq_1, seq_2);
			printf("\n");
		}
		other_percentage = 100 * (1.0 - ((float) total_structures_noted / (float) total_measure_nb));
		printf("Other structures: %.2f%%\n", other_percentage);
		entropy = compute_entropy(n, total_structures_noted, H_strf);
		printf("Estimated entropy over %d structures: %f\n", mcff_structure_nb, entropy);
		printf("Transition frequency: %f\n", transition_frequency);
	}
	if (print_stats) {print_matrix(n, m, total_measure_nb, H_bps, s, t);}
	print_duplex(n, m, final_S, seq_1, seq_2);
	
    return 0;
}