/*  Written in 2019 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

// The ugly names are necessary to avoid polluting the name space.
#define MRG_SIZE 1024
int64_t __MRG32k3a_s10[MRG_SIZE], __MRG32k3a_s11[MRG_SIZE], __MRG32k3a_s12[MRG_SIZE], __MRG32k3a_s20[MRG_SIZE], __MRG32k3a_s21[MRG_SIZE], __MRG32k3a_s22[MRG_SIZE];

static inline uint64_t __MRG32k3a_staffordMix13(uint64_t z) {
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	// This is necessary in Java to guarantee that we return a positive number
	return (z >> 1) ^ (z >> 32);
}

/** Initializes the generator starting from a 64-bit seed. */

void MRG32k3a_init(uint64_t s, int sim_nb) {
	const int64_t m1 = UINT64_C(4294967087);
	const int64_t m2 = UINT64_C(4294944443);

	__MRG32k3a_s10[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m1;
	__MRG32k3a_s11[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m1;
	__MRG32k3a_s12[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m1;
	__MRG32k3a_s20[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m2;
	__MRG32k3a_s21[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m2;
	__MRG32k3a_s22[sim_nb] = __MRG32k3a_staffordMix13(s += 0x9e3779b97f4a7c15) % m2;
}

/** Initializes directly the generator state. The first three values must be
in the range [0..4294967087) and the last three values in the range [0..4294944443). */

void MRG32k3a_init6(const uint32_t s10, const uint32_t s11, const uint32_t s12, const uint32_t s20, const uint32_t s21, const uint32_t s22, int sim_nb) {
	const int64_t m1 = INT64_C(4294967087);
	const int64_t m2 = INT64_C(4294944443);

	if (s10 == 0 && s11 == 0 && s12 == 0) {
		fprintf(stderr, "MRG32k3a error\n");
		abort();
	}
	if (s20 == 0 && s21 == 0 && s22 == 0) {
		fprintf(stderr, "MRG32k3a error\n");
		abort();
	}
	if (s10 >= m1 || s11 >= m1 || s12 >= m1) {
		fprintf(stderr, "MRG32k3a error\n");
		abort();
	}

	if (s20 >= m2 || s21 >= m2 || s22 >= m2) {
		fprintf(stderr, "MRG32k3a error\n");
		abort();
	}

	__MRG32k3a_s10[sim_nb] = s10;
	__MRG32k3a_s11[sim_nb] = s11;
	__MRG32k3a_s12[sim_nb] = s12;
	__MRG32k3a_s20[sim_nb] = s20;
	__MRG32k3a_s21[sim_nb] = s21;
	__MRG32k3a_s22[sim_nb] = s22;
}