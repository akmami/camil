#include "lps.h"
#include <fstream>
#include <iostream>
#include <set>
#include "kseq.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

// g++ -std=c++11 -O3 -Wall -Wextra -Wpedantic -mavx2 -I/home/akmuhammet/programs/lcpan/lcptools/include -c lammiq.cpp
// g++ -std=c++11 -O3 -Wall -Wextra -Wpedantic -o lammiq lammiq.o -L/home/akmuhammet/programs/lcpan/lcptools/lib -llcptools -Wl,-rpath,/home/akmuhammet/programs/lcpan/lcptools/lib -lz
// rm lammiq.o

typedef struct {
	std::set<ulabel> distinct_cores1; // forward
	std::set<ulabel> distinct_cores2; // reverse complement
	uint64_t read_count;
} genome_t;

int check_file(const char *filename) {
	std::ifstream file(filename);
	if (!file.good()) {
		std::cerr << "Error opening: " << filename << std::endl;
		return 0;
	}
	return 1;
}

void process_chr(std::string &sequence, genome_t &genome, int lcp_level) {

	struct lps str;
	// forward
	init_lps(&str, sequence.c_str(), sequence.size());
	lps_deepen(&str, lcp_level);

	for (int index = 0; index < str.size; index++) {
		genome.distinct_cores1.insert(str.cores[index].label);
	}

	free_lps(&str);

	// reverse complement
	init_lps2(&str, sequence.c_str(), sequence.size());
	lps_deepen(&str, lcp_level);

	for (int index = 0; index < str.size; index++) {
		genome.distinct_cores2.insert(str.cores[index].label);
	}

	free_lps(&str);

	sequence.clear();
}

void process_ref(const char *filename, genome_t &genome, int lcp_level) {

	genome.read_count = 0;

	// variables
	std::string line;

	std::fstream file;
	file.open(filename, std::ios::in);

	if (file.is_open()) {

		std::string sequence, id;
		sequence.reserve(250000000);

		while (getline(file, line)) {

			if (line[0] == '>') {

				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					process_chr(sequence, genome, lcp_level);
				}

				id = line.substr(1);
				continue;
			}
			else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			process_chr(sequence, genome, lcp_level);
		}

		file.close();
	}

	printf("[INFO] Processed %s\n", filename);
}

uint64_t process_reads(const char *filename, genome_t &genome1, genome_t &genome2, genome_t &genome3, int lcp_level) {

	genome1.read_count = 0;
	genome2.read_count = 0;
	genome3.read_count = 0;
	uint64_t unassigned = 0;

	gzFile file = gzopen(filename, "r");
	if (file == NULL)
		return 0LU;

	kseq_t *read = kseq_init(file);

	while (kseq_read(read) > 0) {
		struct lps str;
		init_lps(&str, read->seq.s, read->seq.l);
		lps_deepen(&str, lcp_level);

		int cores1_count = 0, cores2_count = 0, cores3_count = 0;

		for (int i = 0; i < str.size; i++) {
			if (genome1.distinct_cores1.find(str.cores[i].label) != genome1.distinct_cores1.end() || genome1.distinct_cores2.find(str.cores[i].label) != genome1.distinct_cores2.end())
				cores1_count++;
			if (genome2.distinct_cores1.find(str.cores[i].label) != genome2.distinct_cores1.end() || genome2.distinct_cores2.find(str.cores[i].label) != genome2.distinct_cores2.end())
				cores2_count++;
			if (genome3.distinct_cores1.find(str.cores[i].label) != genome3.distinct_cores1.end() || genome3.distinct_cores2.find(str.cores[i].label) != genome3.distinct_cores2.end())
				cores3_count++;
		}

		free_lps(&str);

		if (cores1_count && cores2_count == 0 && cores3_count == 0)
			genome1.read_count++;
		else if (cores2_count && cores1_count == 0 && cores3_count == 0)
			genome2.read_count++;
		else if (cores3_count && cores1_count == 0 && cores2_count == 0)
			genome3.read_count++;
		else
			unassigned++;
	}

	kseq_destroy(read);
	gzclose(file);

	return unassigned;
}

int main(int argc, char **argv) {

	if (argc < 6) {
		std::cerr << "Wrong format: " << argv[0] << " [genome-1] [genome-2] [genome-3] [lcp-level] ..." << std::endl;
		return -1;
	}

	if (!check_file(argv[1]) || !check_file(argv[2]) || !check_file(argv[3])) {
		return -1;
	}

	int lcp_level = atoi(argv[4]);

	// section 1
	genome_t genome1, genome2, genome3;

	// initializing coefficients of the alphabet and hash tables
	LCP_INIT();

	process_ref(argv[1], genome1, lcp_level);
	process_ref(argv[2], genome2, lcp_level);
	process_ref(argv[3], genome3, lcp_level);

	printf("Core counts: %lu\t%lu\t%lu\n", genome1.distinct_cores1.size(), genome2.distinct_cores1.size(), genome3.distinct_cores1.size());

	for(int i = 5; i < argc; i++) {
		if (!check_file(argv[i])) {
			std::cerr << "File " << argv[i] << " does not exists..." << std::endl;
			continue;
		}

		uint64_t unassigned = process_reads(argv[i], genome1, genome2, genome3, lcp_level);

		printf("Read Stats for %s: \n%lu\t%lu\t%lu\t%lu\n", argv[i], genome1.read_count, genome2.read_count, genome3.read_count, unassigned);
	}
	return 0;
}