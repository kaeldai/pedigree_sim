#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <random>
#include <cmath>
#include <sstream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace sim{

namespace po = boost::program_options;

po::options_description _desc;
po::positional_options_description p;

typedef std::array<std::array<std::array<double, 4>, 10>, 4> read_prop_table;


// Parameters for Dirichlet Multinomial
typedef struct gamma_params {
    double pi;
    double phi;
    double epsilon;
    double omega;
} gamma_t;


// Stores all the input parameters
struct arg_t {
  std::string prefix;
  std::string ped;
  std::string region;

  size_t read_depths;
  double theta;
  double ref_weight;
  double mu;
  double mu_somatic;
  double mu_library;
  std::array<double, 4> nuc_freqs;
  gamma_t model_a, model_b;

  std::string input;
} arg;

struct params_t {
  double pi;
  double phi;
  double epsilon;
  double omega;

  read_prop_table read_proportions;

  //params_t(){ }

  //params_t(double _pi, double _phi, double _epsilon, double _omega) :
  //         pi(_pi), phi(_phi), epsilon(_epsilon), omega(_omega) {   }
  
} model_a, model_b;

//params_t model_aa, model_bb;




void parse_options(int argc, char *argv[]) {
  std::string nuc_freqs;
  std::vector<std::string> gamma_str;

  _desc.add_options()
    ("read-depths,d", po::value<size_t>(&arg.read_depths)->default_value(30), "Num of calls at each site.")
    ("prefix", po::value<std::string>(&arg.prefix)->default_value(""), "output file names.")
    ("ped,p", po::value<std::string>(&arg.ped)->default_value(""), "pedigree file.")
    ("region,r", po::value<std::string>(&arg.region)->default_value(""), "reference region.")
    ("theta", po::value<double>(&arg.theta)->default_value(0.001), "the population diversity.")
    ("ref-weight,R", po::value<double>(&arg.ref_weight)->default_value(1.0), "bias for reference.")
    ("gamma", po::value<std::vector<std::string> >(&gamma_str)->default_value({"0.98,0.0005,0.0005,1.04", "0.02,0.075,0.005,1.18"}, ""), "distribution parameters.")
    ("mu", po::value<double>(&arg.mu)->default_value(1e-8), "germline mutation rate.")
    ("mu-somatic", po::value<double>(&arg.mu_somatic)->default_value(1e-8), "somatic mutation rate.")
    ("mu-library", po::value<double>(&arg.mu_library)->default_value(1e-8), "library mutation rate.")
    ("nuc-freqs", po::value<std::string>(&nuc_freqs)->default_value("0.3,0.2,0.2,0.3"), "order frequency.")
    ("input", po::value<std::string>(&arg.input)->default_value(""), "reference file.")
    ;

  p.add("input", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(_desc).positional(p).run(), vm);
  po::notify(vm);

  if(arg.input.empty()) {
    throw std::runtime_error("No input fasta file.");
  }
  
  //if(arg.ped.empty()) {
  //  throw std::runtime_error("No pedigree file.");
  //}

  if(arg.region.empty()) {
    throw std::runtime_error("No region specified.");
  }



  // Convert nuc_freq to list
  std::vector<std::string> nfs;
  boost::split(nfs, nuc_freqs, boost::is_any_of(","));
  for(int a = 0; a < 4; ++a) {
    arg.nuc_freqs[a] = std::stod(nfs[a]);
  }

  // Convert gamma parameter into two lists
  std::vector<std::string> plist1;
  boost::split(plist1, gamma_str[0], boost::is_any_of(","));
  arg.model_a.pi = std::stod(plist1[0]);
  arg.model_a.phi = std::stod(plist1[1]);
  arg.model_a.epsilon = std::stod(plist1[2]);
  arg.model_a.omega = std::stod(plist1[3]);

  std::vector<std::string> plist2;
  boost::split(plist2, gamma_str[1], boost::is_any_of(","));
  arg.model_b.pi = std::stod(plist2[0]);
  arg.model_b.phi = std::stod(plist2[1]);
  arg.model_b.epsilon = std::stod(plist2[2]);
  arg.model_b.omega = std::stod(plist2[3]);

}





typedef uint8_t Base;

#define char2base(x) nt4_table[x]
Base A = 0;
Base C = 1;
Base G = 2;
Base T = 3;
Base N = 4;

static Base nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, A, 4, C,  4, 4, 4, G,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  T, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


#define base2char(x) base_lookup[x]
static char base_lookup[4] = {'A', 'C', 'G', 'T'};

#define base2str(x) base_str_lookup[x]
static const char *base_str_lookup[4] = {"A", "C", "G", "T"};


enum Genotype : uint8_t {AA = 0, AC, AG, AT, CC, CG, CT, GG, GT, TT, NN}; 


#define gt2char(x) gtchar_table[x]

static const char *gtchar_table[11] = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "NN"};

#define gt_contains(base, gt) base_gt_table[base][gt]

static int base_gt_table[4][11] = {
  {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, // A
  {0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0}, // C
  {0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0}, // G
  {0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0}  // T
};



#define bases2gt(b1,b2) gtlookup[b1][b2]
Genotype gtlookup[4][4] = {
  {AA, AC, AG, AT},
  {AC, CC, CG, CT},
  {AG, CG, GG, GT},
  {AT, CT, GT, TT}
};


#define gt2bases(x) split_table[x]
Base split_table[10][2] = {
  {A,A}, {A,C}, {A,G}, {A,T}, {C,C}, {C,G}, {C,T}, {G,G}, {G,T}, {T,T}
};


// TODO: Use RowMajor order
typedef Eigen::MatrixXd TransitionMatrix;
typedef Eigen::Matrix4d MutationMatrix;

TransitionMatrix mitosis_haploid_matrix(double mu, std::array<double, 4> nuc_freq) {
  double beta = 1.0/(1.0 - pow(nuc_freq[0], 2) - pow(nuc_freq[1], 2)
		     - pow(nuc_freq[2],2) - pow(nuc_freq[3], 2));
  double p = exp(-beta*mu);
    
  TransitionMatrix P(4, 4);
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      P(i,j) = nuc_freq[j]*(1.0-p);
    }
    P(i,i) += p;
  }

  return P;
}

constexpr int folded_diploid_nucleotides[10][2] = {{0, 0}, {0, 1}, {0, 2}, {0, 3},
    {1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}
};

TransitionMatrix meiosis_haploid_matrix(double mu, std::array<double, 4> nuc_freq) {
  TransitionMatrix T = mitosis_haploid_matrix(mu, nuc_freq);
  TransitionMatrix M(10, 4);
  for(int i = 0; i < 10; ++i) {
    int a = folded_diploid_nucleotides[i][0];
    int b = folded_diploid_nucleotides[i][1];
    for(int j = 0; j < 4; ++j) {
      M(i,j) = 0.5*(T(a,j) + T(b,j));
    }
  }
  
  return M;
}


constexpr int folded_diploid_genotypes[16] = {0, 1, 2, 3, 1, 4, 5, 6, 2, 5, 7, 8, 3, 6, 8, 9};

TransitionMatrix meiosis_diploid_matrix(double mu, std::array<double, 4> nuc_freq) {
  TransitionMatrix P1 = meiosis_haploid_matrix(mu, nuc_freq);
  TransitionMatrix P2 = TransitionMatrix(P1);

  // TODO: If there is no diff between mom & dad we can reduce num of rows
  TransitionMatrix K = kroneckerProduct(P1, P2);
  
  TransitionMatrix M = TransitionMatrix::Zero(100, 10);
  for(int i = 0; i < 100; ++i) {
    for(int j = 0; j < 16; ++j) {
      int k = folded_diploid_genotypes[j];
      M(i,k) += K(i,j);
    }
  }
  return M;
}



size_t findGTIndex(Genotype p1, Genotype p2){
  return p1*10 + p2;

}


TransitionMatrix mitosis_diploid_matrix(double mu_soma, std::array<double, 4> nuc_freq) {
  TransitionMatrix F81 = mitosis_haploid_matrix(mu_soma, nuc_freq);
  TransitionMatrix T(10, 10);

  for(size_t i = 0; i < 10; ++i) {
    Base *row = split_table[i];
    Base i1 = row[0];
    Base i2 = row[1];
    
    for(int j = 0; j < 10; ++j) {
      Base *col = split_table[j];
      Base j1 = col[0];
      Base j2 = col[1];
      
      if(i1 == i2 && j1 == j2)
	T(i,j) = F81(i1,j1)*F81(i1,j1);
      else
	T(i,j) = F81(i1,j1)*F81(i2,j2) + F81(i1,j2)*F81(i2,j1);
    }
  }

  return T;
}

 
const gsl_rng *r;
void initGSL() {
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}


typedef std::array<unsigned int, 4> reads_list; 


void parse_region(std::string &region, std::string &contig, size_t &begin, size_t &end) {

  enum STATE {CONTIG, BEGIN, END};
  STATE cstate = CONTIG;
  std::string contig_txt;
  std::string begin_txt;
  
  for(int a = 0; a < region.size(); ++a) {
    char cchar = region[a];
    switch(cstate){
    case CONTIG:
      if(cchar == ':'){
	cstate = BEGIN;	
      }
      else
	contig += cchar;
      break;
      
    case BEGIN:
      if(cchar == '-'){
	std::stringstream convert(begin_txt);
	convert >> begin;
	cstate = END;
      }
      else if(cchar != ','){
	begin_txt += cchar;
      }
      break;
    }
  }
}


 

} // namespace sim


#endif
