#ifndef MODELS_H
#define MODELS_H




#include "simulator.h"
#include "pedigree.h"


namespace sim {
namespace model {


/**
 * Simulation of germline/somatic/library sequences using models used
 * in dng-call and dng-likilhood.
 */
class DNGModel {

 public:
 DNGModel(double ref_weight, double theta, std::array<double,4> nuc_freqs, double mu, double mu_somatic,
	  gamma_t model_a, gamma_t model_b) : ref_weight_(ref_weight), theta_(theta), nuc_freqs_(nuc_freqs), mu_(mu),
    mu_soma_(mu_somatic), model_a_(model_a), model_b_(model_b) {
      // Create probabilities for founder germline based on ref, freqs, and dispersial (theta)
      genotype_priors_[0] = population_prior(theta_, nuc_freqs_, {ref_weight_, 0, 0, 0});
      genotype_priors_[1] = population_prior(theta_, nuc_freqs_, {0, ref_weight_, 0, 0});
      genotype_priors_[2] = population_prior(theta_, nuc_freqs_, {0, 0, ref_weight_, 0});
      genotype_priors_[3] = population_prior(theta_, nuc_freqs_, {0, 0, 0, ref_weight_});

      // Create a set of distributions based on genotype_priors
      setPopPriorDist(genotype_priors_);

      // Create distribution for germline transitions
      germline_trans_matrix_ = meiosis_diploid_matrix(mu_, nuc_freqs_);
      setGermlineDist(germline_trans_matrix_);

      // Distribution for somatic transitons
      somatic_trans_matrix_ = mitosis_diploid_matrix(mu_soma_, nuc_freqs_);
      setSomaticDist(somatic_trans_matrix_);

      // Create two lookup tables for call probabilities
      createReadProportions(model_a_.epsilon, model_a_.omega, read_proportions_a_);
      createReadProportions(model_b_.epsilon, model_b_.omega, read_proportions_b_);
      initGSL();
      initBernoulliGenerator(model_a_.phi); // model_a_.phi + model_b_.phi = 1
      
  }

  void print_priors() {    
    for(int i : {A, C, G, T}) {
      switch(i) {
      case 0: std::cout << "A: "; break;
      case 1: std::cout << "C: "; break; 
      case 2: std::cout << "G: "; break; 
      case 3: std::cout << "T: "; break;
      }
      
      for(size_t a = 0; a < 8; a++) {
	std::cout << genotype_priors_[i][a] << ", ";
      }
      std::cout << genotype_priors_[i][9] << "]" << std::endl;
    }
  }

  
  Genotype createGameticDNA(ped::Member &m, Base ref) {
    if(ref == N) {
      return NN;
    }
    else {
      return static_cast<Genotype>(floor(genotype_priors_dist_[ref](ran_generator_)));
    }
  }

  Genotype getGermlineTransition(ped::Member &m, Base ref, Genotype mom, Genotype dad) {
    if(ref == N || mom == NN || dad == NN) {
      return NN;
    }
    else {
      // TODO: change germline_dist_ to a 10x10x10 instead of 100x10
      size_t i = findGTIndex(mom, dad);
      return static_cast<Genotype>(floor(germline_dist_[i](ran_generator_)));
    }
  }
  
  Genotype getSomaticTransition(ped::Member &m, Base ref, Genotype gt) {
    if(ref == N || gt == NN) {
      return NN;
    }
    else {
      return static_cast<Genotype>(floor(somatic_dist_[gt](ran_generator_)));
    }
  }

  int call(int n_trials, Base ref, Genotype gt, reads_list &depths) {
    if(ref == N || gt == NN) {
      return 0;
    }

    // Run a Bernoulli trial to determine if model a or model b will be used
    std::array<double, 4> *_read_proportions;
    double _phi = 0.0;
    if(bern_dist_(ran_generator_)) {
      _read_proportions = &read_proportions_a_[ref][gt]; // using parameters from model_a_
      _phi = model_a_.phi;
    }
    else {
      _read_proportions = &read_proportions_b_[ref][gt];
      _phi = model_b_.phi;
    }

    // get a set of random read depths using the given model parameters
    rdm(n_trials, _phi, *_read_proportions, depths);
  }


 protected:
  std::array<double, 10> population_prior(double theta, std::array<double, 4> nuc_freq, std::array<double, 4> prior) {
    double nuc_sum = nuc_freq[0] + nuc_freq[1] + nuc_freq[2] + nuc_freq[3];
    theta = theta/nuc_sum;
    double alpha[4] = {
      theta * nuc_freq[0] + prior[0], theta * nuc_freq[1] + prior[1],
      theta * nuc_freq[2] + prior[2], theta * nuc_freq[3] + prior[3]
    };
    double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
    
    std::array<double, 10> weights;    
    weights[0] = (alpha[0] + alpha[0]*alpha[0]) / alpha_sum / (1.0 + alpha_sum); // AA
    weights[1] = 2*alpha[0]*(alpha[1]) / alpha_sum / (1.0 + alpha_sum); // AC
    weights[2] = 2*alpha[0]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum); // AG
    weights[3] = 2*alpha[0]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum); // AT
    weights[4] = (alpha[1] + alpha[1]*alpha[1]) / alpha_sum / (1.0 + alpha_sum); // CC
    weights[5] = 2*alpha[1]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum); // CG
    weights[6] = 2*alpha[1]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum); // CT
    weights[7] = (alpha[2] + alpha[2]*alpha[2]) / alpha_sum / (1.0 + alpha_sum); // GG
    weights[8] = 2*alpha[2]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum); // GT
    weights[9] = (alpha[3] + alpha[3]*alpha[3]) / alpha_sum / (1.0 + alpha_sum); // TT
    return weights;    
  }


  
  void setPopPriorDist(std::array<double, 10> priors[]) {
    double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    genotype_priors_dist_.emplace_back(std::begin(interval), std::end(interval), priors[0].begin());
    genotype_priors_dist_.emplace_back(std::begin(interval), std::end(interval), priors[1].begin());
    genotype_priors_dist_.emplace_back(std::begin(interval), std::end(interval), priors[2].begin());
    genotype_priors_dist_.emplace_back(std::begin(interval), std::end(interval), priors[3].begin());
  }

  void setGermlineDist(TransitionMatrix &M) {
    double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    for(int i = 0; i < 100; ++i) {
      // TODO: Have been unsuccessfull in passing in an Eigen row into piecewise_constant_dist so
      //  instead had to copy to array. 
      auto r = M.row(i);
      std::array<double, 10> dists;
      for(int a = 0; a < 10; ++a)
	dists[a] = M(i,a);
      
      germline_dist_.emplace_back(std::begin(interval), std::end(interval), dists.begin()); //r.data());
    }
  }

  void setSomaticDist(TransitionMatrix &M) {
    double interval[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    for(int i = 0; i < 10; ++i) {
      auto r = M.row(i);
      std::array<double, 10> dists;
      for(int a = 0; a < 10; ++a)
	dists[a] = M(i,a);
      
      somatic_dist_.emplace_back(std::begin(interval), std::end(interval), dists.begin());//r.data());
    }
  }
  
  void createReadProportions(double epsilon, double omega, read_prop_table &readprops) {
    double u = omega;
    double e = epsilon/3.0;
    double m = 1.0 - 3.0*e;
    double h = 1.0 - 2.0*e;
    
    for(size_t r = 0; r < 4; ++r) {
      for(size_t gt = 0; gt < 10; ++gt){
	readprops[r][gt] = {e, e, e, e};
	Base g1 = gt2bases(gt)[0];
	Base g2 = gt2bases(gt)[1];
	if(g1 == g2) {
	  readprops[r][gt][g1] = m;
	}
	else if(g1 == r) {
	  readprops[r][gt][g1] = h*u/(1.0+u);
	  readprops[r][gt][g2] = h*1.0/(1.0+u);
      }
	else if(g2 == r) {
	  readprops[r][gt][g1] = h*1.0/(1.0+u);
	  readprops[r][gt][g2] = h*u/(1.0+u);
	}
	else {
	  readprops[r][gt][g1] = h/2.0;
	  readprops[r][gt][g2] = h/2.0;
	}      
      }
    }
  }

  
  void initGSL() {
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
  }

  void initBernoulliGenerator(double p) {
    bern_dist_ = std::bernoulli_distribution(p); 
    
  }


  void rdm(int m, double phi, std::array<double, 4> &params, reads_list &depths) {

    //std::cout << "phi = " << phi << std::endl;
    double phi_a = 1.0/(phi + 1.0);
    phi_a = (1.0 - phi_a)/phi_a;
    std::array<double, 4> alphas;
    for(int a = 0; a < 4; ++a) {
      alphas[a] = phi_a*params[a];    
    }
    
    // Run a single trial call
    unsigned int init_vals[4];
    gsl_ran_multinomial(rng, 4, 1, params.begin(), init_vals);
    
    double dirichlet_params[4];
    for(int a = 0; a < 4; ++a) {
      dirichlet_params[a] = alphas[a]+init_vals[a];
    }

    /*
    std::cout << "dirichlet_params: ";
    std::cout << "[" << dirichlet_params[0];
    for(int a = 1; a < 4; a++)
      std::cout << ", " << dirichlet_params[a];
    std::cout << "]" << std::endl;
    */
    
    double thetas[4];
    gsl_ran_dirichlet(rng, 4, dirichlet_params, thetas);
    
    /*
      std::cout << "thetas:" << std::endl;
      std::cout << "[" << thetas[0];
      for(int a = 1; a < 4; a++)
      std::cout << ", " << thetas[a];
      std::cout << "]" << std::endl;
    */
    
    
    //std::array<unsigned int, 4> calls;
    //gsl_ran_multinomial(r, 4, m-1, thetas, calls.begin());
    gsl_ran_multinomial(rng, 4, m-1, thetas, depths.begin());
    for(int a = 0; a < 4; ++a) {
      depths[a] += init_vals[a];
  }


  /*
  //std::cout << "thetas:" << std::endl;
  std::cout << "[" << counts[0];
  for(int a = 1; a < 4; a++)
    std::cout << ", " << counts[a];
  std::cout << "]" << std::endl;
  */
  
}
  
  
  
 protected:
  double ref_weight_; // weight given to reference
  double theta_; // population diversity
  std::array<double, 4> nuc_freqs_; // nucleotide freq's in A,C,G,T order
  double mu_;
  double mu_soma_;
  gamma_t model_a_;
  gamma_t model_b_;
  
  
  // Genotype priors
  std::array<double, 10> genotype_priors_[4];
  std::vector<std::piecewise_constant_distribution<>> genotype_priors_dist_; // used for random generator

  // Distribution for germline transitions
  TransitionMatrix germline_trans_matrix_;
  std::vector<std::piecewise_constant_distribution<>> germline_dist_;

  // Distribution for somatic transitions
  TransitionMatrix somatic_trans_matrix_;
  std::vector<std::piecewise_constant_distribution<>> somatic_dist_;

  // Table for problabilities of the sequencer reading any NT for any given ref+genotype
  read_prop_table read_proportions_a_; // using parameters from model_a_
  read_prop_table read_proportions_b_; // using parameters from model_b_
  std::bernoulli_distribution bern_dist_; // A bernoulli distribution for choosing between models a and b.
  
  // random number generators
  //   TODO: allow users to pass in their own generators to help with parallelization
  std::mt19937 ran_generator_;
  const gsl_rng *rng;
};


} // namespace model
} // namespace sim


#endif
