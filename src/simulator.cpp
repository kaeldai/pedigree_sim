#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <random>

#include "simulator.h"
#include "output.h"
#include "pedigree.h"
#include "models.h"

#include <htslib/faidx.h>
#include <htslib/hts.h>


using namespace sim;

int main(int argc, char *argv[]) {

  // read command line options
  if(parse_options(argc, argv) == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }
 
  // read the pedigree file
  ped::parse_pedigree(arg.ped);
  std::vector<ped::Member> family = ped::family;

  std::cout << "> Reading reference." << std::endl;
    
  // get information about contig
  size_t pos_begin = 0;
  size_t pos_end = 0;
  std::string contig_name;
  parse_region(arg.region, contig_name, pos_begin, pos_end);

  std::string base_fname = arg.prefix;
  if(base_fname.empty()) {
    base_fname = ped::family_name;
  }
  
  // read in reference contig
  int contig_len = 0;
  const faidx_t *contig_ref = fai_load(arg.input.c_str());
  char *contig_str = fai_fetch(contig_ref, arg.region.c_str(), &contig_len);

  // Create a reference
  std::vector<Base> reference(contig_len);  
  for(size_t site = 0; site < contig_len; ++site) {
    reference[site] = char2base(contig_str[site]);
  }

  std::vector<std::vector<Genotype>> gametic_dna(family.size());
  std::vector<std::vector<Genotype>> somatic_dna(family.size());
	     
  model::DNGModel model(arg.ref_weight, arg.theta, arg.nuc_freqs, arg.mu, arg.mu_somatic, arg.model_a, arg.model_b);


  std::cout << "> Creating germline DNA for founders." << std::endl;
  
  // Use the population priors to intialize founder dna
  size_t remaining = family.size();
  for(int mem = 0; mem < family.size() && family[mem].is_founder; ++mem) {
    ped::Member &m = family[mem]; // reference to member
    std::vector<Genotype> &mdna = gametic_dna[mem]; // reference to member's dna
    for(size_t site = 0; site < contig_len; site++) {
      Base ref = reference[site];
      Genotype gt = model.createGameticDNA(m, ref);
      mdna.push_back(gt);
      collect_founder_stats(m, ref, gt);
    }
    m.dna_set = true;
    --remaining;
  }


  std::cout << "> Creating germline DNA for children." << std::endl;
  
  // Create gametic DNA for non-founders
  while(remaining) {
    for(size_t mem; mem < family.size(); ++mem) {
      ped::Member &m = family[mem];
      if(m.dna_set) {
	// member already has gametic DNA
	continue;
      }
      ped::Member &dad = get_dad(m);
      ped::Member &mom = get_mom(m);
      if(!dad.dna_set || !mom.dna_set) {
	// one of member's parent's doesn't have their dna, come back later
	continue;
      }

      // member doesn't have DNA, but both parents do. use a germline tranformation
      // to create new DNA for member.
      std::vector<Genotype> &cdna = gametic_dna[mem]; // child's dna
      std::vector<Genotype> &ddna = gametic_dna[m.dpos]; // dad's dna
      std::vector<Genotype> &mdna = gametic_dna[m.mpos]; // mom's dna
      for(size_t site = 0; site < contig_len; ++site) {
	Base ref = reference[site];
	Genotype dgt = ddna[site];
	Genotype mgt = mdna[site];
	Genotype cgt = model.getGermlineTransition(m, ref, mgt, dgt);
	cdna.push_back(cgt);
	collect_germline_stats(m, cgt, mgt, dgt);
      }
      m.dna_set = true;
      --remaining;
    }
  }

  std::cout << "> Creating somatic DNA." << std::endl;
  
  // Create somatic mutation
  for(size_t mem = 0; mem < family.size(); ++mem) {
    ped::Member &m = family[mem];
    std::vector<Genotype> &dna = gametic_dna[mem];

    for(size_t site = 0; site < contig_len; ++site) {
      Base ref = reference[site];
      Genotype gametic_gt = dna[site];
      Genotype somatic_gt = model.getSomaticTransition(m, ref, gametic_gt);
      somatic_dna[mem].push_back(somatic_gt);
      collect_somatic_stats(m, gametic_gt, somatic_gt);
    }
  }

  
  // Save germline and somatic genotypes in a vcf file.
  io::VCFOutput vcf_file((base_fname + ".vcf").c_str());
  vcf_file.addContig(contig_name.c_str(), contig_len);
  for(ped::Member &m : family) {
   vcf_file.addSample((std::string("GL-") + m.id).c_str());
   vcf_file.addSample((std::string("SL-") + m.id).c_str());
  }
  vcf_file.writeHeader();
  std::vector<Genotype> gts(family.size()*2);
  for(size_t site = 0; site < contig_len; ++site) {
    Base ref = reference[site];
    vcf_file.newSite(contig_name, (pos_begin + site), ref);
    for(size_t mem = 0; mem < family.size(); ++mem) {
      gts[2*mem] = gametic_dna[mem][site];
      gts[2*mem+1] = somatic_dna[mem][site];
    }
    vcf_file.addGenotypes(gts);
    vcf_file.writeSite();
  }
  vcf_file.close();
  
  
  std::cout << "> Creating read calls." << std::endl;
  
  // Create library reads and output to file file
  io::TadOutput tad_file((base_fname + ".tad").c_str());
  tad_file.addContig(contig_name, contig_len);
  std::vector<std::string> lb_names;
  for(ped::Member &m : family) {
    lb_names.emplace_back(m.id);
    tad_file.addLibrary(m.id);
  }
  tad_file.writeHeader();

  
  std::vector<reads_list> read_depths(family.size());
  for(size_t site = 0; site < contig_len; ++site) {
    Base ref = reference[site];
    size_t pos = pos_begin + site + 1;
    tad_file.newSite(contig_name, pos, ref);
    for(size_t mem = 0; mem < family.size(); ++mem) {
      model.call(arg.read_depths, ref, somatic_dna[mem][site], read_depths[mem]); 
    }
    tad_file.addCalls(read_depths); 
    tad_file.writeSite();
  }
  tad_file.close();

  io::CSVOutput statsfile((base_fname + ".csv").c_str());
  statsfile.addStats(family);
  statsfile.close();

}
