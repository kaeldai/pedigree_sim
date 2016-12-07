#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <boost/algorithm/string.hpp>


namespace sim {
namespace ped {


struct Member {
  std::string id;
  int gender;
  bool is_founder;
  int mpos;
  int dpos;
  bool dna_set = false;  
  
  struct stats_t {
    size_t n_unknowns = 0;

    // Statistics for founder
    size_t n_hom_ref_match = 0; // Homozygous, matches ref
    size_t n_het_ref_match = 0; // Heterozygous, one of the alleles matches ref
    size_t n_hom_snp = 0; // Homozgous, not matching ref
    size_t n_het_snp = 0; // Heterozygous, not matching ref

    // Germline dna statistics
    size_t n_mendelians = 0; // GL transition matches medelian inheritence (up to phasing info)
    size_t n_single_mut = 0; // One site mutated
    size_t n_double_mut = 0; // Both alleles mutated

    // Somatic dna statistics
    size_t n_som_matches = 0; // Somatic dna matches 
    size_t n_som_single_mut = 0; // Single mutation 
    size_t n_som_double_mut = 0; // double mutation
    
  } stats;
  
  Member(std::string id, int gender, bool is_founder, int mpos, int dpos) :
      id(id), gender(gender), is_founder(is_founder), mpos(mpos), dpos(dpos) { }
  
};

std::string family_name;
std::vector<Member> family;

#define get_dad(c) family[c.dpos]
#define get_mom(c) family[c.mpos]


void parse_pedigree(std::string &file_name) {
  std::ifstream ped_file(file_name);
  if(!ped_file.is_open()) {
    throw std::runtime_error("Unable to open pedigree file " + file_name + ". Exiting.");
  }

  std::vector<std::array<std::string, 6>> ped_table;
  std::string line;
  while(std::getline(ped_file, line)) {
    std::vector<std::string> cols;
    boost::split(cols, line, boost::is_any_of("\t"));

    if(cols.size() < 5) {
      throw std::runtime_error("All rows in pedigree file must have atleast five columns.");
    }
    
    std::string is_founder = "0";
    if(cols[2] == "0" && cols[3] == "0") {
      is_founder = "1";
    }
    else if(cols[2] == "0" || cols[3] == "0") {
       throw std::runtime_error("Cannot handle partial parentage for " + cols[1] + ".");
    }

    ped_table.emplace_back(std::array<std::string, 6>{cols[0], cols[1], cols[2], cols[3], cols[4], is_founder});   
  }

  ped_file.close();

  // may want to use the family id later
  family_name = ped_table[0][0];
  
  // Check that all members in pedgiree are of the same family
  for(size_t a = 1; a < ped_table.size(); ++a) {
    if(ped_table[a][0] != ped_table[a-1][0]) {
      throw std::runtime_error("Pedigree has multiple families.");
    }
  }

  // First iteration, create family vector
  for(std::array<std::string, 6> row : ped_table) {
    bool is_founder = (row[5] == "1");
    int gender = std::stoi(row[4]);
    if(is_founder) {
      // put founders in the front of the queue
      family.emplace(family.begin(), row[1], gender, is_founder, -1, -1);
    }
    else {
      family.emplace_back(row[1], gender, is_founder, -1, -1);
    }
  }

  // create a map between members and their positions in the family
  std::map<std::string, int> positions;
  for(int pos = 0;  pos < ped_table.size(); ++pos) {
    std::string id = ped_table[pos][1];
    std::map<std::string, int>::iterator it = positions.find(id);
    if(it != positions.end()) {
      throw std::runtime_error("Member " + id + " appears more than once in pedigree.");
    }
    positions.insert(std::make_pair(id, pos));
  }

  // for each child set the positions of their mom/dad
  for(std::array<std::string, 6> row : ped_table) {
    if(row[5] == "1")
      continue;

    std::map<std::string, int>::iterator it = positions.find(row[1]);
    int cpos = it->second;
    it = positions.find(row[2]);
    if(it == positions.end()) {
      throw std::runtime_error("Father for " + row[1] + " missing in pedigree.");
    }
    else{
      family[cpos].dpos = it->second;
    }
      
    it = positions.find(row[3]);
    if(it == positions.end()) {
      throw std::runtime_error("Mother for " + row[1] + " missing in pedigree.");
    }
    else{
      family[cpos].mpos = it->second;
    }    
  }

  /*
  for(Member m : family) {
    std::cout << m.id;
    if(!m.is_founder) {
      std::cout << " (" << get_dad(m).id << ", " << get_mom(m).id << ")";
    }
    std::cout << std::endl;
  }
  */
  
}

void collect_founder_stats(Member &m, Base ref, Genotype gt) {
  if(ref == N) {
    ++m.stats.n_unknowns;
  }
  else {
    bool is_hom = (gt == AA || gt == CC || gt == GG || gt == TT);
    bool contains_ref = (gt_contains(ref, gt) == 1);
    if(contains_ref && is_hom)
      ++m.stats.n_hom_ref_match;
    else if(contains_ref && !is_hom)
      ++m.stats.n_het_ref_match;
    else if(!contains_ref && is_hom)
      ++m.stats.n_hom_snp;
    else
      ++m.stats.n_het_snp;
  }
}

void collect_germline_stats(Member &mem, Genotype cgt, Genotype mgt, Genotype dgt) {
  if(cgt == NN || mgt == NN || dgt == NN) {
    return;
  }  

  // Create a table of all possible mendelian combinations using parents genotypes
  //   TODO: preprocess possible combinations
  Base *m = split_table[mgt];
  Base *d = split_table[dgt];
  std::array<Genotype, 4> mendelian_gt = {
    bases2gt(m[0],d[0]), bases2gt(m[0],d[1]),
    bases2gt(m[1],d[0]), bases2gt(m[1],d[1])
  };

  // Check if child genotype matches with parents
  for(int a = 0; a < 4; a++) {
    if(cgt == mendelian_gt[a]) {
      ++mem.stats.n_mendelians;
      return;
    }
  }

  // Check to see if there was a mutation in one or both of parents during meiosis
  Base *c = split_table[cgt];
  if(c[0] == m[0] || c[0] == m[1] || c[1] == m[0] || c[1] == m[1] ||
     c[0] == d[0] || c[0] == d[1] || c[1] == d[0] || c[1] == d[1]) {
    ++mem.stats.n_single_mut;
  }
  else {
    ++mem.stats.n_double_mut;
  }
}

void collect_somatic_stats(Member &mem, Genotype germ_gt, Genotype soma_gt) {
  if(germ_gt == NN || soma_gt == NN) {
    return;
  }

  if(germ_gt == soma_gt) {
    ++mem.stats.n_som_matches;
  }
  else {
    // Check to see if mitosis had mutations in either one or both alleles
    Base *s = split_table[soma_gt];
    if((gt_contains(s[0], germ_gt) == 1) || (gt_contains(s[1], germ_gt) == 1)) {
      ++mem.stats.n_som_single_mut;
    }
    else {
      ++mem.stats.n_som_double_mut; // very rare
    }     
  } 
}


} // namespace ped
} // namespace sim

#endif
