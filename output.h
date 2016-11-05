#ifndef OUTPUT_H
#define OUTPUT_H

#include "simulator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace sim {
namespace io {


unsigned int allele_orders[4][3] = {{1, 2, 3}, {2, 3, 0}, {3, 0, 1}, {0, 1, 2}};

class VariantOutput {
 protected:
  std::vector<std::string> contigs_;
  std::vector<std::string> libraries_;

  Base c_ref_;
  std::string c_contig_;
  size_t c_pos_;
  std::vector<reads_list> c_depths_;

  std::ofstream outfile_;
  
  
 public:
  VariantOutput(std::string &filename) {
    outfile_.open(filename.c_str());
  }

  VariantOutput(const char *filename) {
    outfile_.open(filename);
  }

  void addContig(std::string &contig_name) {

  }
  
  void addLibrary(std::string &lb_name) {
    libraries_.emplace_back(lb_name);    
    
  }

  void writeHeader() {

  }

  void newSite(std::string &contig_name, size_t pos, Base ref) {
    c_contig_ = contig_name;
    c_pos_ = pos;
    c_ref_ = ref;
  }

  void addCalls(std::vector<reads_list> &depths) {
    c_depths_ = depths;
  }
  
 
  void writeSite() {
    // Get references
    std::string alleles = base2str(c_ref_);
    std::vector<std::string> depths_str(libraries_.size());
    for(int a = 0; a < c_depths_.size(); ++a) {
      depths_str[a] = std::to_string(c_depths_[a][c_ref_]);
    }

    for(Base allele : allele_orders[c_ref_]) {
      int total = 0;
      for(int lb = 0; lb < c_depths_.size(); ++lb) {
	total += c_depths_[lb][allele];
      }
      if(total > 0) {
	alleles += base2str(allele);
	for(int lb = 0; lb < c_depths_.size(); ++lb) {
	  depths_str[lb] += "," + std::to_string(c_depths_[lb][allele]);
	}
      }      
    }
    
    outfile_ << c_contig_ << "\t" << c_pos_ << "\t" << alleles;
    for(std::string &dp : depths_str) {
      outfile_ << "\t" << dp;
    }
    outfile_ << std::endl;
      
  }

  void close() {
    outfile_.close();
  }
  
};


} // namespace io
} // namespace sim

#endif
